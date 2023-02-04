# Script to process global mediation effects
library(tidyverse)
library(bama)

rm(list = ls())

ndat <- 1 # Again, set this 100 if replicating the full study

# All settings to loop through
settings <- c("bl", "hc", "ns","ln","ln_hc","ln_ns")
datasets <- c(1, 2, 3, 4)
seed2 = 1:ndat

# Empty lists of results
mylist_tpr <- list()
mylist_mse0 <- list()
mylist_mse1 <- list()

# Function determine p-value cutoff for empirical FDR correction
test_fdr <- 
  function(score, truth = as.numeric(alpha_a * beta_m != 0), fdr = 0.1, cutoff.lower = T) {
    cutoffs <- sort(unique(score))
    
    
    
    if (cutoff.lower) {
      for (c in rev(cutoffs)) {
        sig <- score <= c
        
        fdp <- sum(sig & !truth) / sum(sig)
        
        if (fdp <= fdr) {
          break
        }
        
      }
      
    } else{
      for (c in cutoffs) {
        sig <- score >= c
        
        fdp <- sum(sig & !truth) / sum(sig)
        
        if (fdp <= fdr) {
          break
        }
        
      }
    }
    
    
    if(fdp <= fdr){
      return(sig)
      
    }else{ #this is for when the loop expires without finding a solution
      return(rep(F,length(sig)))
    }
    
    
    
  }

# Function to determine false discovery proportion for a given
# set of pathway LASSO estimates, since we have estimates for 
# every lambda. We choose the lowest lambda (i.e., the least
# shrunken model) attaining an FDP < 0.10
plasso_fdp <- function(estimates, which_cols, truth){
  # estimates: contribution estimates for a subset of mediators, some zero
  # which_cols: which_entries in "truth" were estimated at all 
  # truth: whether each of the 2000 mediators is active
  selected <- truth & F
  selected[which_cols] <- estimates != 0
  if(!any(selected)){return(1)} # Doesn't really matter what you return here.
  fdp <- sum(selected & !truth) / sum(selected)
  return(fdp)
  
}

for (setting in settings) {
  mylist_tpr[[setting]] <- list()
  mylist_mse0[[setting]] <- list()
  mylist_mse1[[setting]] <- list()
  
  for (dataset in datasets) {
    mylist_tpr[[setting]][[dataset]] <- list()
    mylist_mse0[[setting]][[dataset]] <- list()
    mylist_mse1[[setting]][[dataset]] <- list()
    
    for (seed in seed2) {
      temp <- ls()
      
      results <- list()
      
      # Load dataset
      load(paste0("simulation_datasets/sim_data_",setting,"_d",dataset,"_s",seed,".rda"))
      
      # Make empty dataset of results
      dat <- tibble(mediator = paste0("M", 1:2000),
                    beta_m,
                    alpha_a,
                    ab = alpha_a * beta_m)
      
      # Determine which mediators are active
      if (grepl("ns",setting)) {
        truth <- as.numeric((1:ncol(M)) %in% which_ab)
        
      } else{
        truth <- as.numeric(alpha_a * beta_m != 0)
      }
      
      # Load one-at-a-time results
      load(
        file = paste0(
          "simulation_results/one-at-a-time/sim_out_uni_",
          setting,
          "_d",
          dataset,
          "_s",
          seed,
          ".rda"
        )
      )
      
      results$uni <-
        out_uni |>
        rename(pv = ab_pv) |> 
        select(mediator,beta_hat,alpha_hat,ab_hat, pv) |> 
        left_join(x = dat, by = "mediator") |>
        mutate(method = "Univariate")
      
      # Load PLASSO results
      load(paste0("simulation_results/plasso/RUN_",setting,
                  "_d",dataset,"_s",seed,".RData"))
      
      #determine lambda for fdr <- 0.10
      fdps_by_lambda <- apply(AB.est,2,plasso_fdp,want_cols,truth)
      which_lambda_use <- which(fdps_by_lambda < 0.1)[1]
        # if this is NA, that means nothing was chosen. that is fine. 
      results$plasso <-
        tibble(
          mediator = paste0("M",want_cols),
          beta_hat = B.est[,which_lambda_use],
          alpha_hat = A.est[,which_lambda_use],
          ab_hat = AB.est[,which_lambda_use],
          pv = ifelse(ab_hat == 0, 1, 0.05)
        ) |> 
        left_join(x = dat, by = "mediator") |>
        mutate(across(c(beta_hat, alpha_hat, ab_hat), ~ ifelse(is.na(.x), 0, .x)),
               pv = ifelse(is.na(pv), 1, pv),
               method = "P-LASSO")
    
      # Load HDMA results
      load(
        file = paste0(
          "simulation_results/hdma/sim_out_hdma_",
          setting,
          "_d",
          dataset,
          "_s",
          seed,
          ".rda"
        )
      )
      results$hdma <-
        out_hdma |>
        transmute(
          mediator = str_replace(mediator, "V", "M"),
          beta_hat = beta,
          alpha_hat = alpha,
          ab_hat = alpha_hat * beta_hat,
          pv = ab_pv
        ) |>
        left_join(x = dat, by = "mediator") |>
        mutate(across(c(beta_hat, alpha_hat, ab_hat), ~ ifelse(is.na(.x), 0, .x)),
               pv = ifelse(is.na(pv), 1, pv),
               method = "HDMA")
      
      #Load HIMA results
      load(
        file = paste0(
          "simulation_results/hima/sim_out_hima_",
          setting,
          "_d",
          dataset,
          "_s",
          seed,
          ".rda"
        )
      )
      
      if(is.null(out_hima)){
        out_hima <- 
          data.frame(
            mediator = "V1",
            beta = 0,
            alpha = 0,
            alpha_beta = 0,
            te = 0,
            ab_pv = 1
          )
      }
      
      results$hima <-
        out_hima |>
        transmute(
          mediator = str_replace(mediator, "V", "M"),
          beta_hat = beta,
          alpha_hat = alpha,
          ab_hat = alpha_hat * beta_hat,
          pv = ab_pv
        ) |>
        left_join(x = dat, by = "mediator") |>
        mutate(across(c(beta_hat, alpha_hat, ab_hat), ~ ifelse(is.na(.x), 0, .x)),
               pv = ifelse(is.na(pv), 1, pv),
               method = "HIMA")
      
      
      # Load BSLMM results
      bslmm <- read_csv(paste0("simulation_results/bslmm/sim_out_bslmm_",setting,"_d",
                               dataset,"_s",seed,".csv"))[,-1]
      results$bslmm <-
        bslmm |> 
        select(beta_hat, alpha_hat, pip) |> 
        mutate(
          ab_hat = alpha_hat * beta_hat,
          method = "BSLMM"
        )

      # Load MedFix results
      load(
        file = paste0(
          "simulation_results/med/sim_out_med_",
          setting,
          "_d",
          dataset,
          "_s",
          seed,
          ".rda"
        )
      )
      results$med <-
        out_med |>
        transmute(
          mediator = variable,
          beta_hat = beta,
          alpha_hat = alpha,
          ab_hat = alpha_hat * beta_hat,
          pv = pv
        ) |>
        left_join(x = dat, by = "mediator") |>
        mutate(across(c(beta_hat, alpha_hat, ab_hat), ~ ifelse(is.na(.x), 0, .x)),
               pv = ifelse(is.na(pv), 1, pv),
               method = "MedFix")
      
      d.mse0 <-
        results |>
        bind_rows() |>
        mutate(error = (ab_hat - ab) ^ 2, #error for MSE
               prb = 100 * abs(ab_hat - ab) / ab, #percent relative bias
               method = factor(
                 method,
                 levels = c('Univariate', 'HIMA', 'HDMA', 'MedFix', 'P-LASSO',
                            'BSLMM')
               )) |>
        group_by(method) |>
        filter(!truth) |>
        summarize(
          mse = mean(error)
        ) |>
        ungroup()  |> 
        mutate(
          rmse = mse / mse[method == "Univariate"],
          setting = setting,
          dataset = paste0("Setting ",dataset),
          seed2 = seed
        )
      
      d.mse1 <-
        results |>
        bind_rows() |>
        mutate(error = (ab_hat - ab) ^ 2, #error for MSE
               prb = 100 * abs(ab_hat - ab) / ab, #percent relative bias
               method = factor(
                 method,
                 levels = c('Univariate', 'HIMA', 'HDMA', 'MedFix','P-LASSO',
                            'BSLMM')
               )) |>
        group_by(method) |>
        filter(as.logical(truth)) |>
        summarize(
          mse = mean(error),
        ) |>
        ungroup()  |> 
        mutate(
          rmse = mse / mse[method == "Univariate"],
          setting = setting,
          dataset = paste0("Setting ",dataset),
          seed2 = seed
        )
      
      d.tpr <-
        results |>
        bind_rows() |>
        mutate(score = ifelse(is.na(pv), 1 - pip, pv),
               method = factor(
                 method,
                 levels = c('Univariate', 'HIMA', 'HDMA', 'MedFix','P-LASSO',
                            'BSLMM')
               )) |>
        group_by(method) |>
        mutate(sig = test_fdr(score, truth)) |>
        summarize(
          np = sum(sig),
          ntp = sum(sig & truth),
          nt = sum(truth),
          tpr = sum(sig & truth) / sum(truth),
          nfp = sum(sig & !truth),
          fpr = nfp / sum(!truth)
        ) |>
        ungroup() |> 
        mutate(setting = setting,
               dataset = paste0("Setting ",dataset),
               seed2 = seed)
      
      
      mylist_tpr[[setting]][[dataset]][[seed]] <- d.tpr
      mylist_mse0[[setting]][[dataset]][[seed]] <- d.mse0
      mylist_mse1[[setting]][[dataset]][[seed]] <- d.mse1
      
      #remove temporary objects (but is everything not temporary?)
      rm(list = setdiff(ls(), temp))
      
    }
  }
}

# True positive rate
results |>
  rename(tpr1 = tpr) |> 
  group_by(setting, dataset, method) |>
  summarize(
    tpr = mean(tpr1),
    tpr_l = quantile(tpr1, 0.025),
    tpr_u = quantile(tpr1, 0.975)
  ) |>
  as.data.frame() |> 
  write_csv("simulation_results/tpr.csv")


# False positive rate - looked at this but not in the paper. They were low.
results |>
  rename(fpr1 = fpr) |> 
  group_by(setting, dataset, method) |>
  summarize(
    fpr = mean(fpr1),
    fpr_l = quantile(fpr1, 0.025),
    fpr_u = quantile(fpr1, 0.975)
  ) |>
  as.data.frame() |> 
  write_csv("simulation_results/fpr.csv")


# MSE for inactive mediators
mylist_mse0 |>
  map(bind_rows) |>
  bind_rows() |>
  rename(mse0 = mse, rmse0 = rmse) |> 
  group_by(setting,dataset,method) |> 
  summarize(
    mse = mean(mse0),
    mse_l = quantile(mse0, 0.025),
    mse_u = quantile(mse0, 0.975),
    rmse = mean(rmse0),
    rmse_l = quantile(rmse0, 0.025),
    rmse_u = quantile(rmse0, 0.975)
  ) |>
  as.data.frame() |> 
  write_csv("simulation_results/mse_inactive.csv")


# MSE among active mediators
mylist_mse1 |>
  map(bind_rows) |>
  bind_rows() |>
  rename(mse0 = mse,rmse0 = rmse) |> 
  group_by(setting,dataset,method) |> 
  summarize(
    mse = mean(mse0),
    mse_l = quantile(mse0, 0.025),
    mse_u = quantile(mse0, 0.975),
    rmse = mean(rmse0),
    rmse_l = quantile(rmse0, 0.025),
    rmse_u = quantile(rmse0, 0.975)
  ) |>
  as.data.frame() |> 
  write_csv("simulation_results/mse_active.csv")






