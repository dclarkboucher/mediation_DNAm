#Script for computing the percent relative bias in estimating the global
#mediation effect

library(tidyverse)
library(bama)
ndat <- 100

settings <- c("bl","hc","ns","ln","ln_hc","ln_ns")
datasets <- c(1,2,3,4)
seed2 = 1:ndat

tab <- expand_grid(setting = settings, dataset = datasets,
                   seed2 = 1:ndat)

mylist <- list()


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

for(setting in settings){
  
  mylist[[setting]] <- list()
  
  for(dataset in datasets){
    
    mylist[[setting]][[dataset]] <- list()
    
    for(seed in seed2){
      
      results <- list()
      
      # Load dataset
      load(paste0("simulation_datasets/sim_data_",setting,"_d",dataset,"_s",seed,".rda"))
      
      # Make empty dataset of results
      dat <- tibble(mediator = paste0("M",1:2000), beta_m, alpha_a,
                    ab = alpha_a * beta_m)
      
      # Determine which mediators are active (only needed for PLASSO)
      if (grepl("ns",setting)) {
        truth <- as.numeric((1:ncol(M)) %in% which_ab)
        
      } else{
        truth <- as.numeric(alpha_a * beta_m != 0)
      }
      
      # Load PLASSO results
      load(paste0("simulation_results/plasso/RUN_",setting,
                  "_d",dataset,"_s",seed,".RData"))
      
      #determine lambda for fdr <- 0.10
      fdps_by_lambda <- apply(AB.est,2,plasso_fdp,want_cols,truth)
      which_lambda_use <- which(fdps_by_lambda < 0.1)[1]
      
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
          mediator = str_replace(mediator,"V","M"),
          beta_hat = beta,
          alpha_hat = alpha,
          ab_hat = alpha_hat * beta_hat,
          pv = ab_pv
        ) |> 
        left_join(x = dat, by = "mediator") |> 
        mutate(
          across(c(beta_hat,alpha_hat, ab_hat), ~ ifelse(is.na(.x),0,.x)),
          pv = ifelse(is.na(pv),1,pv),
          method = "HDMA"
        )
      
      
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
          mediator = str_replace(mediator,"V","M"),
          beta_hat = beta,
          alpha_hat = alpha,
          ab_hat = alpha_hat * beta_hat,
          pv = ab_pv
        ) |> 
        left_join(x = dat, by = "mediator") |> 
        mutate(
          across(c(beta_hat,alpha_hat, ab_hat), ~ ifelse(is.na(.x),0,.x)),
          pv = ifelse(is.na(pv),1,pv),
          method = "HIMA"
        )
      
      
      #Load BSLMM results
      bslmm <- read_csv(paste0("simulation_results/bslmm/sim_out_bslmm_",setting,"_d",
                               dataset,"_s",seed,".csv"))[,-1]
      results$bslmm <-
        with(
          bslmm,
          dat |> 
            mutate(
              beta_hat = beta_hat,
              alpha_hat = alpha_hat,
              ab_hat = tie_hat[1] / n(), #this'll work fine.
              pip = pip,
              method = "BSLMM"
            )
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
        mutate(
          across(c(beta_hat,alpha_hat, ab_hat), ~ ifelse(is.na(.x),0,.x)),
          pv = ifelse(is.na(pv),1,pv),
          method = "MedFix"
        )
      
      # Aggregate results for group 1 methods
      d1 <-
        bind_rows(results) |>
        mutate(dataset = paste0("Setting ",dataset)) |>
        select(dataset,method,mediator,everything()) |> 
        group_by(method) |> 
        summarize(dataset = dataset[1], 
                  tie = sum(ab),
                  tie_hat = sum(ab_hat),
                  tie_dif = tie_hat - tie,
                  tie_dif1 = abs(tie_dif),
                  tie_mse = tie_dif^2,
                  tie_prb = tie_dif1 / tie * 100,
                  setting = setting,
                  dataset = dataset,
                  seed2 = seed) |> 
        ungroup()
      
      # PCMA
      load(file = paste0("simulation_results/pcma/sim_out_pcma_",setting,"_d",
                         dataset,"_s",seed,".rda"))
      d2 <-
        tibble(method = "PCMA",
               tie = sum(alpha_a * beta_m),
               tie_hat = out_pcma$IE.total[1],
               tie_dif = tie_hat - tie,
               tie_dif1 = abs(tie_dif),
               tie_mse = tie_dif^2,
               tie_prb = tie_dif1 / tie * 100,
               setting = setting,
               dataset = paste0("Setting ",dataset),
               seed2 = seed)
      
      # HILMA
      load(file = paste0("simulation_results/hilma/sim_out_hilma_",
                         setting,"_d",dataset,"_s",seed,".rda"))
      d3 <-
        tibble(method = "HILMA",
               tie = sum(alpha_a * beta_m),
               tie_hat = out_hilma$beta_hat[1,1] * sd(Y) / sd(A),
               tie_dif = tie_hat - tie,
               tie_dif1 = abs(tie_dif),
               tie_mse = tie_dif^2,
               tie_prb = tie_dif1 / tie * 100,
               setting = setting,
               dataset = paste0("Setting ",dataset),
               seed2 = seed)
      
      # Merge all results
      mylist[[setting]][[dataset]][[seed]] <- bind_rows(list(d1,d2,d3))
    }
  }
}

results <- 
  mylist |> 
  map(bind_rows) |> 
  bind_rows()

results |> 
  group_by(setting,dataset,method) |> 
  summarize(prb = mean(tie_prb), prb_l = quantile(tie_prb,0.025),
            prb_u = quantile(tie_prb,0.975)) |>
  ungroup() |> 
  write_csv("simulation_results/percent_rel_bias.csv")


