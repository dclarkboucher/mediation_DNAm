library(tidyverse)
source("utils/get_data_confounded_effects.R")

out_loc <- "supplementary_analyses/confounded_effects_results/"


filename <- function(method, seed, vu){
  
  if (method %in% c("hilma", "pcma", "plasso", "pmed")){
    app <- ".rda"
  } else{
    app <- ".rds"
  }
  
  paste0(out_loc, method,"_s", seed,"_v",vu, app)
  
}

tab <- expand_grid(seed = seq_len(100),
                   vu = c(1, 2, 3)
)


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



results <- list()
for(r in 1:nrow(tab)){
  
  seed <- tab$seed[r]
  vu <- tab$vu[r]
  
  # Produce data
  dat <-
    get_data_umc(
      seed2 = seed,
      vu = vu
    )
  
  #make empty dataset
  d <- 
    with(dat,
         tibble(
           mediator = paste0("M", 1:2000),
           ab = alpha_a * beta_m,
           true = ab != 0
         )
    )
  
  # Process univariate
  out_uni <- read_rds(filename("uni", seed, vu))
  uni <-
    out_uni |>
    rename(pv = ab_pv) |> 
    select(mediator, ab_hat, pv) |> 
    left_join(x = d, by = "mediator") |>
    mutate(method = "Univariate",
           tie_hat = sum(ab_hat))
  rm(out_uni)
  
  
  # Process HIMA
  out_hima <- read_rds(filename("hima",seed,vu))
  hima <-
    out_hima |>
    rownames_to_column("mediator") |>
    transmute(
      mediator = str_replace(mediator, "V", "M"),
      ab_hat = alpha * beta,
      pv = max.p
    ) |>
    left_join(x = d, by = "mediator") |>
    mutate(across(c(ab_hat), ~ ifelse(is.na(.x), 0, .x)),
           pv = ifelse(is.na(pv), 1, pv),
           method = "HIMA",
           tie_hat = sum(ab_hat))
  rm(out_hima)
  
    
  # Process HDMA
  out_hdma <- read_rds(filename("hdma", seed,vu))
  hdma <-
    out_hdma |>
    rownames_to_column("mediator") |>
    transmute(
      mediator = str_replace(mediator, "V", "M"),
      ab_hat = alpha * beta,
      pv = P.value
    ) |>
    left_join(x = d, by = "mediator") |>
    mutate(across(c(ab_hat), ~ ifelse(is.na(.x), 0, .x)),
           pv = ifelse(is.na(pv), 1, pv),
           method = "HDMA",
           tie_hat = sum(ab_hat))
  rm(out_hdma)
    
    # Process MedFix
    out_med <- read_rds(filename("med", seed,vu))
    med <-
      out_med |>
      transmute(
        mediator = variable,
        ab_hat = alpha * beta,
        pv = pv,
      ) |>
      left_join(x = d, by = "mediator") |>
      mutate(across(c(ab_hat), ~ ifelse(is.na(.x), 0, .x)),
             pv = ifelse(is.na(pv), 1, pv),
             method = "MedFix",
             tie_hat = sum(ab_hat))
    rm(out_med)
    
    # Process BSLMM
    out_bslmm <- read_rds(filename("bslmm", seed, vu))
    bslmm <- 
      out_bslmm |> 
      transmute(
        mediator = paste0("M", 1:n()), 
        ab_hat = beta_hat * alpha_hat,
        tie_hat = tie_hat[1],
        pip = pip,
        method = "BSLMM") |>
      left_join(x = d, by = "mediator")
    rm(out_bslmm)
    
    # Process Pathway LASSO
    load(filename("plasso", seed, vu))
    fdps_by_lambda <- with(out_plasso,
                           apply(AB.est,2,plasso_fdp,want_cols,
                                 d$true))
    which_lambda_use <- which(fdps_by_lambda < 0.1)[1]
    
    if(is.na(which_lambda_use)){
      ab_plasso <- with(out_plasso, 0 * want_cols)
      
    }else{
      ab_plasso <- with(out_plasso, AB.est[,which_lambda_use])
      
    }
    
    plasso <-
      with(out_plasso,
           tibble(
             mediator = paste0("M",want_cols),
             ab_hat = ab_plasso,
             pv = ifelse(ab_hat == 0, 1, 0.05)
           )
      ) |> 
      left_join(x = d, by = "mediator") |>
      mutate(across(c(ab_hat), ~ ifelse(is.na(.x), 0, .x)),
             pv = ifelse(is.na(pv), 1, pv),
             method = "P-LASSO",
             tie_hat = sum(ab_hat))
    rm(out_plasso, ab_plasso)

    
    # Combine methods
    combined <- 
      bind_rows(uni, hima, hdma, med, bslmm, plasso) |> 
      mutate(score = ifelse(is.na(pip), pv, 1 - pip)) |> 
      group_by(method) |> 
      mutate(sig = test_fdr(score, true)) |>
      summarize(
        np = sum(sig),
        ntp = sum(sig & true),
        nt = sum(true),
        tpr = sum(sig & true) / sum(true),
        nfp = sum(sig & !true),
        fpr = nfp / sum(!true),
        mse1 = mean((ab_hat[true] - ab[true])^2),
        mse0 = mean((ab_hat[!true] - ab[!true])^2),
        tie = sum(ab),
        tie_hat = tie_hat[1],
        tie_prb = abs(tie - tie_hat) / tie * 100
      ) |>
      ungroup() 
    
    # HILMA
    load(filename("hilma", seed, vu))
    hilma <-
      tibble(method = "HILMA",
             tie = with(d, sum(ab)),
             tie_hat = as.numeric(with(out_hilma, beta_hat * sdy / sda)),
             tie_prb = as.numeric(abs(tie - tie_hat) / tie * 100))
    rm(out_hilma)
    
    # PCMA
    load(filename("pcma", seed, vu))
    pcma <-
      tibble(method = "PCMA",
             tie = with(d, sum(ab)),
             tie_hat = out_pcma$IE.total[1],
             tie_prb = as.numeric(abs(tie - tie_hat) / tie * 100))
    rm(out_pcma)
    
    # PMED
    load(filename("pmed", seed, vu))
    pmed <-
      tibble(method = "PMED",
             tie = with(d, sum(ab)),
             tie_hat = pmed_tie,
             tie_prb = as.numeric(abs(tie - tie_hat) / tie * 100))
    rm(pmed_tie)
    
    results[[r]] <-
      bind_rows(combined, hilma, pcma, pmed) |> 
      mutate(vu = vu, seed = seed) |> 
      select(vu, seed, method, everything())
  
}

all_results <- bind_rows(results)



# Summarize the results over all replications
# TPR: true positive rate
# FPR: false positive rate
# MSE0: MSE among inactive mediators
# MSE1: MSE among active mediators
# RMSE1, RMSE0: Relative MSE compared to One-at-a-time ("univariate") method
# TIE_PRB: Percent relative bias for the global indirect effect
# _m suffix: mean over 100 replications
# _l, _u suffixes: lower and upper empirical 95% confidence limits

dat.sum <-
  dat |> 
  group_by(vu, seed) |> 
  mutate(
    rmse1 = mse1 / (mse1[method == "Univariate"]),
    rmse0 = mse0 / (mse0[method == "Univariate"])
  ) |> 
  group_by(vu, method) |> 
  summarize(
    across(c(tpr, fpr, rmse0, rmse1, tie_prb, mse1, mse0), 
           list(
             m = mean,
             l = ~ quantile(.x, 0.025, na.rm = T, names = F),
             u = ~ quantile(.x, 0.975, na.rm = T, names = F)))
  )


write_rds(dat.sum, "results/confounded_effects_results.rds")










