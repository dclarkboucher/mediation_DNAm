#Script for implementing a one-at-a-time mediation analysis on simulated data

#For efficiency one should run this on a computing cluster or other system
#allowing parallel operations. The script below assumes instead that there is
#only core responsible for all jobs, which would be relatively slow. To adopt
#the script to enable say, 50 cores or nodes, change the parameter "ncores" to 
#50 and assign "index" to be the index of the job, an integer 1 to 50. 

#Note that, before running this, one must have already generated data using 
#"generate_data.R". Make sure the "ndat" parameter below is not set to be more
#than the number of simulated datasets that were created. 
ndat <- 1 #Change to 100 for full simulation study

library(tidyverse)

tab <- expand_grid(setting = c("bl","hc","ns","ln_ns","ln", "ln_hc"),
                   dataset = c(1,2,3,4),
                   seed2 = 1:ndat)
index <- 1 #index of the current computing core being used
ncores <- 1 #total number of computing cores


which_run <- seq(index,nrow(tab),ncores)

for(row in which_run){
  l1 <- ls()
  
  #Select dataset
  d <- tab[row,"dataset"]
  setting <- tab[row,"setting"]
  seed2 <- tab[row,"seed2"]
  
  #Load dataset
  load(paste0("datasets/sim_data_",setting,"_d",d,"_s",seed2,".rda"))
  
  alphas <- map_dfr(as.data.frame(M), 
                    ~ summary(
                      lm(y ~ ., data = data.frame(y = .x, A = A))
                    )$coefficients[2,c(1,4)])
  alphas <-
    alphas %>% 
    rename(alpha_hat = 1, alpha_pv = 2) %>% 
    mutate(mediator = paste0("M",1:2000),
           alpha_pv_bonf = pmin(alpha_pv * n(),1),
           alpha_pv_fdr = p.adjust(alpha_pv, method = "fdr")) %>% 
    select(mediator,everything())
  
  
  betas <- map_dfr(as.data.frame(M), 
                   ~ summary(
                     lm(y ~ ., data = data.frame(y = Y, x = .x, A = A))
                   )$coefficients[2,c(1,4)])
  betas <-
    betas %>% 
    rename(beta_hat = 1, beta_pv = 2) %>% 
    mutate(mediator = paste0("M",1:2000),
           beta_pv_bonf = pmin(beta_pv * n(),1),
           beta_pv_fdr = p.adjust(beta_pv, method = "fdr")) %>% 
    select(mediator,everything())
  
  out_uni <-
    full_join(betas,alphas, by = "mediator") %>% 
    mutate(ab_hat = alpha_hat * beta_hat,
           ab_pv = pmax(alpha_pv,beta_pv),
           ab_pv_fdr = pmax(alpha_pv_fdr,beta_pv_fdr),
           ab_pv_fdr1 = p.adjust(ab_pv,method = "fdr"),
           ab_pv_bonf = pmin(ab_pv * n(), 1))
  
  save(out_uni, 
       file = paste0("simulation_results/one-at-a-time/sim_out_uni_",
                     setting,"_d",d,"_s",seed2,".rda"))
  
  #Remove unneeded things
  rm(list = setdiff(ls(),l1))
  
}


