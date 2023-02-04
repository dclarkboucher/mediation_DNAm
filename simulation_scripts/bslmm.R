#Script for implementing BSLMM on simulated data

#For efficiency one should run this on a computing cluster or other system
#allowing parallel operations. The script below assumes instead that there is
#only core responsible for all jobs, which would be very slow. To adopt
#the script to enable say, 50 cores or nodes, change the parameter "ncores" to 
#50 and assign "index" to be the index of the job, an integer 1 to 50. 

#Note that, before running this, one must have already generated data using 
#"generate_data.R". Make sure the "ndat" parameter below is not set to be more
#than the number of simulated datasets that were created. 
ndat <- 1 #Change to 100 for full simulation study

library(tidyr) #install.packages("tidyr")
library(bama) #install.packages("bama")

index <- 1 #index of the current computing core being used
ncores <- 1 #total number of computing cores

tab <- as.data.frame(expand_grid(setting = c("bl","hc","ns", "ln","ln_hc","ln_ns"),
                                 dataset = c(1,2,3,4),
                                 seed2 = 1:ndat))

which_run <- seq(index, nrow(tab), ncores)

for(row in which_run){
  l1 <- ls()
  
  #Select dataset
  d <- tab[row,"dataset"]
  setting <- tab[row,"setting"]
  seed2 <- tab[row,"seed2"]
  
  #Load dataset
  load(paste0("simulation_datasets/sim_data_",setting,"_d",d,"_s",seed2,".rda"))
  
  C <- matrix(1,nrow(M),1) # covariate matrix
  
  # Bslmm tends to work better when scaling, re-scaling
  y_sd <- sd(Y)
  m_sds <- apply(M, 2, sd)
  a_sd <- sd(A)
  M <- as.matrix(scale(M))
  Y <- as.numeric(scale(Y))
  A <- as.numeric(scale(A))
  
  # Get one-at-a-time coefficients so that we can input a good prior variance
  alpha_hats <- apply(M, 2, function(x){summary(lm(x ~ A))$coefficients[2,1]})
  alpha_var1 <- var(alpha_hats[abs(alpha_hats) > quantile(abs(alpha_hats),0.9)])
  beta_hats <- apply(M, 2, function(x){summary(lm(Y ~ x + A))$coefficients[2,1]})
  beta_var1 <- var(beta_hats[abs(beta_hats) > quantile(abs(beta_hats),0.9)])
  
  out_bslmm <-
    bama(
      Y = as.vector(Y),
      A = A,
      C1 = C,
      C2 = C,
      M = M,
      burnin = 15000,
      ndraws = 20000,
      method = "BSLMM",
      control = list(k = 2, lma1 = alpha_var1, lm1 = beta_var1,
                     l = 1, lambda0 = 0.04, lambda1 = 0.2, lambda2 = 0.2, 
                     phi0 = 0.01, phi1 = 0.01, a0 = 0.01 * ncol(M),
                     a1 = 0.05 * ncol(M), a2 = 0.05 * ncol(M),
                     a3 = 0.89 * ncol(M)),
      seed = 123
    )
  
   out <- 
    data.frame(
      beta_hat = colMeans(out_bslmm$beta.m) * y_sd / m_sds,
      alpha_hat = colMeans(out_bslmm$alpha.a) * m_sds / a_sd,
      tie_hat =
        c(with(out_bslmm, mean(rowSums(alpha.a * beta.m))) * y_sd / a_sd,
          rep(NA,1999)),
      pip = colMeans(with(out_bslmm,r1 * r3))
    )
  
  write.csv(out, 
       file = paste0("simulation_results/bslmm/sim_out_bslmm_",setting,"_d",d,
                     "_s",seed2,".csv"))
  
 
  
  #Remove unneeded things
  rm(list = setdiff(ls(),l1))
  
}



