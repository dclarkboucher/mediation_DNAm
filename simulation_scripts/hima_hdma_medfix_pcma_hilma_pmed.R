#Script for implementing HIMA, HDMA, MedFix, PCMA, HILMA, and PMED on simulated data

#For efficiency one should run this on a computing cluster or other system
#allowing parallel operations. The script below assumes instead that there is
#only core responsible for all jobs, which would be very slow. To adopt
#the script to enable say, 50 cores or nodes, change the parameter "ncores" to 
#50 and assign "index" to be the index of the job, an integer 1 to 50. 

#Note that, before running this, one must have already generated data using 
#"generate_data.R". Make sure the "ndat" parameter below is not set to be more
#than the number of simulated datasets that were created. 
ndat <- 1 #Change to 100 for full simulation study

#MAKE SURE OUTPUT MATCHES WHAT'S EXPECTED BY LATER FUNCTIONS

library(freebird) #install.packages("freebird")
library(spcma) #devtools::install_github("https://github.com/zhaoyi1026/spcma")
library(tidyr)


source("utils/hima_hdma_utils.R") #may need other packages
source("utils/medfix_functions.R") #may need other packages
source("utils/pmed_functions.R") #may need other packages

index <- 1 #index of the current computing core being used
ncores <- 1 #total number of computing cores



tab <- expand_grid(setting = c("bl","hc","ns","ln_ns","ln", "ln_hc"),
                   dataset = c(1,2,3,4),
                   seed2 = 1:ndat)

#which nodes to run on this core
which_run <- seq(index, nrow(tab), ncores)

for(row in which_run){
  l1 <- ls()
  
  #Select dataset
  d <- tab[row,"dataset"]
  setting <- tab[row,"setting"]
  seed2 <- tab[row,"seed2"]
  
  #Load dataset
  load(paste0("simulation_datasets/sim_data_",setting,"_d",
              d,"_s",seed2,".rda"))
  
  # Number of mediators for SIS screening, where used
  p <- ceiling(nrow(M) / log(nrow(M))) 
  
  #HIMA
  set.seed(123)
  out_hima <-
    hima(
      X = as.vector(A),
      Y = as.vector(Y),
      M = as.data.frame(M),
      topN = p,
      Y.family = "gaussian",
      M.family = "gaussian",
      penalty = "MCP",
      max.iter = 10^5,
      parallel = T
    )
  
  save(out_hima,
       file = paste0("simulation_results/hima/sim_out_hima_",setting,"_d",d,
                     "_s",seed2,".rda"))
  
  #MedFix
  set.seed(123)
  out_med <- medfix(Y,M,A,p=p,screen = "ym")
  save(out_med,
       file = paste0("simulation_results/med/sim_out_med_",setting,"_d",d,
                     "_s",seed2,".rda"))
  
  #HDMA
  set.seed(123)
  out_hdma <-
    hdma(
      X = as.vector(A),
      Y = as.vector(Y),
      M = as.data.frame(M),
      family = "gaussian",
      topN = p,
      method = "lasso",
    )

  save(out_hdma,
       file = paste0("simulation_results/hdma/sim_out_hdma_", setting,"_d",d,
                     "_s",seed2,".rda"))
  
  #PCMA
  set.seed(123)
  try(out_pcma <-
    try(mcma_PCA(
      X = as.vector(A),
      M = as.matrix(M),
      Y = as.vector(Y),
      adaptive = F,
      n.pc = 100,
      boot = TRUE,
      sims = 2,
    ))) 

  try(save(out_pcma,
       file = paste0("simulation_results/pcma/sim_out_pcma_", setting,"_d",d,
                     "_s",seed2,".rda")))
  
  #HILMA
  #Note that HILMA has much better performance when the data are standardized
  #in advance. The resulting global indirect effect is then multiplied by
  #sdy/sda
  set.seed(123)
  Y1 <- scale(Y)
  M1 <- scale(M)
  A1 <- scale(A)
  
  out_hilma <-  hilma(as.vector(Y1), as.matrix(M1), matrix(A1, ncol = 1), center = T)
  
  save(out_hilma,
       file = paste0("simulation_results/hilma/sim_out_hilma_",setting,"_d",d,
                     "_s",seed2,".rda"))
  
  #PMED
  #Note that PMED fails with error if no mediators are selected, which we handle
  #using the "try" function
  set.seed(123)
  out_pmed <- NULL
  C <- matrix(1, nrow(M), 1)
  coef <- glmnet::coef.glmnet
  try({out_pmed <- pmed_wrapper(A = A, M = M, Y = Y, C = C, nlambda = 50)})
  
  if(is.null(out_pmed)){
    pmed_tie <- 0
  }else{
    pmed_tie <- out_pmed$summary_result[2,2]
  }
  
  
  save(pmed_tie,
       file = paste0("simulation_results/pmed/sim_out_pmed_",setting,"_d",d,
                     "_s",seed2,".rda"))
  
  
  #Remove unneeded things
  rm(list = setdiff(ls(),l1))
  
}





