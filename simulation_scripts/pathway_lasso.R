#Script for implementing pathway LASSO on simulated data

#For efficiency one should run this on a computing cluster or other system
#allowing parallel operations. The script below assumes instead that there is
#only core responsible for all jobs, which would be very slow. To adopt
#the script to enable say, 50 cores or nodes, change the parameter "ncores" to 
#50 and assign "index" to be the index of the job, an integer 1 to 50. 

#Note that, before running this, one must have already generated data using 
#"generate_data.R". Make sure the "ndat" parameter below is not set to be more
#than the number of simulated datasets that were created. 
ndat <- 1 #Change to 100 for full simulation study

index <- 1 #index of the current computing core being used
ncores <- 1 #total number of computing cores
library(tidyverse)

source("utils/plasso_functions.R")
dir.data <- paste0("simulation_results/plasso")


tab <- expand_grid(setting = c("bl","hc","ns","ln","ln_hc","ln_ns"),
                   dataset = c(1,2,3,4),
                   seed2 = 1:ndat) 


#which nodes to run on this core
which_run <- seq(index, nrow(tab), ncores)

row <- 1
for(row in which_run){
  l1 <- ls()
  
  #Select dataset
  d <- tab[row,"dataset"]
  setting <- tab[row,"setting"]
  seed2 <- tab[row,"seed2"]
  
  #Load dataset
  load(paste0("simulation_datasets/sim_data_",setting,"_d",d,"_s",seed2,".rda"))
  
  #Path
  run.file <- paste0(dir.data,"/RUN_",setting,"_d",d,"_s",seed2,".RData")
  
  #Screening based on Y~M model
  pvals <- 
    apply(M,2,
          function(x){
            summary(lm(y ~ x1 + a, data.frame(y = Y, x1 = x, a = A))
            )$coefficients[2,4]}
    )
  p <- ceiling((nrow(M))/log((nrow(M)))) 
  colnames(M) <- paste0("M",1:ncol(M))
  want_cols <- sort(order(pvals)[1:p])
  
  M_s <- M[,want_cols]
  
  
  #Input things to pathway LASSO
  Z <- as.vector(A)
  M <- M_s
  R <- Y
  
  #Initial things
  phi<-2
  rho<-1
  max.itr<-5000
  tol<-1e-6
  thred<-1e-6
  thred2<-1e-3
  
  k <- ncol(M)
  colnames(M) <- paste0("M",1:k)
  dd0 <- data.frame(Z=Z, M, R=R)
  Sigma10<-diag(rep(1,k))
  Sigma20<-matrix(1,1,1)
  
  
  # standardize data
  m.Z <- mean(Z)
  m.M <- apply(M,2,mean)
  m.R <- mean(R)
  sd.Z <- sd(Z)
  sd.M <- apply(M,2,sd)
  sd.R <- sd(R)
  
  Z <- scale(Z)
  M <- scale(M)
  R <- scale(R)
  dd <- data.frame(Z=Z, M, R=R)
  
  lambda<-c(10^c(#seq(-5,-3,length.out=5),
    seq(-3,0,length.out=41)[-1],
    seq(0,1,length.out=6)[-1])) # lambda values
  
  d <- 1 
  
  A.pf<-sd.Z/sd.M
  B.pf<-sd.M/sd.R
  C.pf<-sd.Z/sd.R
  AB.pf<-A.pf*B.pf
  
  re <- vector("list",length=length(lambda))
  AB.est = A.est = B.est<-matrix(NA,k,length(lambda))
  C.est<-rep(NA,length(lambda))
  
  #rate parameter for omega vs lambda (using 1-1, best in paper)
  omega.p <- c(0,0.1,1)
  omega.p.idx <- 3
  
  #fit model for each lambda
  for(i in length(lambda):1)
  {
    out<-NULL
    if(i==length(lambda))
    {
      # starting from the largest lambda value
      try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambda[i],
                                     omega=omega.p[omega.p.idx]*lambda[i],
                                     phi=phi,Phi1=NULL,Phi2=NULL,
                                     rho=rho,rho.increase=FALSE,
                                     tol=tol,max.itr=max.itr,thred=thred,
                                     Sigma1=Sigma10,Sigma2=Sigma20,
                                     trace=FALSE))
    }else
    {
      # for smaller lambda (ith lambda), use the (i+1)th lambda results as burn-in 
      try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambda[i],
                                     omega=omega.p[omega.p.idx]*lambda[i],
                                     phi=phi,Phi1=NULL,Phi2=NULL,
                                     rho=rho,rho.increase=FALSE,
                                     tol=tol,max.itr=max.itr,
                                     thred=thred,Sigma1=Sigma10,
                                     Sigma2=Sigma20,trace=FALSE,
                                     Theta0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),
                                     D0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1),
                                     alpha0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),
                                     beta0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1)))
    }
    
    if(is.null(out)==FALSE)
    {
      re[[i]]<-out
      
      # scale the estimate back to the original scale
      B.est[,i]<-out$B*(sd.R/sd.M)
      C.est[i]<-out$C*(sd.R/sd.Z)
      A.est[,i]<-out$A*(sd.M/sd.Z)
      AB.est[,i]<-A.est[,i]*B.est[,i]
    }
    
    print(paste0("lambda index ",i))
  }
  
  save(list=c("re","AB.est","A.est","B.est","C.est","A.pf","B.pf","C.pf","AB.pf",
              "want_cols"), file=run.file)
  
  #Remove unneeded things
  rm(list = setdiff(ls(),l1))
  
}


