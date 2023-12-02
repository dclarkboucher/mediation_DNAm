# Script for generating datasets used in simulations
# To run this file, open the R project "hdmediation" located in the outermost
# directory of this repository. This should put you in the correct location

# To generate data for the whole simulation study, set ndat below to 100 (this
# will take a long time). For an easy and comparatively efficient example, set 
# ndat to 1

#The following functions have 3 main parameters which are also referenced
#in later files:
#"dataset" indexes the 4 possible choices of pve_a, pve_ie, and pve_de
#"setting" denotes the sample size, sparsity, and correlations (see below)
#"seed2" denotes the random seed used for sampling

#"setting" parameter abbreviations:
#bl: baseline correlations and sparsity
#ln: low sample size (n = 1000)
#hc: high-correlation between mediators compared to baseline
#ns: non-sparse instead of sparse

# Make folders where output will be stored.
dir.create("simulation_datasets")
dir.create("simulation_results")
for(name in c("bslmm","hdma","hilma","hima","med","one-at-a-time",
              "pcma","plasso","pmed")){
  dir.create(paste0("simulation_results/",name))
}


ndat <- 1 # 100 for full study

{
  library(mvtnorm)
  library(tidyverse)
  library(Matrix)
}

#load variance-covariance matrix, which is stored as a sparse matrix with the
#lower triangle set to zero
load("vcv_sparse.rda") #loads object "vcv_sparse"
vcv_sparse <- as.matrix(vcv_sparse)
vcv <- vcv_sparse
vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
rm(vcv_sparse)

#generates data for the n = 2500 case
get_data <- function(seed = 5, seed1 = 35, seed2, 
                     setting = c("bl","hc","ns"),
                     dataset = c(1,2,3,4)){
  
  
  setting <- setting[1]
  d <- dataset[1]
  
  if(setting %in% c("bl","hc")){
    
    seed <- seed
    
    n <- 2500
    p <- 2000
    
    na <- 80
    nb <- 80
    nab <- 20
    
    pve_a1 <- 0.2
    pve_de1 <- 0.1
    pve_ie1 <- 0.1
    
    va <- 1
    valpha <- 1
    vbeta <- 1
    vm <- 1 #variance of mediators with alpha = 0. 
    r <- 1 #parameter to regularize variance-covariance
    
    if(setting == "hc"){
      r <- 0.1
    }
    
    
    pve_a <- pve_a1
    pve_de <- pve_de1
    pve_ie <- pve_ie1
    
    if(d == 2){
      pve_a <- pve_a / 2
      
    }else if(d == 3){
      pve_de <- pve_de / 2
      
    }else if(d == 4){
      pve_ie <- pve_ie / 2
      
    }
    
    #SAMPLE POPULATION ATTRIBUTES
    set.seed(seed)
    
    #Sample alpha
    alpha_a <- rep(0,p)
    which_a <- sort(sample(p,na))
    alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
    #Sample beta
    beta_m <- rep(0,p)
    which_ab <- sort(sample(which_a,nab))
    which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
    beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
    
    ab <- alpha_a * beta_m
    tie <- sum(ab)
    
    #Load mediator variance-covariance
    # vcv <- bigreadr::fread2("mediation/vcov2000.csv")
    set.seed(seed2)
    which_mediators <- (sample(2000,p))
    vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
    #standardize mediators
    sds <- sqrt(diag(vcv1))
    vcv1 <- t(vcv1 / sds) / sds
    
    #Regularize variance-covariance matrix so that it's not singular
    diag(vcv1) <- diag(vcv1) + r
    
    #Set desired mediator variance
    scale <- sqrt(diag(vcv1)) / sqrt(vm)
    vcv1 <- t(vcv1 / scale) / scale
    
    #Scale certain mediators to get desired pve_a
    scale1 <- rep(1,p)
    scale1[which_a] <- sqrt( ((va * alpha_a[which_a] ^ 2) / pve_a -
                                va * alpha_a[which_a]^2)/ vm)
    vcv2 <- t(vcv1 * scale1) * scale1
    
    #Compute Y variance components
    var_ie <- va * sum(ab)^2
    var_y <- var_ie / pve_ie #desired variance of Y
    beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
    vc2 <- 2 * beta_a * sum(ab) * va
    vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
    var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
    if(var_res <= 0){
      message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
      break
    }
    
    #SAMPLE DATASET
    A <- as.numeric(scale(rnorm(n))) * sqrt(va)
    
    EM <- rmvnorm(n,rep(0,p),vcv2)
    M <- A %*% t(alpha_a) + EM
    
    ey <- rnorm(n,0,sqrt(var_res))
    Y <- A * beta_a + M %*% beta_m + ey
    
    
    save(Y, M, A, alpha_a, beta_m, beta_a, pve_a, pve_de, pve_ie,
         file = paste0("simulation_datasets/sim_data_",
                       setting,"_d",d,"_s",seed2,".rda"))
    
  }
  
  if(setting == "ns"){
    
    #ASSIGN PARAMETERS
    
    n <- 2500
    p <- 2000
    
    na <- 80
    nb <- 80
    nab <- 20
    
    pve_a1 <- 0.2
    pve_de1 <- 0.1
    pve_ie1 <- 0.1
    
    va <- 1
    valpha <- 1
    vbeta <- 1
    vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
    r <- 1 #parameter to regularize variance-covariance
    
    vbeta0 <- 0.2^2 #difference needs to be small enough for BSLMM to struggle
    valpha0 <- 0.2^2
    
    pve_a <- pve_a1
    pve_de <- pve_de1
    pve_ie <- pve_ie1
    
    if(d == 2){
      pve_a <- pve_a / 2
      
    }else if(d == 3){
      pve_de <- pve_de / 2
      
    }else if(d == 4){
      pve_ie <- pve_ie / 2
      
    }
    
    #SAMPLE POPULATION ATTRIBUTES
    
    set.seed(seed)
    #Sample alpha
    alpha_a <- rep(0,p)
    which_a <- sort(sample(p,na))
    alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
    #Sample beta
    beta_m <- rep(0,p)
    which_ab <- sort(sample(which_a,nab))
    which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
    beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
    
    #Fill in smaller effects
    set.seed(seed1)
    alpha_a[-which_a] <- as.numeric(scale(rnorm(p - na))) * sqrt(valpha0)
    beta_m[-which_b] <- as.numeric(scale(rnorm(p - nb))) * sqrt(vbeta0)
    
    ab <- alpha_a * beta_m
    tie <- sum(ab); tie
    
    # }
    
    sum(alpha_a[which_ab] * beta_m[which_ab])
    
    #Load mediator variance-covariance
    # vcv <- bigreadr::fread2("vcov2000.csv")
    set.seed(seed2)
    which_mediators <- (sample(2000,p))
    vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
    #standardize mediators
    sds <- sqrt(diag(vcv1))
    vcv1 <- t(vcv1 / sds) / sds
    
    #Regularize variance-covariance matrix so that it's not singular
    diag(vcv1) <- diag(vcv1) + r
    
    #Set desired mediator variance
    scale <- sqrt(diag(vcv1)) / sqrt(vm)
    vcv1 <- t(vcv1 / scale) / scale
    
    #Scale certain mediators to get desired pve_a
    scale1 <- sqrt( ((va * alpha_a ^ 2) / pve_a - va * alpha_a^2)/ vm)
    vcv2 <- t(vcv1 * scale1) * scale1
    
    #Compute Y variance components
    var_ie <- va * sum(ab)^2
    var_y <- var_ie / pve_ie #desired variance of Y
    beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
    vc2 <- 2 * beta_a * sum(ab) * va
    vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
    var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
    if(var_res <= 0){
      message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
      break
    }
    
    #SAMPLE DATASET
    A <- as.numeric(scale(rnorm(n))) * sqrt(va)
    
    EM <- rmvnorm(n,rep(0,p),vcv2)
    M <- A %*% t(alpha_a) + EM
    
    ey <- rnorm(n,0,sqrt(var_res))
    Y <- A * beta_a + M %*% beta_m + ey
    
    #Check variance of Y
    var(Y)
    var_y
    
    #Check proportion of variance explained
    alpha_a[which_a[1:5]]^2 / apply(M[,which_a[1:5]],2,var)
    
    save(Y, M, A, alpha_a, beta_m, beta_a, pve_a, pve_de, pve_ie,
         which_a, which_ab,which_b,
         file = paste0("simulation_datasets/sim_data_",
                       setting,"_d",d,"_s",seed2,".rda"))
    
    
    
  }
  
}

#generates data for a flexible choice of N
get_data_flex <- function(seed = 5, seed1 = 35, seed2, 
                          setting = c("ln","ln_hc","ln_ns"),
                          dataset = c(1,2,3,4), n = 1000){
  
  
  setting <- setting[1]
  d <- dataset[1]
  
  if(setting %in% c("ln","ln_hc")){
    
    seed <- seed
    
    p <- 2000
    
    na <- 80
    nb <- 80
    nab <- 20
    
    pve_a1 <- 0.2
    pve_de1 <- 0.1
    pve_ie1 <- 0.1
    
    va <- 1
    valpha <- 1
    vbeta <- 1
    vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
    r <- 1 #parameter to regularize variance-covariance
    
    if(setting == "ln_hc"){
      r <- 0.1
    }
    
    
    pve_a <- pve_a1
    pve_de <- pve_de1
    pve_ie <- pve_ie1
    
    if(d == 2){
      pve_a <- pve_a / 2
      
    }else if(d == 3){
      pve_de <- pve_de / 2
      
    }else if(d == 4){
      pve_ie <- pve_ie / 2
      
    }
    
    #SAMPLE POPULATION ATTRIBUTES
    set.seed(seed)
    
    #Sample alpha
    alpha_a <- rep(0,p)
    which_a <- sort(sample(p,na))
    alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
    #Sample beta
    beta_m <- rep(0,p)
    which_ab <- sort(sample(which_a,nab))
    which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
    beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
    
    ab <- alpha_a * beta_m
    tie <- sum(ab)
    
    #Load mediator variance-covariance
    # vcv <- bigreadr::fread2("vcov2000.csv")
    set.seed(seed2)
    which_mediators <- (sample(2000,p))
    vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
    #standardize mediators
    sds <- sqrt(diag(vcv1))
    vcv1 <- t(vcv1 / sds) / sds
    
    #Regularize variance-covariance matrix so that it's not singular
    diag(vcv1) <- diag(vcv1) + r
    
    #Set desired mediator variance
    scale <- sqrt(diag(vcv1)) / sqrt(vm)
    vcv1 <- t(vcv1 / scale) / scale
    
    #Scale certain mediators to get desired pve_a
    scale1 <- rep(1,p)
    scale1[which_a] <- sqrt( ((va * alpha_a[which_a] ^ 2) / pve_a -
                                va * alpha_a[which_a]^2)/ vm)
    vcv2 <- t(vcv1 * scale1) * scale1
    
    #Compute Y variance components
    var_ie <- va * sum(ab)^2
    var_y <- var_ie / pve_ie #desired variance of Y
    beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
    vc2 <- 2 * beta_a * sum(ab) * va
    vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
    var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
    if(var_res <= 0){
      message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
      break
    }
    
    #SAMPLE DATASET
    A <- as.numeric(scale(rnorm(n))) * sqrt(va)
    
    EM <- rmvnorm(n,rep(0,p),vcv2)
    M <- A %*% t(alpha_a) + EM
    
    ey <- rnorm(n,0,sqrt(var_res))
    Y <- A * beta_a + M %*% beta_m + ey
    
    #Check variance of Y
    var(Y)
    var_y
    
    save(Y, M, A, alpha_a, beta_m, beta_a, pve_a, pve_de, pve_ie,
         file = paste0("simulation_datasets/sim_data_",
                       setting,"_d",d,"_s",seed2,".rda"))
    
  }
  
  if(setting == "ln_ns"){
    
    #ASSIGN PARAMETERS
    seed <- seed
    seed1 <- seed1 #chosen to make residual variance positive
    p <- 2000
    
    na <- 80
    nb <- 80
    nab <- 20
    
    pve_a1 <- 0.2
    pve_de1 <- 0.1
    pve_ie1 <- 0.1
    
    va <- 1
    valpha <- 1
    vbeta <- 1
    vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
    r <- 1 #parameter to regularize variance-covariance
    
    vbeta0 <- 0.2^2 #difference needs to be small enough for BSLMM to struggle
    valpha0 <- 0.2^2
    
    pve_a <- pve_a1
    pve_de <- pve_de1
    pve_ie <- pve_ie1
    
    if(d == 2){
      pve_a <- pve_a / 2
      
    }else if(d == 3){
      pve_de <- pve_de / 2
      
    }else if(d == 4){
      pve_ie <- pve_ie / 2
      
    }
    
    #SAMPLE POPULATION ATTRIBUTES
    
    set.seed(seed)
    #Sample alpha
    alpha_a <- rep(0,p)
    which_a <- sort(sample(p,na))
    alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
    #Sample beta
    beta_m <- rep(0,p)
    which_ab <- sort(sample(which_a,nab))
    which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
    beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
    
    #Fill in smaller effects
    set.seed(seed1)
    alpha_a[-which_a] <- as.numeric(scale(rnorm(p - na))) * sqrt(valpha0)
    beta_m[-which_b] <- as.numeric(scale(rnorm(p - nb))) * sqrt(vbeta0)
    
    ab <- alpha_a * beta_m
    tie <- sum(ab); tie
    
    # }
    
    sum(alpha_a[which_ab] * beta_m[which_ab])
    
    #Load mediator variance-covariance
    # vcv <- bigreadr::fread2("vcov2000.csv")
    set.seed(seed2)
    which_mediators <- (sample(2000,p))
    vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
    #standardize mediators
    sds <- sqrt(diag(vcv1))
    vcv1 <- t(vcv1 / sds) / sds
    
    #Regularize variance-covariance matrix so that it's not singular
    diag(vcv1) <- diag(vcv1) + r
    
    #Set desired mediator variance
    scale <- sqrt(diag(vcv1)) / sqrt(vm)
    vcv1 <- t(vcv1 / scale) / scale
    
    #Scale certain mediators to get desired pve_a
    scale1 <- sqrt( ((va * alpha_a ^ 2) / pve_a - va * alpha_a^2)/ vm)
    vcv2 <- t(vcv1 * scale1) * scale1
    
    #Compute Y variance components
    var_ie <- va * sum(ab)^2
    var_y <- var_ie / pve_ie #desired variance of Y
    beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
    vc2 <- 2 * beta_a * sum(ab) * va
    vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
    var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
    if(var_res <= 0){
      message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
      break
    }
    
    #SAMPLE DATASET
    A <- as.numeric(scale(rnorm(n))) * sqrt(va)
    
    EM <- rmvnorm(n,rep(0,p),vcv2)
    M <- A %*% t(alpha_a) + EM
    
    ey <- rnorm(n,0,sqrt(var_res))
    Y <- A * beta_a + M %*% beta_m + ey
    
    #Check variance of Y
    var(Y)
    var_y
    
    #Check proportion of variance explained
    alpha_a[which_a[1:5]]^2 / apply(M[,which_a[1:5]],2,var)
    
    save(Y, M, A, alpha_a, beta_m, beta_a, pve_a, pve_de, pve_ie,
         which_a, which_ab,which_b,
         file = paste0("simulation_datasets/sim_data_",
                       setting,"_d",d,"_s",seed2,".rda"))
    
    
    
  }
}

#generate N = 2500 data
tab <- expand_grid(setting = c("ln_hc","ln","ln_ns"),dataset = c(1,2,3,4),
                   seed2 = 1:ndat)

pwalk(tab, ~ get_data_flex(setting = ..1, dataset = ..2, 
                           seed2 = ..3))

#generate N = 1000 data
tab <- expand_grid(setting = c("bl","hc","ns"),dataset = c(1,2,3,4),
                   seed2 = 1:ndat)

pwalk(tab, ~ get_data(setting = ..1, dataset = ..2, 
                      seed2 = ..3))


