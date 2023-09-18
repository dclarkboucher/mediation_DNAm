library(mvtnorm)

load("vcv_sparse.rda") #loads object "vcv_sparse"
vcv_sparse <- as.matrix(vcv_sparse)
vcv <- vcv_sparse
vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
rm(vcv_sparse)


get_data_pos <- function(seed = 5, seed2, 
                         dataset = c(1,2,3,4)){
  
  # Inputs
  # seed <- 5
  # dataset <- 1
  # seed2 <- 2
  # 
  d <- dataset[1]
  
  n <- 2500
  p <- 2000
  
  na <- 80
  nb <- 80
  nab <- 20
  
  pve_a1 <- 0.2
  pve_de1 <- 0.1
  pve_ie1 <- 0.1
  
  va <- 1
  vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
  r <- 1 #parameter to regularize variance-covariance
  
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
  alpha_a[which_a] <- as.numeric(scale(rnorm(na)))
  
  #Sample beta
  beta_m <- rep(0,p)
  which_ab <- sort(sample(which_a,nab))
  which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
  beta_m[which_b] <- as.numeric(scale(rnorm(nb)))
  
  # Original effects
  ab <- alpha_a * beta_m
  tie <- sum(ab)
  
  # Absolute effects
  alpha_a1 <- abs(alpha_a)
  beta_m1 <- abs(beta_m)
  ab1 <- alpha_a1 * beta_m1
  tie1 <- sum(ab1)
  
  # Determine the error variance-covariance of the mediators
  var_m <- rep(vm, p) # error variance of Ms
  var_m[which_a] <- va * alpha_a1[which_a]^2 * (1 - pve_a) / pve_a
  scale1 <- sqrt(var_m)
  
  # vcv <- bigreadr::fread2("~/mediation/vcov2000.csv")
  set.seed(seed2)
  which_mediators <- (sample(2000,p))
  vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
  vcv1 <- cov2cor(vcv1)
  diag(vcv1) <- diag(vcv1) + r
  vcv1 <- cov2cor(vcv1)
  vcv2 <- t(vcv1 * scale1) * scale1
  
  # The residual variance should be the same in both situations (positive 
  # effects or mixed effects). This makes the comparison as fair as possible
  
  # First compute the marginal and error variance of Y for the mixed case
  var_ie <- tie^2 * va
  var_y <- var_ie / pve_ie # desired marginal variance of Y
  beta_a <- sqrt(var_y * pve_de / va) * sign(tie) # direct effect
  vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
  var_res <- var_y - (beta_a + tie)^2 * va - vc4 # residual variance of Y

  # Now compute the rest of what we need. 
  var_res1 <- var_res
  beta_a1 <- tie1 * sqrt(pve_de / pve_ie) 
  
  #SAMPLE DATASET
  A <- as.numeric(scale(rnorm(n))) * sqrt(va)
  EM <- mvtnorm::rmvnorm(n, rep(0,p), vcv2)
  M <- A %*% t(alpha_a1) + EM
  ey <- rnorm(n, 0, sqrt(var_res1))
  Y <- A * beta_a1 + M %*% beta_m1 + ey
  
  return(list(
    Y = Y, 
    M = M, 
    A = A,
    alpha_a = alpha_a1, 
    beta_m = beta_m1,
    beta_a = beta_a1
  ))
  
  
}