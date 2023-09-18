library(mvtnorm)

load("vcv_sparse.rda") #loads object "vcv_sparse"
vcv_sparse <- as.matrix(vcv_sparse)
vcv <- vcv_sparse
vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
rm(vcv_sparse)

get_data_umc <- function(seed = 5, seed2, vu = 1){
  
  # Set up
  n <- 2500
  p <- 2000
  
  na <- 80
  nb <- 80
  nab <- 20
  
  # Proportion of variance explained
  pve_a <- 0.2
  pve_de <- 0.1
  pve_ie <- 0.1
  
  # SAMPLE POPULATION ATTRIBUTES
  set.seed(seed)
  
  # Sample alpha
  valpha <- 1
  vbeta <- 1
  alpha_a <- rep(0,p)
  which_a <- sort(sample(p,na))
  alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
  
  # Sample beta
  beta_m <- rep(0,p)
  which_ab <- sort(sample(which_a,nab))
  which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
  beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
  
  ab <- alpha_a * beta_m
  tie <- sum(ab)
  
  # Determine mediator variances
  tau <- (1 - pve_a) * (alpha_a^2) / pve_a # Error variance of Ms
  tau[which(alpha_a == 0)] <- 1
  tau_root <- sqrt(tau)
  var_m <- tau + alpha_a^2 # Variance of M|U
  
  # Load mediator variance-covariance
  # vcv <- bigreadr::fread2("vcov2000.csv")
  set.seed(seed2)
  which_mediators <- (sample(2000,p))
  vcv1 <- as.matrix(vcv[which_mediators,which_mediators]) |> cov2cor()
  diag(vcv1) <- 2
  vcv1 <- cov2cor(vcv1)
  vcv2 <- t(vcv1 * tau_root) * tau_root
  
  # Determine Y variance, effects
  v1 <- tie^2 / pve_ie # Variance of Y
  beta_a <- sqrt(pve_de * v1) * sign(tie) # Direct effect
  C <-  sum( (beta_m %*% t(beta_m)) * vcv2) # Covariance term
  omega <- v1 - (tie + beta_a)^2 - C # Error variance of Y
  # omega
  
  # Confounder effects
  gamma_u <- 1/3
  alpha_u <- rep(0, p) 
  alpha_u[which_a] <- (alpha_a[which_a]) / 2
  alpha_u[setdiff(which_b, which_ab)] <- 1/2
  beta_u <- beta_a / 2
  
  # Sample dataset without confounder
  U <- as.numeric(scale(rnorm(n))) * sqrt(vu)
  A <- as.numeric(scale(rnorm(n))) + U * gamma_u
  Delta <- rmvnorm(n,rep(0,p),vcv2)
  M <- A %*% t(alpha_a) + U %*% t(alpha_u) + Delta
  psi <- rnorm(n, 0, sqrt(omega))
  Y <- A * beta_a + M %*% beta_m + U * beta_u + psi
  
  # Add confounder to each variable
  return(list(
    Y = Y, 
    M = M, 
    A = A,
    U = U,
    alpha_a = alpha_a, 
    beta_m = beta_m,
    beta_a = beta_a
  ))
  
}