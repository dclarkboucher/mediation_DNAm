library(glmnet)


pmed_wrapper <- function(A, M, Y, C = NULL, nlambda = 100){
  
  
  lam_list <- exp(seq(log(0.000001), log(30), length.out = nlambda))
  A <- matrix(A, ncol = 1)
  
  hdMediation(X = A, y = Y, M = M, lam_list = lam_list, S = C)
  
}


# Utilities ---------------------------------------------------------------

hdMediation <-function(X, y, M, lam_list, S = NULL, n_imp = 0){
  # This function implements the high-dimensional mediation analysis
  #' @param X (Required) n by q exposure matrix. q can be 1, and 
  #'                     q < n is required
  #' @param y (Required) n-dimensional outcome vector.
  #' @param M (Required) n by p mediator matrix. p can be larger than n.
  #' @param lam_list (Required) a list of tuning parameter for HBIC
  #' @param S (Optional) n by s confounding variables matrix. s can be 1, 
  #'                     and s < n is required. 
  #' @param n_imp (Optional) an int, specifying the number of unpenalized 
  #'                        mediators in the last n_imp columns of M matrix 
  #'                      
  #' @return
  #'    A list of class `hdMediation'
  #'    - Mediators_imp: Selected important mediators' name`
  #'    - summary_result: a data frame of estimated coefficients, 
  #'                      SE, Test_statistics	p_value
  
  M_A <- HBIC_PPLS(X, y, M, S = S, lam_list = lam_list, n_imp =n_imp)
  # print(str(M_A)) 
  direct_indirect_eff <- mediationInference(X, y, M_A, S)
  
  
  return(list(Mediators_imp = colnames(M_A), 
              summary_result = direct_indirect_eff$summary_result))
}


HBIC_PPLS <- function(X, y, M, S = NULL, lam_list =NULL, n_imp =8, topPotential= 1000){
  #' This function implements solving the partially penalized least squares with the SCAD penalty. 
  #' The tuning parameter lambda for the penalty function is chosen based on the high-dimensional 
  #' BIC (HBIC) method.
  #' @param X The n by q exposure matrix. q can be 1, and q < n is required
  #' @param y The n-dimensional outcome vector.
  #' @param M The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required. 
  #' @param lam_list a list of tuning parameter for HBIC
  #' @param topPotentialtop potential mediators that is retained from screening stage
  #' @return
  #'    A dataframe of selected important mediators (denoted as M_{\hat{A}} in the paper)
  hbic= c()
  
  result =lapply(lam_list, HBIC_calc, xx=X,yy=y, mm=M, S = S, n_imp = n_imp)
  
  for( ii in 1: length(lam_list)){
    hbic[ii] = result[[ii]]$BIC
  }
  # Find minimum HBIC score's corresponding lambda
  id = which(hbic==min(hbic))
  id = tail(id,1)
  lamb = lam_list[id]
  # print(paste("choose lambda:", lamb)) 
  result = result[[id]]
  alpha0_hat = result$alpha0
  alpha1_hat = result$alpha1
  
  A = which(alpha0_hat!=0)
  # Selected mediators 
  M_A = M[,c(A[A>topPotential],A[A<=topPotential])]
  return(M_A)
}

mediationInference<-function(X, Y, M, S){
  #' This function implements statistical inference of mediation model
  #' @param X The n by q exposure matrix. q can be 1, and q < n is required
  #' @param y The n-dimensional outcome vector.
  #' @param M The n by p selected important mediator matrix. p<n
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' 
  
  #' @return
  #'    A list of class `mediationInference`:
  #'    - beta_hat: estimated indirect effect
  #'    - alpha1_hat: estimated direct effect
  #'    - Sn: test statistc value of indirect effect
  #'    - Tn: test statistc value of direct effect
  #'    - summary_result: (dataframe) summary of The estimated coefficients, standard errors, test statistics values and p-values
  
  
  p = ncol(M)
  q = ncol(X)
  n = nrow(X)
  if(length(S) ==0){
    Z = cbind(M,X)
    s = 0
    # print(dim(M)) 
    alpha0_tld = solve(t(M)%*%M)%*%t(M)%*%Y
    RSS02 = t(Y - M%*% alpha0_tld) %*% (Y - M%*%alpha0_tld)}
  else{
    Z = cbind(M,X,S)
    s = ncol(S)
    MS = cbind(M,S)
    alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)}
  
  alpha_rf = solve(t(Z)%*%Z)%*%t(Z)%*%Y
  alpha0_hat = alpha_rf[1:p]
  alpha1_hat = alpha_rf[(p+1):(p+q)]
  alpha2_hat = alpha_rf[(p+q+1):(p+q+s)]
  
  res = Y - Z%*%alpha_rf
  RSS12 = as.numeric(t(res) %*% (res)) # Test direct effect
  
  # RSS01 = t(Y - X%*%alpha1_hat) %*% (Y - X%*%alpha1_hat) # Test indirect effect
  # RSS11 = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)
  
  # degree of freedom
  df = p+q + s
  sigma1_hat = RSS12/(n - df)
  Sigma_MM = t(M)%*%M /n
  if(s == 0){
    gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y 
    sigmaT_hat = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)/(n-q)
    beta_hat = gamma_hat -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    
    invXX = solve(t(X)%*%X/n)
    Sigma_MX =t(M)%*%X/n
    
    B = invXX %*%t(Sigma_MX) %*%solve(Sigma_MM -  Sigma_MX%*%invXX %*% t(Sigma_MX)) %*%Sigma_MX %*%invXX 
    var_alpha1_hat = sigma1_hat*(invXX + B)
    cov_beta_hat = sigma2_hat * invXX + sigma1_hat * B
  }else{
    V = cbind(X,S)
    gamma_hat = solve(t(V)%*%V)%*%t(V)%*%Y 
    sigmaT_hat = t(Y-V%*%gamma_hat) %*% (Y-V%*%gamma_hat)/(n-q-s)
    beta_hat = gamma_hat[1:q] -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    Sigma_VV = t(V)%*%V/n
    invVV = solve(Sigma_VV)
    Sigma_MV =t(M)%*%V/n
    Sigma_VM =t(Sigma_MV)
    B = invVV %*% Sigma_VM %*%solve(Sigma_MM - Sigma_MV %*% invVV %*% Sigma_VM)%*%Sigma_MV%*%invVV
    var_alpha1_hat = sigma1_hat*(invVV + B)[1:q,1:q]
    cov_beta_hat = sigma2_hat * invVV + sigma1_hat * B
    cov_beta_hat = cov_beta_hat[1:q, 1:q]
  }
  # Test for beta
  # Wald's test
  Sn = n*t(beta_hat) %*% solve(cov_beta_hat) %*% beta_hat
  # Test for alpha1
  # LRT
  Tn = (n-df) * (RSS02-RSS12)/RSS12
  
  p_beta = 1 - pchisq(Sn,1)
  p_alpha1 = 1 - pchisq(Tn,1)
  
  std_alpha1 = sqrt(var_alpha1_hat/n)
  std_beta =  sqrt(cov_beta_hat/n)
  options("digits" = 4)
  result.df = data.frame(Coeffcient = c("$\\mathbf{\\alpha}_1$", "$\\mathbf{\\beta}$"),
                         Estimated_Coeffcient = c(alpha1_hat, beta_hat),
                         SE =c(std_alpha1,std_beta),
                         Test_statistics = c(Tn, Sn),
                         p_value = c(p_alpha1, p_beta))
  
  return(list(Sn = Sn, Tn = Tn, sigma1_hat = sigma1_hat,
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B, 
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat,
              p_beta = p_beta, p_alpha1 = p_alpha1,
              summary_result = result.df))
}


HBIC_calc <- function(lamb, xx,yy,mm,S = NULL, n_imp =8){
  # Calculate HBIC for a specific tuning parameter lambda
  #' @param lamb a float value of tuning parameter lambda 
  #' @param xx The n by q exposure matrix. q can be 1, and q < n is required
  #' @param yy The n-dimensional outcome vector.
  #' @param mm The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' @param n_imp an int for important mediators that will not be penalized
  
  #' @return
  #'        A list of class `HBIC_calc'
  #'        - BIC: HBIC score
  #'        - alpha0: estimated alpha0,
  #'        - alpha1: estimated alpha1,
  #'        - alpha2: estimated alpha2,
  #'        - sigma1_hat: estimated sigma1
  #'         
  n = nrow(xx)
  p = ncol(mm)
  q = ncol(xx)
  if(is.null(S)){
    s = 0
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q, n_imp = n_imp)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = NULL
    tmp = yy - mm%*%alpha0 - xx%*% alpha1
  }else{
    s = ncol(S)
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q,n_imp = n_imp, S = S)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = result[(p+q+1):(p+q+s)]
    
    tmp = yy - mm%*%alpha0 - xx%*% alpha1 - S %*% alpha2
  }
  
  df = length(which(alpha0!= 0))+q + s
  sigma_hat = t(tmp)%*%tmp/n
  BIC = log(sigma_hat) + df*log(log(n))*log(p+q + s)/n
  #obj = objective(xx,yy,M,alpha0,alpha1,lamb)
  return(list(BIC=BIC,alpha0=alpha0,alpha1 = alpha1, 
              alpha2 = alpha2, sigma1_hat = sigma_hat))
}

LLA_h1<- function(X,Y,M,lamb,n,p,q, n_imp = 0,S = NULL){
  #' Local linear approximation algorithm
  #' @param X matrix or dataframe, n by q exposure matrix. q can be 1, and requires q < n
  #' @param M matrix or dataframe, n by p mediator matrix and q < n.
  #' @param lamb float, tuning parameter lambda for the penalty
  #' @param n int, the sample size
  #' @param p int, X matrix dimension (number of columns) 
  #' @param q int, M matrix dimension (number of columns)
  #' @param n_imp int, for important mediators that will not be penalized
  
  if(length(S) == 0){
    s = 0
    V = X
  }else{
    s = ncol(S)
    V = cbind(X,S)
  }
  # Step 1 using Lasso
  w = matrix(0,nrow = (p + q + s),ncol=1)
  w[1:(p-n_imp)] = 1
  alpha_int = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                          lambda = lamb, penalty.factor=w,intercept = FALSE))[-1]
  # Step 2 using local linear approximation of SCAD
  w = matrix(0,nrow = (p+q + s),ncol=1)
  for(j in 1:(p- n_imp)){
    w[j] = deSCAD(alpha_int[j],lamb)
  }
  alpha = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                      lambda = lamb, penalty.factor=w, intercept = FALSE))
  return(alpha[-1])
}

deSCAD <- function(z,lamb,a=3.7){
  # First order derivative of SCAD penalty
  # tuning parameter "a" (or "gamma") use the default value 3.7
  return(1*(z<=lamb)+pmax((a*lamb-z),0)/((a-1)*lamb)*(lamb<z))
}