#hdma extra functions

checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel & (ncore > 1)) {
    if (ncore > parallel::detectCores()) {
      message("You requested ", ncore, " cores. There are only ", 
              parallel::detectCores(), " in your machine!")
      ncore <- parallel::detectCores()
    }
    if (verbose) 
      message("    Running ", program.name, " with ", ncore, " cores in parallel...   (", 
              Sys.time(), ")")
    doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              Sys.time(), ")")
    registerDoSEQ()
  }
}



## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", 
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ", 
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text, 
                   "]};invisible(b)}")
  return(script)
}



## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}



## Internal function: scale data (obsolete function)

scaleto <- function(dat) {
  if (is.null(dat)) 
    return(list(dn = NULL, d = NULL, ds = NULL))
  dat_scale <- scale(dat)
  dat_names <- names(dat)
  if (any(class(dat) %in% c("matrix", "data.frame", "data.table"))) {
    dat_names <- colnames(dat)
    dat <- as.matrix(data.frame(dat_scale))
  } else {
    dat_names <- names(dat)
    dat <- as.numeric(dat_scale)
  }
  dat_scale <- as.numeric(attributes(dat_scale)[["scaled:scale"]])
  return(list(dn = dat_names, d = dat, ds = dat_scale))
}



# Internal function: Sure Independent Screening
# Global variables:
globalVariables("n")
globalVariables("M_chunk")

himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel, ncore, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }
  
  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                             glm.family, ", data = datarun))"), "c(1,4)")))
  
  checkParallel(tag, parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}
  
  colnames(results) <- M.names
  return(results)
}



# Internal function: LDPE
# the code of Liu han's JRSSB paper for high-dimensional Cox model
# ID: the index of interested parameter
# X: the covariates matrix with n by p
# OT: the observed time = min(T,C)
# status: the censoring indicator I(T <= C)

LDPE_func <- function(ID, X, OT, status){
  coi <- ID
  x <- X
  d <- dim(x)[2]
  n <- dim(x)[1]
  
  ##Set of tuning parameters
  PF <- matrix(1,1,d)
  PF[ID] <- 1
  fit <- glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  cv.fit <- cv.glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  betas   <-   coef(fit, s = cv.fit$lambda.min)[1:d]  # the semi-penalized initial estimator  # initial estimator
  
  stime = sort(OT)          # Sorted survival/censored times
  otime = order(OT)         # Order of time
  
  Vs  = matrix(rep(0,d*d),nrow = d)
  Hs  = Vs                                 # Hessian
  ind = 0
  
  la  = rep(0,n)                           # Gradient w.r.t parameter of interest
  lb  = matrix(rep(0,(d-1)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
  i   = 1
  
  while( i<=n)
  {
    if (status[otime[i]]==1)
    {
      ind = which(OT >= stime[i])
      S0  = 0
      S1  = rep(0,d)
      S2  = matrix(rep(0,d*d),nrow = d)
      
      if (length(ind)>0)
      {
        for (j in 1:length(ind))
        {
          tmp = exp(x[ind[j],]%*%betas)
          S0  = S0 + tmp
          
          S1  = S1 + tmp %*%t(x[ind[j],])
          
          tmp = apply(tmp,1,as.numeric)     
          S2  = S2 + tmp*x[ind[j],]%*%t(x[ind[j],])
        }
      }
      S0 = apply(S0,1,as.numeric)
      
      la[i]  = -(x[otime[i],coi] - S1[coi]/S0)
      if (coi == 1)
      {
        lb[i,] = -(x[otime[i],c((coi+1):d)] - S1[c((coi+1):d)]/S0)
      } else if (coi == d){
        lb[i,] = -(x[otime[i],c(1:(coi-1))] - S1[c(1:(coi-1))]/S0)
      } else {
        lb[i,] = -(x[otime[i],c(1:(coi-1), (coi+1):d)] - S1[c(1:(coi-1), (coi+1):d)]/S0)
      }
      V   = S0*S2 - t(S1)%*%(S1)
      Hs  = Hs + V/(n*S0^2)          
    }
    i = i + 1
  }
  
  fit <- glmnet(lb,la,alpha = 1, standardize = FALSE,intercept = FALSE,lambda = sqrt(log(d)/n))
  what <- as.numeric(stats::coef(fit)[2:d])   
  
  if (coi == 1)
  {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi]
  } else if (coi == d){
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi]
  } else {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi]
  }
  
  beta_est <- S
  beta_SE  <- sqrt(1/(n*var))
  
  result <- c(beta_est,beta_SE)
  
  return(result)
}


hdma <- function (X, Y, M, COV.XM = NULL, COV.MY = COV.XM, 
                  family = c("gaussian","binomial"), 
                  method = c("lasso", "ridge"), topN = NULL,
                  parallel = FALSE, ncore = 1, 
                  verbose = FALSE, ...){	
  
  
  ####################################################################################################################################
  #########################################                   Function body                ###########################################
  ####################################################################################################################################
  ####### INPUT
  ####### X : Independent variable that is a vector
  ####### Y : Dependent variable that is a vector and can be either continuous or binary variable
  ####### M : High-dimensional mediators that can be either data.frame or matrix. Rows represent samples, columns represent variables
  ####### COV.XM : a data.frame or matrix of covariates dataset for testing the association X ~ M. Default = NULL. 
  #######          If the covariates contain mixed types, please make sure all categorical variables are properly transformed into factor
  #######          type.
  ####### COV.MY : a data.frame or matrix of covariates dataset for testing the association Y ~ M. Using covariates should be careful.
  #######          If the cavariables are not specified, the covariates for Y ~ M are the same with that of M ~ X.
  ####### family : either 'gaussian' or 'binomial', relying on the type of outcome (Y). See hdi package.
  ####### method : either "lasso" or "ridge" to estimate the effect of M -> Y.
  ####### topN : an integer can be used to set the number of top markers by the method of sure independent screening. Default = NULL.
  #######        If topN is NULL, it will be either ceiling(n/log(n)) if family = 'gaussian', or ceiling(n/(2*log(n))) if family = 
  #######	      'binomial', where n is the sample size. If the sample size is greater than topN (pre-specified or calculated), all
  #######        markers will be harbored in the test.
  ####### parallel : logical parameter. The parameter can be employed to enable your computer to do parallel calculation. Default = FALSE.
  ####### ncore : the parameter can be used to set the number of cores to run parallel computing when parallel == TRUE. By default max
  ########	number of cores available in the machine will be utilized.
  ####### verbose : logical. Default = FALSE.
  ####### ... : other arguments passed to hdi.
  ####################################################################################################################################
  ####### Values 
  ####### alpha : the coefficient can reflect the association of X –> M, note that the effect is adjusted by covariables when covariables
  #######	are not NULL.
  ####### beta : the coefficient can reflect the association of M –> Y, note that the effect is adjusted by X. When covariables are not
  #######	NULL, the effect is adjusted by X and covariables. 
  ####### gamma : the coefficient can reflect the linkage of X –> Y, it can represent the total effect.
  ####### alpha*beta : the estimator of mediation effect.
  ####### %total effect : alpha*beta/gamma*100. The proportion of the mediation effect is out of the total effect.
  ####### p-values : joint significant test for mediators.
  ####################################################################################################################################
  ####### checking the necessary packages 
  
  
  pkgs <- list("hdi","MASS","doParallel", "foreach","iterators")
  checking<-unlist(lapply(pkgs, require, character.only = T))
  if(any(checking==F))
    stop("Please install the necessary packages first!")	
  family <- match.arg(family)
  method <- match.arg(method)
  if (parallel & (ncore == 1)) ncore <- parallel::detectCores()
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2*n/log(n)) 
  } else {
    d <- topN      # the number of top mediators that associated with independent variable (X)
  }
  d <- min(p, d)   # if d > p select all mediators
  #############################################################################################################################
  ################################           Step-1 Sure Independent Screening (SIS)          #################################
  #############################################################################################################################
  message("Step 1: Sure Independent Screening ...", " (", Sys.time(), ")")
  if(family == "binomial") 
  {
    if(verbose) message("Screening M using the association between X and M: ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", parallel = parallel, 
                                   ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", parallel = parallel,
                           ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)						 
  #################################################################################################################################
  ####################################          Step-2  High-dimensional Inference (HDI)         ##################################
  #################################################################################################################################
  message("Step 2: High-dimensional inference (", method, ") ...", "     (", Sys.time(), ")")
  ## Based on the SIS results in step 1. We will find the most influential M on Y.	
  if (is.null(COV.MY)) {
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM, Y, family = family) else fit <- ridge.proj(XM, Y, family = family)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM_COV, Y, family = family) else fit <- ridge.proj(XM_COV, Y, family = family)
  }
  P_hdi<-fit$pval[1:length(ID)]
  index<-which(P_hdi<=0.05)
  if(verbose)  message("Non-zero ",method, " beta estimate(s) of mediator(s) found: ", paste0(names(index), collapse = ","))
  if(length(index)==0)
    stop("No mediatior is identified !")
  ID_test <- ID[index]  
  if(family == "binomial")
  {
    ## This has been done in step 1 (when Y is binary)
    alpha_est <- alpha[,ID_test, drop = FALSE]
  } else {
    if(verbose) message(" Estimating alpha (effect of X on M): ", appendLF = FALSE)
    alpha_est <- himasis(NA, M[,ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                         parallel = FALSE, ncore = ncore, verbose, tag = "site-by-site ordinary least squares estimation")
  }   
  beta_P<-P_hdi[index]
  beta_hat<-fit$bhat[index]              # the estimator for beta
  alpha_hat<-as.numeric(alpha_est[1, ])
  ab_est<-beta_hat*alpha_hat
  alpha_P<-alpha_est[2,]
  PA <- rbind(beta_P,alpha_P)
  P_value <- apply(PA,2,max)
  ###################################################################################################################################    
  ###############################          STEP 3   Computing the proportion of mediation effect          ############################
  ###################################################################################################################################
  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  gamma_est <- stats::coef(glm(Y ~ ., family = family, data = YX))[2]
  results <- data.frame(alpha = alpha_hat, 
                        alpha_pv = alpha_P,
                        alpha_pv_fdr = p.adjust(alpha_P, method = "fdr"),
                        beta = beta_hat, 
                        beta_pv = beta_P,
                        beta_pv_fdr = p.adjust(beta_P, method = "fdr"),
                        gamma = gamma_est, 
                        `alpha*beta` = ab_est, 
                        `%total effect` = ab_est/gamma_est*100, 
                        `P.value` = P_value,
                        check.names = FALSE)
  results$ab_pv_fdr = with(results,pmax(beta_pv_fdr,alpha_pv_fdr))
  
  message("Done!", " (", Sys.time(), ")")
  doParallel::stopImplicitCluster()
  return(results)
}


