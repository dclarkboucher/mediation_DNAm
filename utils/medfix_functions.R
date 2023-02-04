#MedFix script

library(gcdnet)


adlasso <- function(X,Y){

  Y <- as.vector(Y)

  colnames_original <- colnames(X)

  colnames(X) <- paste0("M",1:ncol(X))

  #initial fit with elastic net
  n <- nrow(X)

  lambda2_list <- exp(seq(log(0.0001), log(0.02), length.out = 50))

  nlambda2 <- length(lambda2_list)
  mincv <- numeric(nlambda2)

  message("   Obtaining initial fit\n")
  for(i in 1:nlambda2){

    cv_enet <-
      cv.gcdnet(
        X,
        Y,
        lambda2 = lambda2_list[i],
        standardize = T,
        method = "ls")

    if(i %% 10 == 0){
      message(paste0("     ",round(100 * i / (nlambda2  + 1)),"%\n"))
    }


    mincv[i] <- min(cv_enet$cvm)

  }

  lambda2_use <- lambda2_list[which.min(mincv)]


  mod1 <- cv.gcdnet(X,Y, lambda2 = lambda2_use, standardize = T, method = "ls")
  message(paste0("     ",100,"%\n"))


  #weight computation and second fit
  coefs <- coef(mod1, s = "lambda.min")
  v_star <- log(sum(coefs!=0)) / log(n)
  gamma <- ceiling(2*v_star/(1 - v_star)) + 1
  X_sd <- apply(X, 2, sd)
  weights_ad <- (abs(coefs[-1] * X_sd) + 1/n)^(-gamma)
  X1 <- X[, which(coefs[-1]!=0)]
  weights_ad_non0 <- weights_ad[which(coefs[-1] != 0)]
  message("   Obtaining weighted fit\n")
  mod2 <- cv.gcdnet(X1,Y, standardize = T, pf = weights_ad_non0,
                    method = "ls")

  #get p-values
  coefs1 <- coef(mod2, s = "lambda.min")
  y_preds <- predict(mod2,X1)
  s2 <- crossprod(Y - y_preds) / (n - sum(coefs1!=0) - 1)
  X2 <- X1[, which(coefs1[-1] != 0)]
  variance <- solve(crossprod(X2))*c(s2)

  #variances
  var1 <-
    diag(1 / apply(X2, 2, sd)) %*%
    variance %*%
    diag(1 / apply(X2, 2, sd))

  #variance for the intercept
  vari <-
    apply(X2, 2, mean) %*%
    variance %*%
    apply(X2, 2, mean)

  adlasso_tab <- matrix(0, sum(coefs1!=0), 3)
  rownames(adlasso_tab) <- rownames(coefs1)[which(coefs1!=0)]
  colnames(adlasso_tab) <- c("beta", "sd", "p-value")
  adlasso_tab[, 1] <- coefs1[which(coefs1!=0), ] #beta
  adlasso_tab[, 2] <- c(sqrt(vari), sqrt(diag(var1))) #sd
  adlasso_tab[, 3] <- sapply(abs(adlasso_tab[, 1]) / adlasso_tab[, 2],
                             function(x) 2*(1 - pnorm(x)))

  adlasso_tab1 <- data.frame(variable = rownames(adlasso_tab),adlasso_tab,
                             index = c(0,which(colnames(X) %in% rownames(adlasso_tab))),
                             row.names = NULL)

  #add in variables that were not selected
  zeros <- setdiff(colnames(X),adlasso_tab1$variable)
  df_zeros <- data.frame(variable = zeros, beta = 0, sd = 0 , `p-value` = 1,
                         index = which(colnames(X) %in% zeros))
  outdat <- rbind(adlasso_tab1,df_zeros)
  outdat <- outdat[order(outdat$index),][-1,]

  outdat$variable <- colnames_original

  rownames(outdat) <- NULL

  return(outdat)

}

screen_ym <- function(M,Y,A,C = rep(1,nrow(M)),p){
  
  M <- as.data.frame(M)
  
  pvs <- apply(M,2,
               function(x,Y,A,C){
                 return(summary(lm(Y ~ .,as.data.frame(cbind(x,A,C))))$coefficients[2,4])
                 },
               Y,A,C
               )
  
  
  #top p mediators
  keep <- rep(0,ncol(M))
  want_m <- order(pvs)[1:p]
  keep[want_m] <- 1
  names(keep) <- colnames(M)
  
  return(keep)
}


screen_ma <- function(A,M, C = rep(1,nrow(M)),p){
  
  M <- as.data.frame(M)
  pvs <- apply(M,2,
               function(x,A,C){
                 return(summary(lm(x~.,as.data.frame(cbind(A,C))))$coefficients[2,4]
                        )
               },
               A,C)
  
  #top p mediators
  keep <- rep(0,ncol(M))
  keep[order(pvs)[1:p]] <- 1
  names(keep) <- colnames(M)
  
  return(keep)
}

  
  
medfix <- function(Y,M,A,C1 = rep(1,nrow(M)),C2 = rep(1,nrow(M)),
                   p = ceiling(nrow(M) / log(nrow(M))), screen = c("ym","ma","none")){
  
  if(is.null(colnames(M))){
    colnames(M) <- paste0("M",1:ncol(M))
  }
  
  #Step 1: Pre-Screening
  screen <- screen[1]
  if(screen == "ym"){
    
    message("Screening based on outcome model")
    want_vars <- screen_ym(M,Y,A,C1,p=p)
    
    
  }else if(screen == "ma"){
    
    message("Screening based on mediator models")
    want_vars <- screen_ma(A,M1,C2,p=p)
    
  }else{
    
    want_vars <- rep(1,ncol(M))
    
  }
  
  M1 <- M[,as.logical(want_vars)]
  
  #Step 2: Fit outcome model with adaptive LASSO 
  
  ##2.a: Regress out A and C1 from Y and M
  AC1 <- as.data.frame(cbind(A,C1))
  Y_r <- lm(Y ~ ., AC1)$residuals
  M_r <- apply(M1,2,function(x){return(lm(x ~ ., AC1)$residuals)})
  
  sd_Y_r <- sd(Y_r)
  sd_M_r <- apply(M_r,2,sd)
  Y_r <- as.numeric(scale(Y_r))
  M_r <- as.matrix(scale(M_r))
  
  ##2.b: Adaptive lasso
  message("Fitting adaptive LASSO...")
  mod_out <- adlasso(M_r,Y_r)
  ym_dat <- mod_out[, c(1,2,4)]
  colnames(ym_dat) <- c("variable","beta","beta_pv")
  ym_dat$beta <- ym_dat$beta * sd_Y_r / sd_M_r #original order is maintained
  
  ym_dat <- ym_dat[ym_dat$beta != 0, ]
  
  #Step 3: Fit mediator models
  M_want <- M1[,ym_dat$variable]
  
  #get mediator model results
  message("Fitting mediator models")
  ma_results <- 
    lapply(as.data.frame(M_want),
           function(x,A,C){
             return(
               summary(lm(x~.,as.data.frame(cbind(A,C))))$coefficients[2,c(1,4)]
             )
           },
           A,C2)
  ma_results <- do.call(rbind,ma_results)
  
  ma_dat <- data.frame(variable = rownames(ma_results),
                       ma_results, row.names = NULL)
  colnames(ma_dat) <- c("variable","alpha","alpha_pv")
  
  dmerge <- merge(ym_dat,ma_dat, by = "variable", sort = F)
  dmerge$alpha_beta <- with(dmerge,alpha*beta)
  dmerge$pv <- pmax(dmerge$alpha_pv,dmerge$beta_pv)
  dmerge$pv_bonf <- pmin(dmerge$pv * nrow(dmerge), 1)
  dmerge$pv_fdr <- pmax(p.adjust(dmerge$alpha_pv, method = "fdr"),
                        p.adjust(dmerge$beta_pv, method = "fdr"))
  
  dmerge$gamma <- summary(lm(Y ~ ., data = AC1))$coefficients[2,1] # total effect
  
  return(dmerge[,-c(3,5)])
  
}




