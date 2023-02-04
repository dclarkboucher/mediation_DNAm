# Script for implementing pathway lasso on DNAm data

# Load libraries
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
source("utils/plasso_functions.R")

# Load MESA data
dat <- readRDS("dnam_files/mesa_data_sim.rds")

# Standardize things and simplify dataset
dat <- 
  dat |> 
  transmute(
    race1c = race1c,
    female = female,
    chip_meth = chip_meth,
    pos_meth = pos_meth,
    educ1 = educ1,
    age5c = as.numeric(scale(age5c)),
    bcell = as.numeric(scale(bcell)),
    tcell = as.numeric(scale(tcell)),
    nkcell = as.numeric(scale(nkcell)),
    neutro = as.numeric(scale(neutro)),
    hba1c = as.numeric(scale(hba1c))
  )

# Load mediators with random effects regressed out and scale them
cpg_dat <- readRDS("dnam_files/dnam_use.rds")
cpg_dat <- scale(cpg_dat)

# Subset mediators further. This is necessary because applying HDMM when there
# are more mediators than observations is much more challenging and requires
# lengthy, time-consuming computations. See the original HDMM paper by Chen
# et. al for more details 
n <- nrow(cpg_dat)
p <- ceiling(n / log(n)) # same p used in HIMA, HDMA
top_2k_mediators <- readRDS("dnam_files/single_mediator/single_med_results.rds")
top_cpgs <- with(top_2k_mediators, cpg[order(alpha_pv)][1:p])
cpg_dat <- cpg_dat[, top_cpgs]

# Identify covariates
cnames <- c("age5c","female","race1c","bcell","tcell","nkcell","neutro")
covariates <- dat |> select(all_of(cnames))

# Regress covariates out of M since methods cannot handle them directly
cpg_dat <- apply(cpg_dat, 2, function(x){lm(x ~ ., covariates)$residuals})

# Regress out covariates and random effects from Y
y_formula <- 
  as.formula(
    paste0("hba1c ~", paste0(cnames, collapse = " + "), " + (1|chip_meth)+(1|pos_meth)")
  )

# Create final X variable and Y variable
yv <- 
  summary(
    lmer(y_formula, data = dat, control = lmerControl(optimizer = "Nelder_Mead"))
  )$residuals
xv <- dat |> pull(educ1)
rm(dat)


# Implement Pathway LASSO -------------------------------------------------

# Change to pathway lasso naming scheme 
Z <- xv
R <- yv
M <- cpg_dat
k <- ncol(M)
dd0 <- data.frame(Z = Z, M, R = R)
Sigma10 <- diag(rep(1, k))
Sigma20 <- matrix(1, 1, 1)

# Record SDs and standardize
sd.Z <- sd(Z)
sd.M <- apply(M, 2, sd)
sd.R <- sd(R)

Z <- scale(Z)
M <- scale(M)
R <- scale(R)
dd <- data.frame(Z = Z, M, R = R)


# Set candidate lambdas
lambdas <-
  c(10 ^ c(
    seq(-5,-3, length.out = 5),
    seq(-3, 0, length.out = 21)[-1],
    seq(0, 2, length.out = 11)[-1],
    seq(2, 4, length.out = 6)[-1]
  ))

nlambda <- length(lambdas)

# Define a few parameters
rho <- 1
thred <- 1e-6
thred2 <- 1e-3
zero.cutoff <- 1e-3
omega_ratio <- 1
phi <- 2
maxit <- 5000
tol <- 1e-6

# Output objects to fill in
re <- vector("list", length = length(lambdas))
AB.est = A.est = B.est <- matrix(NA, k, length(lambdas))
C.est <- rep(NA, length(lambdas))

# Fit pathway lasso for each lambdas
for (i in nlambda:1)
{
  out <- NULL
  if (i == nlambda)
  {
    # starting from the largest lambda value
    try(out <- mediation_net_ADMM_NC(
      Z,
      M,
      R,
      lambda = lambdas[i],
      omega = omega_ratio * lambdas[i],
      phi = phi,
      Phi1 = NULL,
      Phi2 = NULL,
      rho = rho,
      rho.increase = FALSE,
      tol = tol,
      max.itr = maxit,
      thred = thred,
      Sigma1 = Sigma10,
      Sigma2 = Sigma20,
      trace = FALSE
    ))
  } else
  {
    # for smaller lambda (ith lambda), use the (i+1)th lambda results as burn-in
    try(out <- mediation_net_ADMM_NC(
      Z,
      M,
      R,
      lambda = lambdas[i],
      omega = omega_ratio * lambdas[i],
      phi = phi,
      Phi1 = NULL,
      Phi2 = NULL,
      rho = rho,
      rho.increase = FALSE,
      tol = tol,
      max.itr = maxit,
      thred = thred,
      Sigma1 = Sigma10,
      Sigma2 = Sigma20,
      trace = FALSE,
      Theta0 = matrix(c(1, A.est[, i + 1] *
                          (sd.Z / sd.M)), nrow = 1),
      D0 = matrix(c(
        C.est[i + 1] * (sd.Z / sd.R), B.est[, i + 1] * (sd.M / sd.R)
      ), ncol = 1),
      alpha0 = matrix(c(1, A.est[, i + 1] *
                          (sd.Z / sd.M)), nrow = 1),
      beta0 = matrix(c(
        C.est[i + 1] * (sd.Z / sd.R), B.est[, i + 1] * (sd.M / sd.R)
      ), ncol = 1)
    ))
  }
  
  if (!is.null(out))
  {
    re[[i]] <- out
    
    # scale the estimate back to the original scale
    B.est[, i] <- out$B * (sd.R / sd.M)
    C.est[i] <- out$C * (sd.R / sd.Z)
    A.est[, i] <- out$A * (sd.M / sd.Z)
    AB.est[, i] <- A.est[, i] * B.est[, i]
  }
  message(paste("Lambda", i, "finished"))
}

plasso_varnames
# save(
#   list = c("re","AB.est","A.est","B.est","C.est"),
#   file = "dnam_files/multivariate_results/plasso_fit.rda"
#   )

# Now apply Variable Selection Stability Criterion to choose lambda
# Note: this might take a while. 
vss_rep <- 5 # Number of repetitions
vss_cutoff <- 0.1
vss_results <-
  mediation_net_ADMM_NC_KSC(Z,M,R,zero.cutoff=zero.cutoff,n.rep=vss_rep,
                            vss.cut=vss_cutoff,lambda=lambdas,
                            omega=omega_ratio*lambdas,
                            phi=phi,Phi1=NULL,Phi2=NULL,rho=rho,
                            rho.increase=FALSE,tol=tol,max.itr=maxit,
                            thred=thred, Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,
                            Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)
# save(vss_results, file = "dnam_files/multivariate_results/plasso_vssc.rda")
which_lambda <- vss_results$lambda.idx[1]

out <- 
  data.frame(
    cpg = colnames(cpg_dat),
    alpha = A.est[,which_lambda],
    beta = B.est[,which_lambda],
    ab = AB.est[,which_lambda],
    de = C.est[which_lambda]
  )

saveRDS(out, file = "dnam_files/multivariate_results/out_plasso.rds")



