# Script for running HDMM on MESA data

# Load libraries
library(PDM) # devtools::install_github("https://github.com/oliverychen/PDM")
  #please make sure you have the correct PDM package loaded. Do not use the one on
  #CRAN, which is unrelated. 
library(mediation) # install.packages("mediation")
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)

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
covariates <- dat |> dplyr::select(all_of(cnames))

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

# Use PDM_1 function to get the first directions of mediation. 
set.seed(123)
hdmm_out <- 
  PDM_1(
    x = xv, 
    y = yv, 
    m = cpg_dat,
    theta = rep(1, 5),
    w1 = rep(1, p),
    interval = 10^6, 
    step = 10^4,
    imax = 100,
    tol = 10^-5
  )

# Use R package "mediation" to perform mediation analysis with latent mediator
est_w <- hdmm_out$w1
dm1 <- cpg_dat %*% as.matrix(est_w, ncol = 1)
m <- lm(dm1 ~ xv) 
y <- lm(yv ~ xv + dm1)

set.seed(633788)
med <- mediate(m, y, sims = 2000, treat = "xv", mediator = "dm1")

outdat <- data.frame(cpg = top_cpgs, wt = est_w)
rownames(outdat) <- NULL

out_hdmm <- c()
out_hdmm$pdms <- outdat
out_hdmm$med <- med
save(out_hdmm, file = "dnam_files/multivariate_results/out_hdmm.rda")
