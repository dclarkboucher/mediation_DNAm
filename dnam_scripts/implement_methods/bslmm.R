
library(bama)
library(lme4)
library(dplyr)
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

# Identify covariates
cnames <- c("age5c","female","race1c","bcell","tcell","nkcell","neutro")
covariates <- dat |> select(all_of(cnames))

# Regress out random effects from Y
y_formula <- hba1c ~ (1|chip_meth) + (1|pos_meth)
yv <- 
  summary(
    lmer(y_formula, data = dat, control = lmerControl(optimizer = "Nelder_Mead"))
  )$residuals

xv <- dat |> pull(educ1)

# Load single-mediator results and get variances to use for priors
top2km <- readRDS("dnam_files/single_mediator/single_med_results.rds")
var_alpha1 <- var(top2km$alpha[abs(top2km$alpha) > 
                                 quantile(abs(top2km$alpha),0.9)])
var_beta1 <- var(top2km$beta[abs(top2km$beta) > 
                               quantile(abs(top2km$beta),0.9)])

out_bslmm <- 
  bama(
    Y = yv, 
    A = xv,
    C1 = model.matrix(~., covariates),
    C2 = model.matrix(~., covariates),
    M = as.matrix(cpg_dat),
    method = "BSLMM",
    burnin = 35000, 
    ndraws = 40000,
    control = 
      list(
        lm0 = 10^-4,
        lm1 = var_beta1,
        lma1 = var_alpha1),
    seed = 123
  ) 


save(out_bslmm, file = "dnam_files/multivariate_results/out_bslmm.rda")



