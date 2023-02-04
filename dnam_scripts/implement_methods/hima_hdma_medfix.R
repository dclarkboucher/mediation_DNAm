# Script for running HIMA, HDMA, and MedFix on MESA data

# Load libraries
library(dplyr)
library(readr)
library(gcdnet)
library(lme4)
library(lmerTest)

source("utils/hima_hdma_utils.R")
source("utils/medfix_functions.R")

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

# Select desired mediators (SIS screening)
n <- nrow(cpg_dat)
p <- ceiling(n / log(n)) 
top_2k_mediators <- readRDS("dnam_files/single_mediator/single_med_results.rds")
top_cpgs <- with(top_2k_mediators, cpg[order(alpha_pv)][1:p])
cpg_dat <- cpg_dat[, top_cpgs]

# Identify covariates
cnames <- c("age5c","female","race1c","bcell","tcell","nkcell","neutro")
covariates <- dat |> dplyr::select(all_of(cnames))

# Regress out random effects from Y
y_formula <- hba1c ~ (1|chip_meth) + (1|pos_meth)
yv <- 
  summary(
    lmer(y_formula, data = dat, control = lmerControl(optimizer = "Nelder_Mead"))
  )$residuals

xv <- dat |> pull(educ1)
rm(dat)

# Zhang 2016 - HIMA -------------------------------------------------------

out_hima <-
  hima(
    X = xv,
    Y = yv,
    M = cpg_dat,
    COV.XM = covariates,
    COV.MY = covariates,
    topN = p,
    Y.family = "gaussian",
    M.family = "gaussian",
    penalty = "MCP",
    parallel = F
  )

saveRDS(out_hima, file = "dnam_files/multivariate_results/out_hima.rds")

# Gao 2019 - HDMA ---------------------------------------------------------
out_hdma <-
  hdma(
    X = xv,
    Y = yv,
    M = cpg_dat,
    COV.XM = covariates,
    COV.MY = covariates,
    family = "gaussian",
    topN = p,
    method = "lasso"
  ); out_hdma

saveRDS(out_hdma, file = "dnam_files/multivariate_results/out_hdma.rds")


# Zhang 2021 - MedFix -----------------------------------------------------
set.seed(1)
out_med <-
  medfix(
    Y = yv, 
    A = xv, 
    M = cpg_dat,
    C1 = covariates,
    C2 = covariates,
    p = p,
    screen = "none"
  ); out_med
saveRDS(out_med, file = "dnam_files/multivariate_results/out_med.rds")



