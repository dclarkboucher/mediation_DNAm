# Script for running PCMA and SPCMA on MESA data

# Load libraries
library(spcma) # devtools::install_github("https://github.com/zhaoyi1026/spcma")
library(dplyr)
library(readr)
library(stringr)
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

# Note that since SPCMA uses the fused LASSO penalty, which attempts to 
# give adjacent variables similar regression coefficients, it is preferable
# to order the CpGs by their position in the genome. One can do this very easily
# by sorting the CpGs in 'cpg_dat' by chromosome and position. Since we have 
# not assigned chromosomes and positions to our fake CpG data, we will simply pretend
# that they are numbered by position and sort by name. This is specific to SPCMA
# and for the other methods the order of the variables should not matter. 

cpgs <- colnames(cpg_dat)
cpgs <- cpgs[order(as.numeric(str_remove(cpgs,"cg")))]
cpg_dat <- cpg_dat[, cpgs]

# At last we can apply SPCMA.
out_spcma <- 
  spcma(
    X = xv, 
    M = cpg_dat,
    Y = yv, 
    adaptive = F, 
    n.pc = 100,
    gamma = 2,
    sims = 100
  )

save(out_spcma, file = "dnam_files/multivariate_results/out_spcma.rda")

