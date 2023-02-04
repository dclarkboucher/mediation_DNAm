# Script for running HILMA on MESA data

# Load libraries
library(freebird) #install.packages("freebird")
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

# Run HILMA
set.seed(123) 
# Hilma performs slightly better, in our experience, if the mediators are
# standardized prior to analysis. Since mediation pathways (and therefore
# the global mediation effect) are independent of the scale of the mediators,
# this does not change the interpretation of our results
out_hilma <-  hilma(yv, scale(cpg_dat), matrix(xv, ncol = 1), center = T)

save(out_hilma, file = "dnam_files/multivariate_results/out_hilma.rda")



