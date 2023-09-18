
library(lme4)
library(dplyr)
library(lmerTest)

source("~/utils/pmed_functions.R")

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

# Implement PMED

out <- NULL
try({ # Errors happen when no mediators are selected
  out <- pmed_wrapper(A = xv, M = as.matrix(cpg_dat), Y = yv, 
                      C = model.matrix(~ ., covariates), nlambda = 100)
})

if(!is.null(out)){
  tie <- out$summary[2,2]
  de <- out$summary[1,2]
  
}else{
  tie <- 0
  de <- coef(lm(yv ~ ., as.data.frame(cbind(xv, covariates))))[2]
}

pmed_tie <- tie
pmed_de <- de
save(pmed_tie, pmed_de, file = "dnam_files/multivariate_results/out_pmed.rda")


