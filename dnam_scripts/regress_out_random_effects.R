# Script to regress out random effects from mediators

# Load libraries
library(bigreadr)
library(lme4)
library(lmerTest)

# Read in data
mesa <- readRDS("dnam_files/mesa_data_sim.rds")

# Read list of desired CpGs and filter
cpg_list <- read.csv("dnam_files/want_cpgs.csv", header = T)$cpg

# Make objects used in function
variables <- c('chip_meth','pos_meth')
formula_right <- '(1|chip_meth) + (1|pos_meth)'

# Function to get model residuals
mediator_model <- function(cpg) {
  temp <- mesa[, c(cpg, variables)]

  # random effect model
  form <- as.formula(paste(cpg, ' ~ ', formula_right, sep=''))
  
  resids <-
    summary(
      lmer(form, data = temp, control = lmerControl(optimizer = "Nelder_Mead"))
    )$residuals
  
  output <- data.frame(cpg = resids)
  colnames(output) <- cpg
  
  output
  
}

# Loop through CpGs to get residuals
mylist <- list()
for (cpg in cpg_list){
  mylist[[cpg]] <- mediator_model(cpg)
}

outdat <- do.call(cbind, mylist)

# Write output
path <- paste0("dnam_files/dnam_use.rds")
saveRDS(outdat, file = path)



