# Script to fit single-mediator models for initial mediator screening

# With a large number of mediators this script will need to be parallelized.
# In particular, to break this script into multiple jobs, define index to be 
# the index of each job and set ncores to the total number of jobs. The example
# here is simple so we pretend only 1 job is needed. 

rm(list=ls())
ncores <- 1
index <- 1

# Load libraries
library(lme4)
library(lmerTest)

# Read in data
mesa <- readRDS("dnam_files/mesa_data_sim.rds")
cpg_list <- colnames(mesa)[grepl("cg", colnames(mesa))]
n_cpgs <- length(cpg_list)
cpgs_per_job <- ceiling(n_cpgs / ncores)

# Make start and end indices. 
start <- (index - 1) * cpgs_per_job 
end <- min(start + cpgs_per_job, n_cpgs)

# Make objects used in function
variables <- c('educ1','age5c','female','race1c','bcell','tcell',
               'nkcell','neutro','chip_meth','pos_meth')
formula_right <- 
  'educ1 + age5c + female + race1c + bcell + tcell + nkcell + neutro + (1|chip_meth)+(1|pos_meth)'

# Standardize things (though on example data we already did)
scale <- function(x) as.numeric(base::scale(x))
mesa$age5c <- scale(mesa$age5c)
mesa$neutro <- scale(mesa$neutro)
mesa$bcell <- scale(mesa$bcell)
mesa$tcell <- scale(mesa$tcell)
mesa$nkcell <- scale(mesa$nkcell)
mesa$hba1c <- scale(mesa$hba1c)


# Mediator models ---------------------------------------------------------

# Function to fit M ~ A + C model
mediator_model <- function(cpg) {
  temp <- mesa[, c(cpg, variables)]
  
  # Standardize cpg
  temp[[cpg]] = scale(temp[[cpg]])

  # random effect model
  form <- as.formula(paste0(cpg, ' ~ ', formula_right))
  
  lmer_result <-
    summary(
      lmer(form, data = temp, control = lmerControl(optimizer = "Nelder_Mead"))
    )
  
  # Output result
  unname(lmer_result$coefficients['educ1', c("Estimate", "Std. Error", "Pr(>|t|)")])

}

# Make empty dataset for results
mediator_results <- as.data.frame(matrix(NA, end - start, 4))
colnames(mediator_results) = c('cpg', 'alpha', 'alpha_se', 'alpha_pv')
mediator_results$cpg <- cpg_list[(start + 1):end]

# Loop through CpG chunk
for (i in 1:(end - start)) {
  cpg <- as.character(cpg_list[i + start])
  mediator_results[i, 2:4] <- mediator_model(cpg)
} 

# Write output
path <- paste0("dnam_files/single_mediator/ma_out", index, ".csv")
write.csv(mediator_results, file = path, quote = F, row.names = F)
rm(mediator_results)

# Outcome models ----------------------------------------------------------


# Function to fit Y ~ M + A + C model
outcome_model <- function(cpg) {
  temp <- mesa[, c(cpg, variables, "hba1c")]
  
  # Standardize cpg
  temp[[cpg]] = scale(temp[[cpg]])
  
  # random effect model
  form <- as.formula(paste("hba1c ~", cpg, "+", formula_right))
  
  lmer_result <-
    summary(
      lmer(form, data = temp, control = lmerControl(optimizer = "Nelder_Mead"))
    )
  
  # Output result
  unname(lmer_result$coefficients[cpg, c("Estimate", "Std. Error", "Pr(>|t|)")])
  
}


# Make empty dataset for results
mediator_results <- as.data.frame(matrix(NA, end - start, 4))
colnames(mediator_results) = c('cpg', 'beta', 'beta_se', 'beta_pv')
mediator_results$cpg <- cpg_list[(start + 1):end]

# Loop through CpG chunk
for (i in 1:(end - start)) {
  cpg <- as.character(cpg_list[i + start])
  mediator_results[i, 2:4] <- outcome_model(cpg)
} 

# Write output
path <- paste0("dnam_files/single_mediator/ym_out", index, ".csv")
write.csv(mediator_results, file = path, quote = F, row.names = F)




