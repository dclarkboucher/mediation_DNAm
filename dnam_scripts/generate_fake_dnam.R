# The purpose of this script is to generate fake DNAm data resembling the 
# observed data from our study. That data can then be used in our analysis 
# pipeline as a basic example. We use the same variable names so that our
# analysis can be understood by those who have access to data from MESA

set.seed(1)
n <- 963 # samples with DNAm data who were not on diabetes meds (diabins5 == 0)
total_cpgs <- 2500 # actual number was 402,339

# Simulate low education exposure. 
# 1 is low education (below a 4 year degree), 0 is 4 year degree or more
educ1 <- sample(c(1,0), n, replace = T)

# Simulate genotyping information
chip_meth <- sample(108, n, replace = T) * (8.2 - 5.8) * 1e9 / 108 + 5.8e9
    # Methylation chip (values chosen roughly based on observed variable.)

pos_meth <- as.factor(paste0("c", sample(6, n, replace = T))) # methyl position

# Simulate covariates
age5c <- as.numeric(scale(rnorm(n)))  # Standardized age variable
bcell <- as.numeric(scale(rnorm(n)))  # b cells
tcell <- as.numeric(scale(rnorm(n)))  # t cells
nkcell <- as.numeric(scale(rnorm(n)))  # natural killer cells
neutro <- as.numeric(scale(rnorm(n)))  # neutrophils
race1c <- as.factor(paste0("r", sample(c(1, 3, 4), n, replace = T))) # race
female <- sample(2, n, replace = T) - 1

# Simulate DNAm
cpgs <- mvtnorm::rmvnorm(n, sigma = diag(nrow = total_cpgs))
alpha_a <- rnorm(total_cpgs, 0, 0.5) # alpha_a effects
cpgs <- cpgs + educ1 %*% t(alpha_a)
colnames(cpgs) <- paste0("cg", seq_len(total_cpgs))

# Simulate Y
beta_m <- rnorm(total_cpgs, 0, 1)
beta_a <- 1
hba1c <- as.numeric(scale((cpgs %*% beta_m + beta_a * educ1) + rnorm(n, 0, 2)))

mesa_sim <- data.frame(educ1, hba1c, chip_meth, pos_meth, age5c, bcell,
                       tcell, nkcell, neutro, race1c, female, cpgs)

saveRDS(mesa_sim, "dnam_files/mesa_data_sim.rds")




