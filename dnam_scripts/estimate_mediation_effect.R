# Estimated global mediation effect (Manuscript Table 3)

#read hba1c output for CpGs

library(tidyverse)
library(bama)
rm(list = ls())

cdat <- read_csv("dnam_files/want_cpgs.csv", show_col_types = F)

results <- list()
l1 <- ls()

# BSLMM
load("dnam_files/multivariate_results/out_bslmm.rda")

results$bslmm <-
  out_bslmm |> 
  with( 
    tibble(
      method = "BSLMM",
      tie = mean(rowSums(alpha.a * beta.m)), 
      tie_lower = quantile(rowSums(alpha.a * beta.m),0.025),
      tie_upper = quantile(rowSums(alpha.a * beta.m),0.975),
      de = mean(beta.a), 
      de_lower = quantile(beta.a,0.025),
      de_upper = quantile(beta.a,0.975),
      te = mean(beta.a + rowSums(alpha.a * beta.m)),
      te_lower = quantile((beta.a + rowSums(alpha.a * beta.m)),0.025),
      te_upper = quantile((beta.a + rowSums(alpha.a * beta.m)),0.975),
      pmed = tie / te
    )
  ) 

# HIMA
hima <- read_rds("dnam_files/multivariate_results/out_hima.rds"); hima
results$hima <-
  hima |> 
  summarize(method = "HIMA", tie = sum(alpha_beta), te = te[1], 
            de = te - tie, pmed = tie / te)
  
# HDMA
hdma <- read_rds("dnam_files/multivariate_results/out_hdma.rds"); head(hdma)

results$hdma <-
  hdma |> 
  summarize(method = "HDMA", tie = sum(alpha_beta), te = te[1],
            de = te - tie, pmed = tie / te)

# MedFix
med <- read_rds("dnam_files/multivariate_results/out_med.rds"); head(med)

results$med <-
  med |> 
  summarize(method = "MedFix", tie = sum(alpha_beta), te = gamma[1],
            de = te - tie, pmed = tie / te)

# Pathway LASSO
plasso <- readRDS("dnam_files/multivariate_results/out_plasso.rds") 
results$plasso <-
  plasso |> 
  summarize(method = "Pathway LASSO", tie = sum(ab), de = de[1],
            te = de + tie, pmed = tie / te)

# PCMA and SPCMA
load("dnam_files/multivariate_results/out_spcma.rda")
results$pcma <-
  out_spcma |> 
  with(
    tibble(method = "PCMA", 
           tie = PCA$IE.total[1,1], 
           tie_lower = PCA$IE.total[1,3],
           tie_upper = PCA$IE.total[1,4],
           tie_pv = PCA$IE.total[1,2],
           de = PCA$DE[1],
           de_lower = PCA$DE[3],
           de_upper = PCA$DE[4],
           te = de + tie,
           pmed = tie / te
           )
  )

results$spcma <-
  out_spcma |> 
  with(
    tibble(method = "SPCMA", 
           tie = SPCA$IE.total[1,1], 
           tie_lower = SPCA$IE.total[1,3],
           tie_upper = SPCA$IE.total[1,4],
           tie_pv = SPCA$IE.total[1,2],
           de = SPCA$DE[1],
           de_lower = SPCA$DE[3],
           de_upper = SPCA$DE[4],
           te = de + tie,
           pmed = tie / te
    )
  )

# HILMA
load("dnam_files/multivariate_results/out_hilma.rda")
results$hilma <-
  out_hilma |> 
  with(
    tibble(method = "HILMA", 
           tie = beta_hat, 
           tie_pv = pvalue_beta_hat,
           de = alpha1_hat,
           te = tie + de,
           pmed = tie / te
    )
  )

rm(list = setdiff(ls(),l1)) # delete extra stuff we loaded along the way

# Combine results
dat <- 
  bind_rows(results) |> 
  select(-contains("lower"),-contains("upper"),-contains("pv")) |> 
  mutate(method = factor(method, levels = c("HIMA","HDMA","MedFix","BSLMM",
                                            "Pathway LASSO",
                                            "PCMA","SPCMA","HILMA")),
         across(starts_with("tie"), ~ round(.x,3)),
         across(starts_with("de"), ~ round(.x,3)),
         across(starts_with("te"), ~ round(.x,3)),
         across(starts_with("pmed"), ~ round(.x,3))
         #interval = paste0("(",tie_lower,", ",tie_upper,")")
  ) |> 
  as.data.frame() |> 
  rename(
    `Global Mediation Effect` = tie,
    `Direct Effect` = de,
    `Total Effect` = te, 
    `Proportion Mediated` = pmed
  )


dat |> 
  arrange(method) |> 
  select(method,tie,de,te,pmed) |> 
  write.csv("dnam_files/table3_mediation_effects.csv")

