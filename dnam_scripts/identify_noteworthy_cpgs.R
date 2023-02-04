

# Identify Noteworthy CpGs ------------------------------------------------

library(tidyverse)
library(bama)
library(qvalue) # install.packages("BiocManager"); BiocManager::install("qvalue")

rm(list = ls())

results <- list()
l1 <- ls()

cdat <- read_csv("dnam_files/want_cpgs.csv", show_col_types = F)

# Univariate (One-at-a-time)
uni <- read_rds("dnam_files/single_mediator/single_med_results.rds")

results$uni <-
  uni |> 
  transmute(cpg = cpg, 
            alpha,
            beta,
            ab = alpha * beta,
            pv = qvalue(pmax(alpha_pv, beta_pv))$qvalues,
            noteworthy = pv < 0.1,
            method = "Univariate") |> 
  left_join(x = cdat, by = "cpg")

# BSLMM
load("dnam_files/multivariate_results/out_bslmm.rda")
ls()

results$bslmm <-
  out_bslmm |> 
  with( 
    tibble(
      method = "BSLMM",
      cpg = colnames(beta.m), 
      alpha = colMeans(alpha.a),
      beta = colMeans(beta.m),
      ab = colMeans(alpha.a * beta.m),
      pip = colMeans(r1 * r3),
      noteworthy = pip > 0
    )
  )

# HIMA
replace_na <- function(x,val = 0){ # helper function to replace omitted mediators
  logic <- is.logical(x)
  x[is.na(x)] <- val
  if(logic) return(as.logical(x))
  return(x)
}

hima <- read_rds("dnam_files/multivariate_results/out_hima.rds"); dim(hima)
results$hima <-
  hima |> 
  transmute(cpg = mediator,
            alpha,
            beta,
            ab = alpha_beta, 
            pv = ab_pv * n(),
            noteworthy = ab != 0) |>  # noteworthy?
  left_join(x = cdat, by = "cpg") |> 
  mutate(across(c(alpha,beta,ab,noteworthy),replace_na),
         method = "HIMA"
  )

# HDMA
hdma <- read_rds("dnam_files/multivariate_results/out_hdma.rds"); dim(hdma)

results$hdma <-
  hdma |> 
  transmute(cpg = mediator, alpha, beta, ab = alpha_beta, 
            pv = ab_pv * n(),
            noteworthy = ab != 0) |> 
  left_join(x = cdat, by = "cpg") |> 
  mutate(across(c(alpha,beta,ab,noteworthy),replace_na),
         method = "HDMA"
  )

# MedFix
med <- read_rds("dnam_files/multivariate_results/out_med.rds"); dim(med)
head(med)

results$med <-
  med |> 
  transmute(cpg = variable, alpha,beta,ab = alpha_beta, 
            pv = pv_bonf,
            noteworthy = ab != 0) |> 
  left_join(x = cdat, by = "cpg") |> 
  mutate(across(c(alpha,beta,ab,noteworthy),replace_na),
         method = "MedFix"
  )

#PLASSO
plasso <- readRDS("dnam_files/multivariate_results/out_plasso.rds") 
results$plasso <-
  plasso |> 
  mutate(
    noteworthy = ab != 0
  ) |> 
  left_join(x = cdat, by = "cpg") |> 
  mutate(across(c(alpha,beta,ab,noteworthy),replace_na),
         method = "Pathway LASSO"
  )

rm(list=setdiff(ls(),l1))


# Combine results
dat <- bind_rows(results)


# First task: All noteworthy CpGs (Supplementary file 1)
want_cpgs <- # noteworthy CpGs
  dat |> 
  group_by(cpg) |> 
  summarize(keep = any(noteworthy)) |> 
  filter(keep) |> 
  pull(cpg); length(want_cpgs)

out <-
  dat |> 
  mutate(noteworthy = ifelse(noteworthy,"Yes","No")) |> 
  pivot_wider(id_cols = cpg, names_from = c(method), 
              values_from = c(alpha,beta,ab,pv,pip,noteworthy),
              names_glue = "{method}_{.value}",
              names_var = "slowest"
  ) |>
  filter(apply(across(contains("noteworthy"), ~ .x == "Yes"),1,any)) |> 
  #very complex select. just trying to make sure the ordering is ok
  select(-contains("BSLMM"),starts_with("BSLMM"),-BSLMM_noteworthy,
         -BSLMM_pv,-contains("pip"),BSLMM_pip,BSLMM_noteworthy) |> 
  select(-`Pathway LASSO_pv`)

write_csv(out, file = "dnam_files/supp_file1_all_noteworthy_cpgs.csv")


# Second task: CpGs selected by >= 2 methods (Manuscript Table 2)
noteworthy2 <-
  dat |> 
  group_by(cpg) |> 
  summarize(n_noteworthy = sum(noteworthy)) |> 
  filter(n_noteworthy >= 2) |> 
  pull(cpg); length(noteworthy2)


number_noteworthy <- # number noteworthy, by method.
  out |> 
  select(contains("noteworthy")) |> 
  map_dbl(.f = ~sum(.x=="Yes"))

table2 <-
  out |> 
  filter(cpg %in% noteworthy2) |> 
  select(cpg,contains("ab"),contains("noteworthy")) |>
  arrange(desc(abs(Univariate_ab)),desc(abs(HIMA_ab)),desc(abs(HDMA_ab)),
          desc(abs(MedFix_ab)),desc(abs(`Pathway LASSO_ab`))) |> 
  pivot_longer(cols = contains("ab"), 
               names_to = "method", values_to = "ab") |> 
  pivot_longer(cols = contains("noteworthy"), 
               names_to = "method1", values_to = "noteworthy") |> 
  filter(method==str_replace(method1,"noteworthy","ab")) |> 
  mutate(
    ab = paste0(ifelse(ab==0,0,paste(round(100*ab,2),"x10-2")),
                ifelse(noteworthy=="Yes","*",""))
  ) |> 
  select(-noteworthy,-method1) |> 
  pivot_wider(id_cols = c(cpg), names_from = method, values_from = ab) |> 
  rename_with(
    ~ paste0(
      str_remove(.x,"_ab"),
      " (",
      number_noteworthy[str_replace(.x,"ab","noteworthy")],
      " sites identified)"
    ),
    .cols = contains("_ab")
  ) |> 
  rename(
    `CpG Site Name` = cpg
  )

write_csv(table2, file = "dnam_files/table2_cpgs.csv")









