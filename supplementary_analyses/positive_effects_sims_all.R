library(tidyverse)
library(spcma)
library(bama)
library(freebird)

source("utils/hima_hdma_utils.R")
source("utils/bslmm_function.R")
source("utils/medfix_functions.R")
source("utils/plasso_functions.R")
source("utils/univariate_function.R")
source("utils/pmed_functions.R")
source("utils/get_data_positive_effects.R")

out_loc <- "supplementary_analyses/positive_effects_results/"

filename <- function(method, seed, dataset){
  
  if (method %in% c("hilma", "pcma", "plasso", "pmed")){
    app <- ".rda"
  } else{
    app <- ".rds"
  }
  
  paste0(out_loc,method,"_s",seed,"_d",dataset, app)
  
}

tab <- expand_grid(seed = seq_len(100), dataset = 1:4)

# index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
index <- 1 # Job index
tasks <- 400 # Number of tasks (i.e., number of rows of "tab" data frame)
cores <- 1 # Number of parallel jobs
tasks_here <- rev(seq(index, tasks, cores)) # Tasks to do in this job
tasks_here <- 1
for(t in tasks_here){
  # print(t)
  dataset <- tab$dataset[t]
  seed <- tab$seed[t]
  
  dat <-
    get_data_pos(
      seed2 = seed,
      dataset = dataset
    )

  Y <- dat$Y
  A <- dat$A
  M <- dat$M

  p <- ceiling(nrow(M) / log(nrow(M)))

  # HIMA
  if(!file.exists(filename("hima",seed, dataset))){
    set.seed(123)
  out_hima <-
    hima(
      X = as.vector(A),
      Y = as.vector(Y),
      M = as.data.frame(M),
      topN = p,
      family = "gaussian",
      penalty = "MCP",
      max.iter = 10^5
    )

    saveRDS(out_hima, file = filename("hima",seed, dataset))

  }

  # MedFix
  if (!file.exists(filename("med", seed, dataset))){

    dat <-
      get_data_pos(
        seed2 = seed,
        dataset = dataset
      )

    Y <- dat$Y
    A <- dat$A
    M <- dat$M

    p <- ceiling(nrow(M) / log(nrow(M)))



    set.seed(123)
    out_med <- medfix(Y,M,A,p=p,screen = "ym")
    saveRDS(out_med, file = filename("med", seed, dataset))
  }

  # HDMA
  if (!file.exists(filename("hdma", seed, dataset))){
    set.seed(123)
    out_hdma <-
      hdma(
        X = as.vector(A),
        Y = as.vector(Y),
        M = as.data.frame(M),
        family = "gaussian",
        topN = p,
        method = "lasso",
      )
    saveRDS(out_hdma, file = filename("hdma", seed, dataset))

  }

  # PCMA
  if (!file.exists(filename("pcma", seed, dataset))){
    set.seed(123)
    try(out_pcma <-
          try(mcma_PCA(
            X = as.vector(A),
            M = as.matrix(M),
            Y = as.vector(Y),
            adaptive = F,
            n.pc = 100,
            boot = TRUE,
            sims = 2,
          )))
    try(save(out_pcma, file = filename("pcma", seed, dataset)))

  }

  # HILMA
  if (!file.exists(filename("hilma", seed, dataset))){
    set.seed(123)
    Y1 <- scale(Y)
    M1 <- scale(M)
    A1 <- scale(A)

    out_hilma <-  hilma(as.vector(Y1),M1,matrix(A1,ncol = 1), center = T)
    out_hilma$sdy <- sd(Y)
    out_hilma$sda <- sd(A)
    save(out_hilma, file = filename("hilma", seed, dataset))

  }
  
  # PMED
  if (!file.exists(filename("pmed", seed, dataset))){

    set.seed(123)
    out_pmed <- NULL
    coef <- glmnet::coef.glmnet
    try({out_pmed <- pmed_wrapper(A = A, M = M, Y = Y,
                                  C = matrix(1,nrow(M),1),
                                  nlambda = 50)})
    if(is.null(out_pmed)){
      pmed_tie <- 0
    }else{
      pmed_tie <- out_pmed$summary_result[2,2]
    }

    save(pmed_tie, file = filename("pmed", seed, dataset))

  }
  
  # One at a time
  if (!file.exists(filename("uni", seed, dataset))){
    set.seed(123)
    out_uni <- single_mediation(A,M,Y)
    saveRDS(out_uni, file = filename("uni", seed, dataset))
  }

  # BSLMM
  if (!file.exists(filename("bslmm", seed, dataset))){
    set.seed(123)
    out_bslmm <- bslmm_wrapper(A, M, Y)
    saveRDS(out_bslmm, file = filename("bslmm", seed, dataset))

  }

  # Pathway LASSO
  if (!file.exists(filename("plasso", seed, dataset))){
    set.seed(123)
    out_plasso <- plasso_wrapper(A, M, Y)
    save(out_plasso, file = filename("plasso", seed, dataset))

  }


  
}
