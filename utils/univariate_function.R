require(dplyr)
require(purrr)
require(tidyr)
require(tibble)

single_mediation <- function(A,M,Y){
  
  alphas <- map_dfr(as.data.frame(M), 
                    ~ summary(
                      lm(y ~ ., data = data.frame(y = .x, A = A))
                    )$coefficients[2,c(1,4)])
  alphas <-
    alphas %>% 
    rename(alpha_hat = 1, alpha_pv = 2) %>% 
    mutate(mediator = paste0("M",1:2000),
           alpha_pv_bonf = pmin(alpha_pv * n(),1),
           alpha_pv_fdr = p.adjust(alpha_pv, method = "fdr")) %>% 
    dplyr::select(mediator,everything())
  
  
  betas <- map_dfr(as.data.frame(M), 
                   ~ summary(
                     lm(y ~ ., data = data.frame(y = Y, x = .x, A = A))
                   )$coefficients[2,c(1,4)])
  betas <-
    betas %>% 
    rename(beta_hat = 1, beta_pv = 2) %>% 
    mutate(mediator = paste0("M",1:2000),
           beta_pv_bonf = pmin(beta_pv * n(),1),
           beta_pv_fdr = p.adjust(beta_pv, method = "fdr")) %>% 
    dplyr::select(mediator,everything())
  
  out_uni <-
    full_join(betas,alphas, by = "mediator") %>% 
    mutate(ab_hat = alpha_hat * beta_hat,
           ab_pv = pmax(alpha_pv,beta_pv),
           ab_pv_fdr = pmax(alpha_pv_fdr,beta_pv_fdr),
           ab_pv_fdr1 = p.adjust(ab_pv,method = "fdr"),
           ab_pv_bonf = pmin(ab_pv * n(), 1))
  
  return(out_uni)
  
}

