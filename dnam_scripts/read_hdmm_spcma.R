# One final script to view the results for SPCMA and HDMM, as described in the 
# paper. This file does not produce output tables but simply prints out the
# results that we ended up commenting on. 


# SPCMA -------------------------------------------------------------------
l1 <- ls()

load("dnam_files/multivariate_results/out_spcma.rda")
ls()
attach(out_spcma$SPCA)

# Are any sparse PCs significant?
IE <- as.data.frame(IE)
summary(IE)
sum(IE$adjpv < 0.1) # Number of significant PCs (10% FDR)
sig_pcs <- which(IE$adjpv < 0.1)
IE[sig_pcs,]

#How many variables contribute to the significant sparse PCs?
n_cpgs <- apply(W,2,function(x){return(sum(x!=0))})
summary(n_cpgs)
n_cpgs[sig_pcs]

rm(list=ls())


# HDMM --------------------------------------------------------------------
library(mediation)
load("dnam_files/multivariate_results/out_hdmm.rda")
ls()

# Is the latent mediator based on the first direction of mediation significant?
# What are the estimated mediation effects?
summary(out_hdmm$med)

# What CpGs contribute most to the latent mediator and what are their 
# contributions?
pdms <- out_hdmm$pdms$wt
out_hdmm$pdms[order(-abs(pdms))[1:5], ]

