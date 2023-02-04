# Script to run all high-dimensional mediation functions sequentially

# Either run this script (which will take a long time, potentially >24 hours)
# or run the scripts in "dnam_scripts/implement_methods" one by one. 
# We recommend the latter so that you can see how the methods were implemenented
# and potentially reduce parameters to make them run faster. For example, 
# in pathway LASSO you should reduce the length of the vector "lambdas" (by 
# taking out some of the larger values); in SPCMA you should reduce the number
# of PCs; and in BSLMM you should reduce the number "ndraws" and "burnin". 
# These parameters are set to what we used in the analysis, not what is fastest.


files <- list.files("dnam_scripts/implement_methods")
for(f in files){
  source(f)
}