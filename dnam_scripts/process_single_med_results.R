# Script to read CpG screening output

# This script identifies the 2000 mediators most closely-related to the
# exposure based on the mediator-model p-values produced by running the
# script "dnam_scripts/cpg_screening_models.R"

library(bigreadr)

# Function to get file path based on index
path <- function(i) paste0("dnam_files/single_mediator/ma_out", i, ".csv")

index <- 1
mylist <- c()
while(file.exists(path(index))){ 
  mylist[[index]] <- fread2(path(index))
  index <- index + 1
}

# Merge results
results <- as.data.frame(do.call(rbind, mylist))

# Identify top 2000 CpGs
p <- 2000

results1 <- results[order(results$alpha_pv)[1:p],]

want_cpgs <- data.frame(cpg = results1$cpg)

write.csv(want_cpgs, file = "dnam_files/want_cpgs.csv", 
          row.names = F, quote = F)

# Now get outcome model results
path_ym <- function(i) paste0("dnam_files/single_mediator/ym_out", i, ".csv")
index <- index - 1
mylist <- c()
for(i in 1:index){
  mylist[[index]] <- fread2(path_ym(i))
}
results_ym <- do.call(rbind, mylist)


single_mediator_results <- 
  merge(results1, results_ym, by = "cpg", all.x = T, all.y = F, sort = F)

saveRDS(single_mediator_results, 
        file = "dnam_files/single_mediator/single_med_results.rds")





