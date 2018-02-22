args <- commandArgs(TRUE)
root <- args[1]
options(warn=1)
model <- "gen"
#model <- "add"

dosage_add <- readRDS(paste0(root,"_rounded_",model,".rds"))
#dosage_add <- readRDS(root)


sink(paste0(root, "_rounded_",model,".txt"))
#sink("results_BMA/bma_HC79.txt")

dosage_add
sink()

