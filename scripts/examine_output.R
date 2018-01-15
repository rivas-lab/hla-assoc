args <- commandArgs(TRUE)
root <- args[1]
options(warn=1)

dosage_add <- readRDS(paste0(root,"_rounded_add.rds"))
#dosage_add <- readRDS(root)


sink(paste0(root, "_rounded_add.txt"))
#sink("results_BMA/bma_HC79.txt")

dosage_add
sink()

