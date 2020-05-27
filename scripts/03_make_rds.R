#!/bin/R

# file names
dosage_name <- "/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v2.txt"
rounded_dosage_name <- "/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v2_rounded.txt"
fam <- "/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.fam"
covar <- "/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe"

white_british <- "/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe"
redacted <- "/oak/stanford/groups/mrivas/ukbb24983/sqc/w24983_20200204.csv"

print("Reading in dosage and covariate files...")
dosage <- as.data.frame(read.table(paste(dosage_name), header=TRUE))
rounded_dosage <- as.data.frame(read.table(paste(rounded_dosage_name), header=TRUE))
covar <- as.data.frame(read.table(covar, header=TRUE, row.names=1))
print("Subsetting to covariates of interest...")
covar <- covar[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "N_CNV", "LEN_CNV")]

row.names(dosage) <- row.names(covar)
row.names(rounded_dosage) <- row.names(covar)

# get order of sample ids
ids <- read.table(fam, colClasses=c(rep("integer", 1), rep("NULL", 5)))

# get sample ids to keep
print("Subsetting to white British individuals...")
to_keep <- read.table(white_british, colClasses=c(rep("integer", 1), rep("NULL", 1)))
inds <- which(ids[,1] %in% to_keep[,1])
dosage <- dosage[inds,]
rounded_dosage <- rounded_dosage[inds,]
covar <- covar[inds,]

# remove redacted individuals
print("Removing redacted individuals...")
redacted_phe <- read.table(redacted, colClasses=c("integer"))
to_remove <- which(ids[,1] %in% redacted_phe[,1])
dosage <- dosage[-to_remove]
rounded_dosage <- rounded_dosage[-to_remove]
covar <- covar[-to_remove]
rounded_dosage[rounded_dosage == -1] <- NA

# save
print("Saving results to file...")
saveRDS(covar, "../hla_genotype_data/ukb_hla_v3_covar_24983_wb.rds")
saveRDS(dosage, "../hla_genotype_data/ukb_hla_v3_24983_wb.rds")
saveRDS(rounded_dosage, "../hla_genotype_data/ukb_hla_v3_24983_rounded_wb.rds")
