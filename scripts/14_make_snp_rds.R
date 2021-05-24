#!/bin/R

# file names
dosage_name <- "../output/hla_snp_dosage.ped"
fam <- "/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.fam"
white_british <- "/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe"

print("Reading in dosage and covariate files...")
dosage <- as.data.frame(read.table(paste(dosage_name), header=TRUE, row.names=1))

# get order of sample ids
ids <- read.table(fam, colClasses=c(rep("integer", 1), rep("NULL", 5)))

# get sample ids to keep
print("Subsetting to white British individuals...")
to_keep <- read.table(white_british, colClasses=c(rep("integer", 1), rep("NULL", 1)))
inds <- which(ids[,1] %in% to_keep[,1])
dosage <- dosage[inds,]

print(head(dosage))
# save
print("Saving results to file...")
saveRDS(dosage, "../hla_genotype_data/hla_snp_wb.rds")
