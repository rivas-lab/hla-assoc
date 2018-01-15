save_RDS <- TRUE
covar_RDS <- TRUE
test <- TRUE

# file names
if (test) {
    dosage_name = "../data/ukbb_files/test_100.txt"
    root <- "test_100"
    fam <- "/output/create_fam/test_100.fam"
    covar <- "../data/ukbb_files/test_100_covar.phe"
} else {
    dosage_name <- "/scratch/PI/mrivas/ukbb/24983/hla/ukb_hla_v2.txt"
    root <- "ukb_hla_v2"
    fam <- "ukb_hla_v2.fam"
    covar <- "/scratch/PI/mrivas/ukbb/24983/phe/covar_16698.phe"
}

remove <- "/scratch/PI/mrivas/ukbb/24983/phe/sampleqc_16698.phe"

dosage <- as.data.frame(read.table(paste(dosage_name), header=TRUE))
rounded_dosage <- as.data.frame(read.table(paste(root, "_rounded.txt", sep=""), header=TRUE))
covar <- as.data.frame(read.table(covar, header=TRUE, row.names=1))

row.names(dosage) <- row.names(covar)
row.names(rounded_dosage) <- row.names(covar)

# get order of sample ids
ids <- read.table(fam, colClasses=c(rep("integer", 1), rep("NULL", 5)))
#print(typeof(ids))
#print(length(row.names(dosage)))
#print(length(ids))
#print(ids)
#rownames(dosage) <- ids
#rownames(rounded_dosage) <- ids
#print("here")



# get sample ids to remove
to_remove <- read.table(remove, colClasses=c(rep("integer", 1), rep("NULL", 1)))

# find indices of rows to be removed in dosage
inds <- which(ids[,1] %in% to_remove[,1])

# remove rows
dosage <- dosage[-inds,]
rounded_dosage <- rounded_dosage[-inds,]
covar <- covar[-inds,]
#dosage <- dosage[-to_remove]
#rounded_dosage <- rounded_dosage[-to_remove]
#covar <- covar[-to_remove]
rounded_dosage[rounded_dosage == -1] <- NA
#print(rounded_dosage)
#print(dosage)

if (covar_RDS) {
    saveRDS(covar, paste("output/make_dosage_rds/",root,"_covar_16698_remove.rds", sep=""))
}

if (save_RDS) {
    saveRDS(dosage, paste("output/make_dosage_rds/",root,"_remove.rds", sep=""))
    saveRDS(rounded_dosage, paste("output/make_dosage_rds/",root, "_rounded_remove.rds", sep=""))
}
