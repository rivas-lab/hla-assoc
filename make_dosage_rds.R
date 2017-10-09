# file names
#dosage_name <- "/scratch/PI/mrivas/ukbb/24983/hla/ukb_hla_v2.txt"
dosage_name = "test_100.txt"
root <- "test_100"
#root <- "ukb_hla_v2"
fam <- "ukb_hla_v2.fam"
#fam <- "test_100.fam"
remove <- "/scratch/PI/mrivas/ukbb/24983/phe/sampleqc_16698.phe"
save_RDS <- FALSE
covar_RDS <- TRUE
covar <- "/scratch/PI/mrivas/ukbb/24983/phe/covar_16698.phe"
phe <- "HC303.phe"

dosage <- read.table(paste(dosage_name), header=TRUE)
rounded_dosage <- read.table(paste(root, "_rounded.txt", sep=""), header=TRUE)
covar <- read.table(covar, header=TRUE, row.names=1)
phe <- read.table(phe, header=FALSE)

print(nrow(covar))

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
phe <- phe[-inds,]
#dosage <- dosage[-to_remove]
#rounded_dosage <- rounded_dosage[-to_remove]
#covar <- covar[-to_remove]
rounded_dosage[rounded_dosage == -1] <- NA
#print(rounded_dosage)
#print(dosage)
print(nrow(dosage))
print(nrow(covar))

if (covar_RDS) {
    saveRDS(phe, "HC303_remove.rds")
    saveRDS(covar, "covar_16698_remove.rds")
}
#write.table(dosage, file=paste(root,"_remove.txt",sep=""))
#write.table(rounded_dosage, file=paste(root,"_rounded_remove.txt",sep=""))
if (save_RDS) {
#    dosage <- read.table(paste(root, "txt", sep="."), header=TRUE)
    saveRDS(dosage, paste(root,"_remove.rds", sep=""))
    
#    rounded_dosage <- read.table(paste(root, "_rounded.txt", sep=""), header=TRUE)
    saveRDS(rounded_dosage, paste(root, "_rounded_remove.rds", sep=""))
}
