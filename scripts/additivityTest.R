args <- commandArgs(TRUE)
name <- args[1]
chrom <- args[2]
phe <- args[3]

dy <- "/oak/stanford/groups/mrivas/users/cdeboeve/repos/rivas-lab/ukb-ptv-phewas/private_output/process_data/"
counts <- paste(dy, "chr", chrom, "_ptv_counts.rds", sep="")

# counts <- read.table(gzfile(counts), header=T, row.names=1, check.names=F, nrows=337209)
counts <- readRDS(counts)
counts[counts > 2] <- 2
phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))
phe <- phe[rownames(counts),]
covars <- as.data.frame(read.table("/scratch/PI/mrivas/ukbb/24983/phe/covar_16698.phe", 
				   header=T, row.names=1))
covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]
covars <- covars[rownames(counts),]
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

if (dim(covars)[1] != dim(counts)[1]) {
	counts <- counts[rownames(covars),]
}

results.gen <- vector("list", dim(counts)[2])
names(results.gen) <- colnames(counts)
results.add <- vector("list", dim(counts)[2])
names(results.add) <- colnames(counts)

for (gene in colnames(counts)) {
	covars["gcounts"] <- as.numeric(counts[,gene])
	fit <- summary(glm(status ~ as.factor(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.gen[[gene]] <- fit
	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.add[[gene]] <- fit
}

saveRDS(results.gen, paste(name, "_", chrom, "_gen.rds", sep=""))
saveRDS(results.add, paste(name, "_", chrom, "_add.rds", sep=""))
