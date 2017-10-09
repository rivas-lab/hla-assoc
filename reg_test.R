
dosage <- as.data.frame(readRDS("ukb_hla_v2.rds"))
rounded_dosage <- as.data.frame(readRDS("ukb_hla_v2_rounded.rds"))
covar <- readRDS("covar_16698_remove.rds")
#phe <- readRDS("HC303_remove.rds")
phe <- as.data.frame(read.table("HC303.phe", header=FALSE, row.names=1))
outname <- "dosage_out"

print(nrow(dosage))
print(nrow(rounded_dosage))
print(nrow(covar))
phe <- phe[rownames(covar),]
print(nrow(phe))
print("done loading")

print(phe[1:10,])

covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

if (dim(covars)[1] != dim(dosage)[1]) {
	dosage <- dosage[rownames(covars),]
        rounded_dosage <- dosage[rownames(covars),]
}

results.gen <- vector("list", dim(dosage)[2])
names(results.gen) <- colnames(dosage)
results.add <- vector("list", dim(dosage)[2])
names(results.add) <- colnames(dosage)

for (gene in colnames(dosage)) {
	covars["gcounts"] <- as.numeric(dosage[,gene])
	fit <- summary(glm(status ~ as.factor(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.gen[[gene]] <- fit
	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.add[[gene]] <- fit
}

saveRDS(results.gen, paste(outname, "_gen.rds", sep=""))
saveRDS(results.add, paste(outname, "_add.rds", sep=""))
