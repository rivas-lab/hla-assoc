print("hello world")
#dosage <- readRDS("ukb_hla_v2_remove.rds")
rounded_dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_rounded_remove.rds")
dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_remove.rds")

print(dosage[1:10,])
#print(rounded_dosage[1:10,])

rounded_dosage[is.na(rounded_dosage)] <- -1
print(rounded_dosage[1:10,])

print(paste0("dosage nonzeros: ", sum(dosage != 0)))
print(paste0("rounded dosage nonzeros: ", sum(rounded_dosage != 0)))
print(paste0("dosage missing: ", sum(rounded_dosage == -1)))


#print(rounded_dosage)
