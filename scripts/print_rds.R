args <- commandArgs(TRUE)
argsLen <- length(args)
fn <- args[1]

options(warn=1)

df <- readRDS(fn)
print(df)

#if (argsLen > 1) {
create_out <- args[2]

#write.csv(df,create_out)
#print(create_out)
sink(create_out)
df
sink()

#}


