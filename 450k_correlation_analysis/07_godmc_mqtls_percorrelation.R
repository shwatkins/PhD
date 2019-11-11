library(data.table)

# list files with the correlating pairs (files represent 0.1 bands of correlation, 
# some correlation bands may be spread over several files):
files = list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

godmc.mqtl <- fread("/path/to/assoc_meta_all.csv")
dim(godmc.mqtl)
head(godmc.mqtl)

load("/path/to/probeDetails_meffil.Rdata")

mqtl_investigation <- function(i){
  load(i)
  print(i)
  print(head(dat))
  dat <- setDT(dat)
  print(dim(dat))
  dat <- dat[complete.cases(dat)]
  print(head(dat))
  print(dim(dat))
  dat$Var1.chr <- probeDetails$chromosome[match(dat$Var1, probeDetails$name)]
  dat$Var2.chr <- probeDetails$chromosome[match(dat$Var2, probeDetails$name)]
  dat$Var1.position <- probeDetails$position[match(dat$Var1, probeDetails$name)]
  dat$Var2.position <- probeDetails$position[match(dat$Var2, probeDetails$name)]
  dat$absdistance <- abs(dat$Var1.position - dat$Var2.position)
  dat$cis_or_trans <- as.character("trans")
  dat$cis_or_trans[dat$Var1.chr==dat$Var2.chr & dat$absdistance < 1000000] <- "cis"
  dat$Var1.mqtl <- as.numeric(c(0))
  dat[, Var1.mqtl := Var1.mqtl][Var1 %in% godmc.mqtl$cpg, Var1.mqtl := 1]
  dat$Var1.mqtl <- as.numeric(dat$Var1.mqtl)
  print(dat[1:10,])
  dat$Var2.mqtl <- as.numeric(c(0))
  dat[, Var2.mqtl := Var2.mqtl][Var2 %in% godmc.mqtl$cpg, Var2.mqtl := 1]
  dat$Var2.mqtl <- as.numeric(dat$Var2.mqtl)
  print(dat[1:10,])
  dat$correlation_total <- dat$Var1.mqtl + dat$Var2.mqtl
  print(dat[1:10,])
  F7_mqtls <- as.data.frame(table(dat$correlation_total,dat$cis_or_trans))
  print(F7_mqtls)
  return(F7_mqtls)
}

F7_mqtls_correlations <- lapply(files, FUN = mqtl_investigation)

names(F7_mqtls_correlations) <- substring(files,41)
print(F7_mqtls_correlations)

save(F7_mqtls_correlations, file="/path/to/F7_mQTL_results_percorrelation.Rdata")