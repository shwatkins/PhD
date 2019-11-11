library(data.table)
############ files1

# path to correlation band files. 
files <- list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

# split the file list into 3, because for the 450k the dataframes will be
# too big to work in R.
files1 <- files[1:25]
files2 <- files[26:48]
files3 <- files[49:length(files)]

# create list of 450k probes inc gene and position
load("/path/to/probeDetails_meffil.Rdata")
class(probeDetails$chromosome)
# character
table(probeDetails$chromosome)
# Reduce to the columns we need:
genelist <- probeDetails[,c(1,2,3)]
head(genelist)
dim(genelist)
genelist <- na.omit(genelist)
dim(genelist)
genelist[1:5,]
#genelist <- genelist[,-1]
genelist[1:5,]

resultslist_1 <- list()

# function identifies and extracts all cis correlations
cismqtl_investigation <- function(i){
  load(i)
  print(i)
  print(head(dat))
  dat <- setDT(dat)
  print(dim(dat))
  
  Var1 <- as.character(dat$Var1)
  Var1 <- unique(Var1)
  Var2 <- as.character(dat$Var2)
  Var2 <- unique(Var2)
  genelistVar1 <- genelist[Var1,]
  genelistVar1.sort <- genelistVar1[match(Var1, rownames(genelistVar1)), ]
  genelistVar2 <- genelist[Var2,]
  genelistVar2.sort <- genelistVar2[match(Var2, rownames(genelistVar2)), ]
  
  dat$Var1.gene <- genelistVar1$gene.symbol[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.gene <- genelistVar2$gene.symbol[match(dat$Var2, rownames(genelistVar2))]
  dat$Var1.position <- genelistVar1$position[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.position <- genelistVar2$position[match(dat$Var2, rownames(genelistVar2))]
  dat$distance <- (dat$Var1.position - dat$Var2.position)
  dat$Var1.chromosome <- genelistVar1$chromosome[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.chromosome <- genelistVar2$chromosome[match(dat$Var2, rownames(genelistVar2))]
  dat$absdistance <- abs(dat$distance)
  dat$cis_or_trans <- as.character(c("trans"))
  dat[, cis_or_trans := cis_or_trans][absdistance <= 1000000 & Var1.chromosome == Var2.chromosome, cis_or_trans := "cis"]
  print(head(dat))
  print(summary(dat$cis_or_trans))
  print(table(dat$cis_or_trans))
  # new table for results
  resultslist_1[[i]] <- table(dat$cis_or_trans)
  cis_dat <- dat[dat$cis_or_trans == "cis",]
  print(dim(cis_dat))
  cis_dat <- cis_dat[,c("Var1","Var2","value", "absdistance", "Var1.chromosome")]
  return(cis_dat)
}

F7_cis_correlations_1 <- lapply(files1, FUN = cismqtl_investigation)

names(F7_cis_correlations_1) <- substring(files1,41)
print(F7_cis_correlations_1)

save(F7_cis_correlations_1, file="/path/to/F7_cis_fordecayplots_1.Rdata")


############ files2

files <- list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

files1 <- files[1:25]
files2 <- files[26:48]
files3 <- files[49:length(files)]

# create list of 450k probes inc gene and position
load("/path/to/probeDetails_meffil.Rdata")
class(probeDetails$chromosome)
# character
table(probeDetails$chromosome)
# Reduce to the columns we need:
genelist <- probeDetails[,c(1,2,3)]
head(genelist)
dim(genelist)
genelist <- na.omit(genelist)
dim(genelist)
genelist[1:5,]
#genelist <- genelist[,-1]
genelist[1:5,]

resultslist_2 <- list()

cismqtl_investigation <- function(i){
  load(i)
  print(i)
  print(head(dat))
  dat <- setDT(dat)
  print(dim(dat))
  
  Var1 <- as.character(dat$Var1)
  Var1 <- unique(Var1)
  Var2 <- as.character(dat$Var2)
  Var2 <- unique(Var2)
  genelistVar1 <- genelist[Var1,]
  genelistVar1.sort <- genelistVar1[match(Var1, rownames(genelistVar1)), ]
  genelistVar2 <- genelist[Var2,]
  genelistVar2.sort <- genelistVar2[match(Var2, rownames(genelistVar2)), ]
  
  dat$Var1.gene <- genelistVar1$gene.symbol[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.gene <- genelistVar2$gene.symbol[match(dat$Var2, rownames(genelistVar2))]
  dat$Var1.position <- genelistVar1$position[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.position <- genelistVar2$position[match(dat$Var2, rownames(genelistVar2))]
  dat$distance <- (dat$Var1.position - dat$Var2.position)
  dat$Var1.chromosome <- genelistVar1$chromosome[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.chromosome <- genelistVar2$chromosome[match(dat$Var2, rownames(genelistVar2))]
  dat$absdistance <- abs(dat$distance)
  dat$cis_or_trans <- as.character(c("trans"))
  dat[, cis_or_trans := cis_or_trans][absdistance <= 1000000 & Var1.chromosome == Var2.chromosome, cis_or_trans := "cis"]
  print(head(dat))
  print(summary(dat$cis_or_trans))
  print(table(dat$cis_or_trans))
  # new table for results
  resultslist_2[[i]] <- table(dat$cis_or_trans)
  cis_dat <- dat[dat$cis_or_trans == "cis",]
  print(dim(cis_dat))
  cis_dat <- cis_dat[,c("Var1","Var2","value", "absdistance", "Var1.chromosome")]
  return(cis_dat)
}

F7_cis_correlations_2 <- lapply(files2, FUN = cismqtl_investigation)

names(F7_cis_correlations_2) <- substring(files2,41)
print(F7_cis_correlations_2)

save(F7_cis_correlations_2, file="/path/to/F7_cis_fordecayplots_2.Rdata")


############ files3

files <- list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

files1 <- files[1:25]
files2 <- files[26:48]
files3 <- files[49:length(files)]

# create list of 450k probes inc gene and position
load("/path/to/probeDetails_meffil.Rdata")
class(probeDetails$chromosome)
# character
table(probeDetails$chromosome)
# Reduce to the columns we need:
genelist <- probeDetails[,c(1,2,3)]
head(genelist)
dim(genelist)
genelist <- na.omit(genelist)
dim(genelist)
genelist[1:5,]
#genelist <- genelist[,-1]
genelist[1:5,]

resultslist_3 <- list()

cismqtl_investigation <- function(i){
  load(i)
  print(i)
  print(head(dat))
  dat <- setDT(dat)
  print(dim(dat))
  
  Var1 <- as.character(dat$Var1)
  Var1 <- unique(Var1)
  Var2 <- as.character(dat$Var2)
  Var2 <- unique(Var2)
  genelistVar1 <- genelist[Var1,]
  genelistVar1.sort <- genelistVar1[match(Var1, rownames(genelistVar1)), ]
  genelistVar2 <- genelist[Var2,]
  genelistVar2.sort <- genelistVar2[match(Var2, rownames(genelistVar2)), ]
  
  dat$Var1.gene <- genelistVar1$gene.symbol[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.gene <- genelistVar2$gene.symbol[match(dat$Var2, rownames(genelistVar2))]
  dat$Var1.position <- genelistVar1$position[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.position <- genelistVar2$position[match(dat$Var2, rownames(genelistVar2))]
  dat$distance <- (dat$Var1.position - dat$Var2.position)
  dat$Var1.chromosome <- genelistVar1$chromosome[match(dat$Var1, rownames(genelistVar1))]
  dat$Var2.chromosome <- genelistVar2$chromosome[match(dat$Var2, rownames(genelistVar2))]
  dat$absdistance <- abs(dat$distance)
  dat$cis_or_trans <- as.character(c("trans"))
  dat[, cis_or_trans := cis_or_trans][absdistance <= 1000000 & Var1.chromosome == Var2.chromosome, cis_or_trans := "cis"]
  print(head(dat))
  print(summary(dat$cis_or_trans))
  print(table(dat$cis_or_trans))
  # new table for results
  resultslist_3[[i]] <- table(dat$cis_or_trans)
  cis_dat <- dat[dat$cis_or_trans == "cis",]
  print(dim(cis_dat))
  cis_dat <- cis_dat[,c("Var1","Var2","value", "absdistance", "Var1.chromosome")]
  return(cis_dat)
}

F7_cis_correlations_3 <- lapply(files3, FUN = cismqtl_investigation)

names(F7_cis_correlations_3) <- substring(files3,41)
print(F7_cis_correlations_3)

save(F7_cis_correlations_3, file="/path/to/F7_cis_fordecayplots_3.Rdata")


###################################################

# now for each of the 3 files of cis correlations, go through and create a file
# containing correlations on each chromosome

library(data.table)
load("/path/to/F7_cis_fordecayplots_1.Rdata")

chr1 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr1",])
save(chr1, file="/path/to/cis_correlations/chr1_1.Rdata")
print(head(chr1[1]))
rm(chr1)

chr2 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr2",])
save(chr2, file="/path/to/cis_correlations/chr2_1.Rdata")
print(head(chr2[1]))
rm(chr2)

chr3 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr3",])
save(chr3, file="/path/to/cis_correlations/chr3_1.Rdata")
print(head(chr3[1]))
rm(chr3)

chr4 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr4",])
save(chr4, file="/path/to/cis_correlations/chr4_1.Rdata")
print(head(chr4[1]))
rm(chr4)

chr5 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr5",])
save(chr5, file="/path/to/cis_correlations/chr5_1.Rdata")
print(head(chr5[1]))
rm(chr5)

chr6 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr6",])
save(chr6, file="/path/to/cis_correlations/chr6_1.Rdata")
print(head(chr6[1]))
rm(chr6)

chr7 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr7",])
save(chr7, file="/path/to/cis_correlations/chr7_1.Rdata")
print(head(chr7[1]))
rm(chr7)

chr8 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr8",])
save(chr8, file="/path/to/cis_correlations/chr8_1.Rdata")
print(head(chr8[1]))
rm(chr8)

chr9 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr9",])
save(chr9, file="/path/to/cis_correlations/chr9_1.Rdata")
print(head(chr9[1]))
rm(chr9)

chr10 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr10",])
save(chr10, file="/path/to/cis_correlations/chr10_1.Rdata")
print(head(chr10[1]))
rm(chr10)

chr11 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr11",])
save(chr11, file="/path/to/cis_correlations/chr11_1.Rdata")
print(head(chr11[1]))
rm(chr11)

chr12 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr12",])
save(chr12, file="/path/to/cis_correlations/chr12_1.Rdata")
print(head(chr12[1]))
rm(chr12)

chr13 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr13",])
save(chr13, file="/path/to/cis_correlations/chr13_1.Rdata")
print(head(chr13[1]))
rm(chr13)

chr14 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr14",])
save(chr14, file="/path/to/cis_correlations/chr14_1.Rdata")
print(head(chr14[1]))
rm(chr14)

chr15 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr15",])
save(chr15, file="/path/to/cis_correlations/chr15_1.Rdata")
print(head(chr15[1]))
rm(chr15)

chr16 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr16",])
save(chr16, file="/path/to/cis_correlations/chr16_1.Rdata")
print(head(chr16[1]))
rm(chr16)

chr17 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr17",])
save(chr17, file="/path/to/cis_correlations/chr17_1.Rdata")
print(head(chr17[1]))
rm(chr17)

chr18 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr18",])
save(chr18, file="/path/to/cis_correlations/chr18_1.Rdata")
print(head(chr18[1]))
rm(chr18)

chr19 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr19",])
save(chr19, file="/path/to/cis_correlations/chr19_1.Rdata")
print(head(chr19[1]))
rm(chr19)

chr20 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr20",])
save(chr20, file="/path/to/cis_correlations/chr20_1.Rdata")
print(head(chr20[1]))
rm(chr20)

chr21 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr21",])
save(chr21, file="/path/to/cis_correlations/chr21_1.Rdata")
print(head(chr21[1]))
rm(chr21)

chr22 <- lapply(F7_cis_correlations_1, function(x) x[x$Var1.chromosome=="chr22",])
save(chr22, file="/path/to/cis_correlations/chr22_1.Rdata")
print(head(chr22[1]))
rm(chr22)

print("done")

library(data.table)
load("/path/to/F7_cis_fordecayplots_2.Rdata")

chr1_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr1",])
save(chr1_2, file="/path/to/cis_correlations/chr1_2.Rdata")
print(head(chr1_2[1]))
rm(chr1_2)

chr2_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr2",])
save(chr2_2, file="/path/to/cis_correlations/chr2_2.Rdata")
print(head(chr2_2[1]))
rm(chr2_2)

chr3_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr3",])
save(chr3_2, file="/path/to/cis_correlations/chr3_2.Rdata")
print(head(chr3_2[1]))
rm(chr3_2)

chr4_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr4",])
save(chr4_2, file="/path/to/cis_correlations/chr4_2.Rdata")
print(head(chr4_2[1]))
rm(chr4_2)

chr5_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr5",])
save(chr5_2, file="/path/to/cis_correlations/chr5_2.Rdata")
print(head(chr5_2[1]))
rm(chr5_2)

chr6_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr6",])
save(chr6_2, file="/path/to/cis_correlations/chr6_2.Rdata")
print(head(chr6_2[1]))
rm(chr6_2)

chr7_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr7",])
save(chr7_2, file="/path/to/cis_correlations/chr7_2.Rdata")
print(head(chr7_2[1]))
rm(chr7_2)

chr8_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr8",])
save(chr8_2, file="/path/to/cis_correlations/chr8_2.Rdata")
print(head(chr8_2[1]))
rm(chr8_2)

chr9_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr9",])
save(chr9_2, file="/path/to/cis_correlations/chr9_2.Rdata")
print(head(chr9_2[1]))
rm(chr9_2)

chr10_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr10",])
save(chr10_2, file="/path/to/cis_correlations/chr10_2.Rdata")
print(head(chr10_2[1]))
rm(chr10_2)

chr11_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr11",])
save(chr11_2, file="/path/to/cis_correlations/chr11_2.Rdata")
print(head(chr11_2[1]))
rm(chr11_2)

chr12_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr12",])
save(chr12_2, file="/path/to/cis_correlations/chr12_2.Rdata")
print(head(chr12_2[1]))
rm(chr12_2)

chr13_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr13",])
save(chr13_2, file="/path/to/cis_correlations/chr13_2.Rdata")
print(head(chr13_2[1]))
rm(chr13_2)

chr14_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr14",])
save(chr14_2, file="/path/to/cis_correlations/chr14_2.Rdata")
print(head(chr14_2[1]))
rm(chr14_2)

chr15_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr15",])
save(chr15_2, file="/path/to/cis_correlations/chr15_2.Rdata")
print(head(chr15_2[1]))
rm(chr15_2)

chr16_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr16",])
save(chr16_2, file="/path/to/cis_correlations/chr16_2.Rdata")
print(head(chr16_2[1]))
rm(chr16_2)

chr17_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr17",])
save(chr17_2, file="/path/to/cis_correlations/chr17_2.Rdata")
print(head(chr17_2[1]))
rm(chr17_2)

chr18_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr18",])
save(chr18_2, file="/path/to/cis_correlations/chr18_2.Rdata")
print(head(chr18_2[1]))
rm(chr18_2)

chr19_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr19",])
save(chr19_2, file="/path/to/cis_correlations/chr19_2.Rdata")
print(head(chr19_2[1]))
rm(chr19_2)

chr20_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr20",])
save(chr20_2, file="/path/to/cis_correlations/chr20_2.Rdata")
print(head(chr20_2[1]))
rm(chr20_2)

chr21_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr21",])
save(chr21_2, file="/path/to/cis_correlations/chr21_2.Rdata")
print(head(chr21_2[1]))
rm(chr21_2)

chr22_2 <- lapply(F7_cis_correlations_2, function(x) x[x$Var1.chromosome=="chr22",])
save(chr22_2, file="/path/to/cis_correlations/chr22_2.Rdata")
print(head(chr22_2[1]))
rm(chr22_2)

print("done")

library(data.table)
load("/path/to/F7_cis_fordecayplots_3.Rdata")

chr1_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr1",])
save(chr1_3, file="/path/to/cis_correlations/chr1_3.Rdata")
print(head(chr1_3[1]))
rm(chr1_3)

chr2_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr2",])
save(chr2_3, file="/path/to/cis_correlations/chr2_3.Rdata")
print(head(chr2_3[1]))
rm(chr2_3)

chr3_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr3",])
save(chr3_3, file="/path/to/cis_correlations/chr3_3.Rdata")
print(head(chr3_3[1]))
rm(chr3_3)

chr4_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr4",])
save(chr4_3, file="/path/to/cis_correlations/chr4_3.Rdata")
print(head(chr4_3[1]))
rm(chr4_3)

chr5_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr5",])
save(chr5_3, file="/path/to/cis_correlations/chr5_3.Rdata")
print(head(chr5_3[1]))
rm(chr5_3)

chr6_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr6",])
save(chr6_3, file="/path/to/cis_correlations/chr6_3.Rdata")
print(head(chr6_3[1]))
rm(chr6_3)

chr7_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr7",])
save(chr7_3, file="/path/to/cis_correlations/chr7_3.Rdata")
print(head(chr7_3[1]))
rm(chr7_3)

chr8_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr8",])
save(chr8_3, file="/path/to/cis_correlations/chr8_3.Rdata")
print(head(chr8_3[1]))
rm(chr8_3)

chr9_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr9",])
save(chr9_3, file="/path/to/cis_correlations/chr9_3.Rdata")
print(head(chr9_3[1]))
rm(chr9_3)

chr10_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr10",])
save(chr10_3, file="/path/to/cis_correlations/chr10_3.Rdata")
print(head(chr10_3[1]))
rm(chr10_3)

chr11_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr11",])
save(chr11_3, file="/path/to/cis_correlations/chr11_3.Rdata")
print(head(chr11_3[1]))
rm(chr11_3)

chr12_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr12",])
save(chr12_3, file="/path/to/cis_correlations/chr12_3.Rdata")
print(head(chr12_3[1]))
rm(chr12_3)

chr13_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr13",])
save(chr13_3, file="/path/to/cis_correlations/chr13_3.Rdata")
print(head(chr13_3[1]))
rm(chr13_3)

chr14_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr14",])
save(chr14_3, file="/path/to/cis_correlations/chr14_3.Rdata")
print(head(chr14_3[1]))
rm(chr14_3)

chr15_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr15",])
save(chr15_3, file="/path/to/cis_correlations/chr15_3.Rdata")
print(head(chr15_3[1]))
rm(chr15_3)

chr16_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr16",])
save(chr16_3, file="/path/to/cis_correlations/chr16_3.Rdata")
print(head(chr16_3[1]))
rm(chr16_3)

chr17_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr17",])
save(chr17_3, file="/path/to/cis_correlations/chr17_3.Rdata")
print(head(chr17_3[1]))
rm(chr17_3)

chr18_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr18",])
save(chr18_3, file="/path/to/cis_correlations/chr18_3.Rdata")
print(head(chr18_3[1]))
rm(chr18_3)

chr19_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr19",])
save(chr19_3, file="/path/to/cis_correlations/chr19_3.Rdata")
print(head(chr19_3[1]))
rm(chr19_3)

chr20_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr20",])
save(chr20_3, file="/path/to/cis_correlations/chr20_3.Rdata")
print(head(chr20_3[1]))
rm(chr20_3)

chr21_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr21",])
save(chr21_3, file="/path/to/cis_correlations/chr21_3.Rdata")
print(head(chr21_3[1]))
rm(chr21_3)

chr22_3 <- lapply(F7_cis_correlations_3, function(x) x[x$Var1.chromosome=="chr22",])
save(chr22_3, file="/path/to/cis_correlations/chr22_3.Rdata")
print(head(chr22_3[1]))
rm(chr22_3)

print("done")

#####################################

# now plot the cis decay graphs for each chromosome

library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(Hmisc)

###### chr1 #####

load("/path/to/cis_correlations/chr1_1.Rdata")
load("/path/to/cis_correlations/chr1_2.Rdata")
load("/path/to/cis_correlations/chr1_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr1_1_dim <- lapply(chr1, FUN = getdim)
chr1_1_dim_total <- lapply(chr1_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr1_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr1_2_dim <- lapply(chr1_2, FUN = getdim)
chr1_2_dim_total <- lapply(chr1_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr1_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr1_3_dim <- lapply(chr1_3, FUN = getdim)
chr1_3_dim_total <- lapply(chr1_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr1_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr1_1 <- rbindlist(chr1)
dim(chr1_1)
chr1_2 <- rbindlist(chr1_2)
dim(chr1_2)
chr1_3 <- rbindlist(chr1_3)
dim(chr1_3)

chr1 <- rbind(chr1_1, chr1_2)
dim(chr1)
chr1 <- rbind(chr1, chr1_3)
dim(chr1)
F7_cis_dim <- list()
F7_cis_dim$chr1_F7 <- dim(chr1)

chr1_10k <- chr1[chr1$absdistance <= 10000]
dim(chr1_10k)
head(chr1_10k)
chr1_10k_sort <- chr1_10k[order(absdistance)]
head(chr1_10k_sort)

# split to positive and negative correlations
chr1_10k_pos <- chr1_10k_sort[chr1_10k_sort$value >= 0,]
dim(chr1_10k_pos)
chr1_10k_neg <- chr1_10k_sort[chr1_10k_sort$value < 0,]
dim(chr1_10k_neg)

# bin the positive correlations on distance
chr1_10k_pos$group <- as.numeric(cut2(chr1_10k_pos$absdistance, m=100))
summary(chr1_10k_pos$group)
chr1_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr1_10k_pos_valuemean <- aggregate(chr1_10k_pos[, 3], list(chr1_10k_pos$group), mean)
class(chr1_10k_pos_valuemean)
head(chr1_10k_pos_valuemean)
# calculate the median distance of each bin
chr1_10k_pos_mediandistance <- aggregate(chr1_10k_pos[, 4], list(chr1_10k_pos$group), median)
class(chr1_10k_pos_mediandistance)
head(chr1_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr1_10k_pos_df <- merge(chr1_10k_pos_valuemean, chr1_10k_pos_mediandistance, by=c("Group.1"))
head(chr1_10k_pos_df)

# bin the negative correlations on distance
chr1_10k_neg$group <- as.numeric(cut2(chr1_10k_neg$absdistance, m=100))
summary(chr1_10k_neg$group)
chr1_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr1_10k_neg_valuemean <- aggregate(chr1_10k_neg[, 3], list(chr1_10k_neg$group), mean)
class(chr1_10k_neg_valuemean)
head(chr1_10k_neg_valuemean)
# calculate the median distance of each bin
chr1_10k_neg_mediandistance <- aggregate(chr1_10k_neg[, 4], list(chr1_10k_neg$group), median)
class(chr1_10k_neg_mediandistance)
head(chr1_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr1_10k_neg_df <- merge(chr1_10k_neg_valuemean, chr1_10k_neg_mediandistance, by=c("Group.1"))
head(chr1_10k_neg_df)

# calculate the SD in each positive bin
chr1_10k_pos_sd <- aggregate(chr1_10k_pos[, 3], list(chr1_10k_pos$group), sd)
class(chr1_10k_pos_sd)
head(chr1_10k_pos_sd)
chr1_10k_pos_df <- merge(chr1_10k_pos_df, chr1_10k_pos_sd, by=c("Group.1"))
head(chr1_10k_pos_df)
chr1_10k_pos_df$n_in_group <- as.numeric(table(chr1_10k_pos$group))
head(chr1_10k_pos_df)
table(chr1_10k_pos$group)

# calculate the SD in each negative bin
chr1_10k_neg_sd <- aggregate(chr1_10k_neg[, 3], list(chr1_10k_neg$group), sd)
class(chr1_10k_neg_sd)
head(chr1_10k_neg_sd)
chr1_10k_neg_df <- merge(chr1_10k_neg_df, chr1_10k_neg_sd, by=c("Group.1"))
head(chr1_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr1_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr1_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr1_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr1_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr1_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 1: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr1_1k <- chr1_10k_sort[chr1_10k_sort$absdistance <= 1000,]
chr1_1k$group <- "group"
chr1_1k$group[chr1_1k$value >=-1 & chr1_1k$value < -0.9] <- "-1 to -0.9"
chr1_1k$group[chr1_1k$value >=-0.9 & chr1_1k$value < -0.8] <- "-0.9 to -0.8"
chr1_1k$group[chr1_1k$value >=-0.8 & chr1_1k$value < -0.7] <- "-0.8 to -0.7"
chr1_1k$group[chr1_1k$value >=-0.7 & chr1_1k$value < -0.6] <- "-0.7 to -0.6"
chr1_1k$group[chr1_1k$value >=-0.6 & chr1_1k$value < -0.5] <- "-0.6 to -0.5"
chr1_1k$group[chr1_1k$value >=-0.5 & chr1_1k$value < -0.4] <- "-0.5 to -0.4"
chr1_1k$group[chr1_1k$value >=-0.4 & chr1_1k$value < -0.3] <- "-0.4 to -0.3"
chr1_1k$group[chr1_1k$value >=-0.3 & chr1_1k$value < -0.2] <- "-0.3 to -0.2"
chr1_1k$group[chr1_1k$value >=-0.2 & chr1_1k$value < -0.1] <- "-0.2 to -0.1"
chr1_1k$group[chr1_1k$value >=-0.1 & chr1_1k$value < -0] <- "-0.1 to 0"
chr1_1k$group[chr1_1k$value >=0 & chr1_1k$value < 0.1] <- "0 to 0.1"
chr1_1k$group[chr1_1k$value >=0.1 & chr1_1k$value < 0.2] <- "0.1 to 0.2"
chr1_1k$group[chr1_1k$value >=0.2 & chr1_1k$value < 0.3] <- "0.2 to 0.3"
chr1_1k$group[chr1_1k$value >=0.3 & chr1_1k$value < 0.4] <- "0.3 to 0.4"
chr1_1k$group[chr1_1k$value >=0.4 & chr1_1k$value < 0.5] <- "0.4 to 0.5"
chr1_1k$group[chr1_1k$value >=0.5 & chr1_1k$value < 0.6] <- "0.5 to 0.6"
chr1_1k$group[chr1_1k$value >=0.6 & chr1_1k$value < 0.7] <- "0.6 to 0.7"
chr1_1k$group[chr1_1k$value >=0.7 & chr1_1k$value < 0.8] <- "0.7 to 0.8"
chr1_1k$group[chr1_1k$value >=0.8 & chr1_1k$value < 0.9] <- "0.8 to 0.9"
chr1_1k$group[chr1_1k$value >=0.9 & chr1_1k$value < 1] <- "0.9 to 1"

chr1_1k$group <- factor(chr1_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr1_1k$group)
head(chr1_1k)
cor.group <- as.data.frame(table(chr1_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr1_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 1: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()

# how many probes are there for each correlation band? (ie is there a network)
# calculate the mean correlation in each bin
groups <- as.character(c("-1 to -0.9", "-0.9 to -0.8", 
                         "-0.8 to -0.7", "-0.7 to -0.6",
                         "-0.6 to -0.5", "-0.5 to -0.4",
                         "-0.4 to -0.3", "-0.3 to -0.2",
                         "-0.2 to -0.1", "-0.1 to 0",
                         "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                         "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                         "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                         "0.9 to 1"))

get_probes <- function(x){
  print(x)
  var1 <- as.character(chr1_1k$Var1[chr1_1k$group == i])
  print(length(var1))
  var2 <- as.character(chr1_1k$Var2[chr1_1k$group == i])
  print(length(var2))
  var1var2 <- as.character(c(var1,var2))
  print(length(var1var2))
  var1var2 <- unique(var1var2)
  print(length(var1var2))
  return(var1var2)
}

n_of_probes <- lapply(groups, FUN=get_probes)
summary(n_of_probes)
names(n_of_probes) <- groups
summary(n_of_probes)
n_of_probes[20]

n_of_probes <- list()
for (i in groups){
  print(i)
  var1 <- as.character(chr1_1k$Var1[chr1_1k$group == i])
  print(length(var1))
  var2 <- as.character(chr1_1k$Var2[chr1_1k$group == i])
  print(length(var2))
  var1var2 <- as.character(c(var1,var2))
  print(length(var1var2))
  var1var2 <- unique(var1var2)
  print(length(var1var2))
  if (length(var1var2) == 0) {
    n_of_probes[i] <- 0
  } else {
    n_of_probes[i] <- var1var2}
  names(n_of_probes[i]) <- i
}
n_of_probes[1]
n_of_probes[15]
summary(n_of_probes)
summary(n_of_probes[20])
var1 <- as.character(chr1_1k$Var1[chr1_1k$group == "0.9 to 1"])
print(length(var1))
var2 <- as.character(chr1_1k$Var2[chr1_1k$group == "0.9 to 1"])
print(length(var2))
var1var2 <- as.character(c(var1,var2))
print(length(var1var2))
var1var2 <- unique(var1var2)
print(length(var1var2))

chr1_1k_var1 <- aggregate(chr1_1k[, 1], list(chr1_1k$group), mean)
class(chr1_10k_neg_valuemean)
head(chr1_10k_neg_valuemean)


###### chr2 #####

load("/path/to/cis_correlations/chr2_1.Rdata")
load("/path/to/cis_correlations/chr2_2.Rdata")
load("/path/to/cis_correlations/chr2_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr2_1_dim <- lapply(chr2, FUN = getdim)
chr2_1_dim_total <- lapply(chr2_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr2_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr2_2_dim <- lapply(chr2_2, FUN = getdim)
chr2_2_dim_total <- lapply(chr2_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr2_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr2_3_dim <- lapply(chr2_3, FUN = getdim)
chr2_3_dim_total <- lapply(chr2_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr2_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr2_1 <- rbindlist(chr2)
dim(chr2_1)
chr2_2 <- rbindlist(chr2_2)
dim(chr2_2)
chr2_3 <- rbindlist(chr2_3)
dim(chr2_3)

chr2 <- rbind(chr2_1, chr2_2)
dim(chr2)
chr2 <- rbind(chr2, chr2_3)
dim(chr2)
F7_cis_dim <- list()
F7_cis_dim$chr2_F7 <- dim(chr2)

chr2_10k <- chr2[chr2$absdistance <= 10000]
dim(chr2_10k)
head(chr2_10k)
chr2_10k_sort <- chr2_10k[order(absdistance)]
head(chr2_10k_sort)

# split to positive and negative correlations
chr2_10k_pos <- chr2_10k_sort[chr2_10k_sort$value >= 0,]
dim(chr2_10k_pos)
chr2_10k_neg <- chr2_10k_sort[chr2_10k_sort$value < 0,]
dim(chr2_10k_neg)

# bin the positive correlations on distance
chr2_10k_pos$group <- as.numeric(cut2(chr2_10k_pos$absdistance, m=100))
summary(chr2_10k_pos$group)
chr2_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr2_10k_pos_valuemean <- aggregate(chr2_10k_pos[, 3], list(chr2_10k_pos$group), mean)
class(chr2_10k_pos_valuemean)
head(chr2_10k_pos_valuemean)
# calculate the median distance of each bin
chr2_10k_pos_mediandistance <- aggregate(chr2_10k_pos[, 4], list(chr2_10k_pos$group), median)
class(chr2_10k_pos_mediandistance)
head(chr2_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr2_10k_pos_df <- merge(chr2_10k_pos_valuemean, chr2_10k_pos_mediandistance, by=c("Group.1"))
head(chr2_10k_pos_df)

# bin the negative correlations on distance
chr2_10k_neg$group <- as.numeric(cut2(chr2_10k_neg$absdistance, m=100))
summary(chr2_10k_neg$group)
chr2_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr2_10k_neg_valuemean <- aggregate(chr2_10k_neg[, 3], list(chr2_10k_neg$group), mean)
class(chr2_10k_neg_valuemean)
head(chr2_10k_neg_valuemean)
# calculate the median distance of each bin
chr2_10k_neg_mediandistance <- aggregate(chr2_10k_neg[, 4], list(chr2_10k_neg$group), median)
class(chr2_10k_neg_mediandistance)
head(chr2_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr2_10k_neg_df <- merge(chr2_10k_neg_valuemean, chr2_10k_neg_mediandistance, by=c("Group.1"))
head(chr2_10k_neg_df)

# calculate the SD in each positive bin
chr2_10k_pos_sd <- aggregate(chr2_10k_pos[, 3], list(chr2_10k_pos$group), sd)
class(chr2_10k_pos_sd)
head(chr2_10k_pos_sd)
chr2_10k_pos_df <- merge(chr2_10k_pos_df, chr2_10k_pos_sd, by=c("Group.1"))
head(chr2_10k_pos_df)
chr2_10k_pos_df$n_in_group <- as.numeric(table(chr2_10k_pos$group))
head(chr2_10k_pos_df)
table(chr2_10k_pos$group)

# calculate the SD in each negative bin
chr2_10k_neg_sd <- aggregate(chr2_10k_neg[, 3], list(chr2_10k_neg$group), sd)
class(chr2_10k_neg_sd)
head(chr2_10k_neg_sd)
chr2_10k_neg_df <- merge(chr2_10k_neg_df, chr2_10k_neg_sd, by=c("Group.1"))
head(chr2_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr2_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr2_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr2_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr2_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr2_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 2: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr2_1k <- chr2_10k_sort[chr2_10k_sort$absdistance <= 1000,]
chr2_1k$group <- "group"
chr2_1k$group[chr2_1k$value >=-1 & chr2_1k$value < -0.9] <- "-1 to -0.9"
chr2_1k$group[chr2_1k$value >=-0.9 & chr2_1k$value < -0.8] <- "-0.9 to -0.8"
chr2_1k$group[chr2_1k$value >=-0.8 & chr2_1k$value < -0.7] <- "-0.8 to -0.7"
chr2_1k$group[chr2_1k$value >=-0.7 & chr2_1k$value < -0.6] <- "-0.7 to -0.6"
chr2_1k$group[chr2_1k$value >=-0.6 & chr2_1k$value < -0.5] <- "-0.6 to -0.5"
chr2_1k$group[chr2_1k$value >=-0.5 & chr2_1k$value < -0.4] <- "-0.5 to -0.4"
chr2_1k$group[chr2_1k$value >=-0.4 & chr2_1k$value < -0.3] <- "-0.4 to -0.3"
chr2_1k$group[chr2_1k$value >=-0.3 & chr2_1k$value < -0.2] <- "-0.3 to -0.2"
chr2_1k$group[chr2_1k$value >=-0.2 & chr2_1k$value < -0.1] <- "-0.2 to -0.1"
chr2_1k$group[chr2_1k$value >=-0.1 & chr2_1k$value < -0] <- "-0.1 to 0"
chr2_1k$group[chr2_1k$value >=0 & chr2_1k$value < 0.1] <- "0 to 0.1"
chr2_1k$group[chr2_1k$value >=0.1 & chr2_1k$value < 0.2] <- "0.1 to 0.2"
chr2_1k$group[chr2_1k$value >=0.2 & chr2_1k$value < 0.3] <- "0.2 to 0.3"
chr2_1k$group[chr2_1k$value >=0.3 & chr2_1k$value < 0.4] <- "0.3 to 0.4"
chr2_1k$group[chr2_1k$value >=0.4 & chr2_1k$value < 0.5] <- "0.4 to 0.5"
chr2_1k$group[chr2_1k$value >=0.5 & chr2_1k$value < 0.6] <- "0.5 to 0.6"
chr2_1k$group[chr2_1k$value >=0.6 & chr2_1k$value < 0.7] <- "0.6 to 0.7"
chr2_1k$group[chr2_1k$value >=0.7 & chr2_1k$value < 0.8] <- "0.7 to 0.8"
chr2_1k$group[chr2_1k$value >=0.8 & chr2_1k$value < 0.9] <- "0.8 to 0.9"
chr2_1k$group[chr2_1k$value >=0.9 & chr2_1k$value < 1] <- "0.9 to 1"

chr2_1k$group <- factor(chr2_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr2_1k$group)
cor.group <- as.data.frame(table(chr2_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr2_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 2: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()



###### chr3 #####

load("/path/to/cis_correlations/chr3_1.Rdata")
load("/path/to/cis_correlations/chr3_2.Rdata")
load("/path/to/cis_correlations/chr3_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr3_1_dim <- lapply(chr3, FUN = getdim)
chr3_1_dim_total <- lapply(chr3_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr3_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr3_2_dim <- lapply(chr3_2, FUN = getdim)
chr3_2_dim_total <- lapply(chr3_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr3_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr3_3_dim <- lapply(chr3_3, FUN = getdim)
chr3_3_dim_total <- lapply(chr3_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr3_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr3_1 <- rbindlist(chr3)
dim(chr3_1)
chr3_2 <- rbindlist(chr3_2)
dim(chr3_2)
chr3_3 <- rbindlist(chr3_3)
dim(chr3_3)

chr3 <- rbind(chr3_1, chr3_2)
dim(chr3)
chr3 <- rbind(chr3, chr3_3)
dim(chr3)
F7_cis_dim <- list()
F7_cis_dim$chr3_F7 <- dim(chr3)

chr3_10k <- chr3[chr3$absdistance <= 10000]
dim(chr3_10k)
head(chr3_10k)
chr3_10k_sort <- chr3_10k[order(absdistance)]
head(chr3_10k_sort)

# split to positive and negative correlations
chr3_10k_pos <- chr3_10k_sort[chr3_10k_sort$value >= 0,]
dim(chr3_10k_pos)
chr3_10k_neg <- chr3_10k_sort[chr3_10k_sort$value < 0,]
dim(chr3_10k_neg)

# bin the positive correlations on distance
chr3_10k_pos$group <- as.numeric(cut2(chr3_10k_pos$absdistance, m=100))
summary(chr3_10k_pos$group)
chr3_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr3_10k_pos_valuemean <- aggregate(chr3_10k_pos[, 3], list(chr3_10k_pos$group), mean)
class(chr3_10k_pos_valuemean)
head(chr3_10k_pos_valuemean)
# calculate the median distance of each bin
chr3_10k_pos_mediandistance <- aggregate(chr3_10k_pos[, 4], list(chr3_10k_pos$group), median)
class(chr3_10k_pos_mediandistance)
head(chr3_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr3_10k_pos_df <- merge(chr3_10k_pos_valuemean, chr3_10k_pos_mediandistance, by=c("Group.1"))
head(chr3_10k_pos_df)

# bin the negative correlations on distance
chr3_10k_neg$group <- as.numeric(cut2(chr3_10k_neg$absdistance, m=100))
summary(chr3_10k_neg$group)
chr3_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr3_10k_neg_valuemean <- aggregate(chr3_10k_neg[, 3], list(chr3_10k_neg$group), mean)
class(chr3_10k_neg_valuemean)
head(chr3_10k_neg_valuemean)
# calculate the median distance of each bin
chr3_10k_neg_mediandistance <- aggregate(chr3_10k_neg[, 4], list(chr3_10k_neg$group), median)
class(chr3_10k_neg_mediandistance)
head(chr3_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr3_10k_neg_df <- merge(chr3_10k_neg_valuemean, chr3_10k_neg_mediandistance, by=c("Group.1"))
head(chr3_10k_neg_df)

# calculate the SD in each positive bin
chr3_10k_pos_sd <- aggregate(chr3_10k_pos[, 3], list(chr3_10k_pos$group), sd)
class(chr3_10k_pos_sd)
head(chr3_10k_pos_sd)
chr3_10k_pos_df <- merge(chr3_10k_pos_df, chr3_10k_pos_sd, by=c("Group.1"))
head(chr3_10k_pos_df)
chr3_10k_pos_df$n_in_group <- as.numeric(table(chr3_10k_pos$group))
head(chr3_10k_pos_df)
table(chr3_10k_pos$group)

# calculate the SD in each negative bin
chr3_10k_neg_sd <- aggregate(chr3_10k_neg[, 3], list(chr3_10k_neg$group), sd)
class(chr3_10k_neg_sd)
head(chr3_10k_neg_sd)
chr3_10k_neg_df <- merge(chr3_10k_neg_df, chr3_10k_neg_sd, by=c("Group.1"))
head(chr3_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr3_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr3_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr3_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr3_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr3_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 3: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr3_1k <- chr3_10k_sort[chr3_10k_sort$absdistance <= 1000,]
chr3_1k$group <- "group"
chr3_1k$group[chr3_1k$value >=-1 & chr3_1k$value < -0.9] <- "-1 to -0.9"
chr3_1k$group[chr3_1k$value >=-0.9 & chr3_1k$value < -0.8] <- "-0.9 to -0.8"
chr3_1k$group[chr3_1k$value >=-0.8 & chr3_1k$value < -0.7] <- "-0.8 to -0.7"
chr3_1k$group[chr3_1k$value >=-0.7 & chr3_1k$value < -0.6] <- "-0.7 to -0.6"
chr3_1k$group[chr3_1k$value >=-0.6 & chr3_1k$value < -0.5] <- "-0.6 to -0.5"
chr3_1k$group[chr3_1k$value >=-0.5 & chr3_1k$value < -0.4] <- "-0.5 to -0.4"
chr3_1k$group[chr3_1k$value >=-0.4 & chr3_1k$value < -0.3] <- "-0.4 to -0.3"
chr3_1k$group[chr3_1k$value >=-0.3 & chr3_1k$value < -0.2] <- "-0.3 to -0.2"
chr3_1k$group[chr3_1k$value >=-0.2 & chr3_1k$value < -0.1] <- "-0.2 to -0.1"
chr3_1k$group[chr3_1k$value >=-0.1 & chr3_1k$value < -0] <- "-0.1 to 0"
chr3_1k$group[chr3_1k$value >=0 & chr3_1k$value < 0.1] <- "0 to 0.1"
chr3_1k$group[chr3_1k$value >=0.1 & chr3_1k$value < 0.2] <- "0.1 to 0.2"
chr3_1k$group[chr3_1k$value >=0.2 & chr3_1k$value < 0.3] <- "0.2 to 0.3"
chr3_1k$group[chr3_1k$value >=0.3 & chr3_1k$value < 0.4] <- "0.3 to 0.4"
chr3_1k$group[chr3_1k$value >=0.4 & chr3_1k$value < 0.5] <- "0.4 to 0.5"
chr3_1k$group[chr3_1k$value >=0.5 & chr3_1k$value < 0.6] <- "0.5 to 0.6"
chr3_1k$group[chr3_1k$value >=0.6 & chr3_1k$value < 0.7] <- "0.6 to 0.7"
chr3_1k$group[chr3_1k$value >=0.7 & chr3_1k$value < 0.8] <- "0.7 to 0.8"
chr3_1k$group[chr3_1k$value >=0.8 & chr3_1k$value < 0.9] <- "0.8 to 0.9"
chr3_1k$group[chr3_1k$value >=0.9 & chr3_1k$value < 1] <- "0.9 to 1"

chr3_1k$group <- factor(chr3_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr3_1k$group)
cor.group <- as.data.frame(table(chr3_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr3_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 3: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()



###### chr4 #####

load("/path/to/cis_correlations/chr4_1.Rdata")
load("/path/to/cis_correlations/chr4_2.Rdata")
load("/path/to/cis_correlations/chr4_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr4_1_dim <- lapply(chr4, FUN = getdim)
chr4_1_dim_total <- lapply(chr4_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr4_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr4_2_dim <- lapply(chr4_2, FUN = getdim)
chr4_2_dim_total <- lapply(chr4_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr4_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr4_3_dim <- lapply(chr4_3, FUN = getdim)
chr4_3_dim_total <- lapply(chr4_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr4_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr4_1 <- rbindlist(chr4)
dim(chr4_1)
chr4_2 <- rbindlist(chr4_2)
dim(chr4_2)
chr4_3 <- rbindlist(chr4_3)
dim(chr4_3)

chr4 <- rbind(chr4_1, chr4_2)
dim(chr4)
chr4 <- rbind(chr4, chr4_3)
dim(chr4)
F7_cis_dim <- list()
F7_cis_dim$chr4_F7 <- dim(chr4)

chr4_10k <- chr4[chr4$absdistance <= 10000]
dim(chr4_10k)
head(chr4_10k)
chr4_10k_sort <- chr4_10k[order(absdistance)]
head(chr4_10k_sort)

# split to positive and negative correlations
chr4_10k_pos <- chr4_10k_sort[chr4_10k_sort$value >= 0,]
dim(chr4_10k_pos)
chr4_10k_neg <- chr4_10k_sort[chr4_10k_sort$value < 0,]
dim(chr4_10k_neg)

# bin the positive correlations on distance
chr4_10k_pos$group <- as.numeric(cut2(chr4_10k_pos$absdistance, m=100))
summary(chr4_10k_pos$group)
chr4_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr4_10k_pos_valuemean <- aggregate(chr4_10k_pos[, 3], list(chr4_10k_pos$group), mean)
class(chr4_10k_pos_valuemean)
head(chr4_10k_pos_valuemean)
# calculate the median distance of each bin
chr4_10k_pos_mediandistance <- aggregate(chr4_10k_pos[, 4], list(chr4_10k_pos$group), median)
class(chr4_10k_pos_mediandistance)
head(chr4_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr4_10k_pos_df <- merge(chr4_10k_pos_valuemean, chr4_10k_pos_mediandistance, by=c("Group.1"))
head(chr4_10k_pos_df)

# bin the negative correlations on distance
chr4_10k_neg$group <- as.numeric(cut2(chr4_10k_neg$absdistance, m=100))
summary(chr4_10k_neg$group)
chr4_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr4_10k_neg_valuemean <- aggregate(chr4_10k_neg[, 3], list(chr4_10k_neg$group), mean)
class(chr4_10k_neg_valuemean)
head(chr4_10k_neg_valuemean)
# calculate the median distance of each bin
chr4_10k_neg_mediandistance <- aggregate(chr4_10k_neg[, 4], list(chr4_10k_neg$group), median)
class(chr4_10k_neg_mediandistance)
head(chr4_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr4_10k_neg_df <- merge(chr4_10k_neg_valuemean, chr4_10k_neg_mediandistance, by=c("Group.1"))
head(chr4_10k_neg_df)

# calculate the SD in each positive bin
chr4_10k_pos_sd <- aggregate(chr4_10k_pos[, 3], list(chr4_10k_pos$group), sd)
class(chr4_10k_pos_sd)
head(chr4_10k_pos_sd)
chr4_10k_pos_df <- merge(chr4_10k_pos_df, chr4_10k_pos_sd, by=c("Group.1"))
head(chr4_10k_pos_df)
chr4_10k_pos_df$n_in_group <- as.numeric(table(chr4_10k_pos$group))
head(chr4_10k_pos_df)
table(chr4_10k_pos$group)

# calculate the SD in each negative bin
chr4_10k_neg_sd <- aggregate(chr4_10k_neg[, 3], list(chr4_10k_neg$group), sd)
class(chr4_10k_neg_sd)
head(chr4_10k_neg_sd)
chr4_10k_neg_df <- merge(chr4_10k_neg_df, chr4_10k_neg_sd, by=c("Group.1"))
head(chr4_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr4_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr4_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr4_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr4_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr4_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 4: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr4_1k <- chr4_10k_sort[chr4_10k_sort$absdistance <= 1000,]
chr4_1k$group <- "group"
chr4_1k$group[chr4_1k$value >=-1 & chr4_1k$value < -0.9] <- "-1 to -0.9"
chr4_1k$group[chr4_1k$value >=-0.9 & chr4_1k$value < -0.8] <- "-0.9 to -0.8"
chr4_1k$group[chr4_1k$value >=-0.8 & chr4_1k$value < -0.7] <- "-0.8 to -0.7"
chr4_1k$group[chr4_1k$value >=-0.7 & chr4_1k$value < -0.6] <- "-0.7 to -0.6"
chr4_1k$group[chr4_1k$value >=-0.6 & chr4_1k$value < -0.5] <- "-0.6 to -0.5"
chr4_1k$group[chr4_1k$value >=-0.5 & chr4_1k$value < -0.4] <- "-0.5 to -0.4"
chr4_1k$group[chr4_1k$value >=-0.4 & chr4_1k$value < -0.3] <- "-0.4 to -0.3"
chr4_1k$group[chr4_1k$value >=-0.3 & chr4_1k$value < -0.2] <- "-0.3 to -0.2"
chr4_1k$group[chr4_1k$value >=-0.2 & chr4_1k$value < -0.1] <- "-0.2 to -0.1"
chr4_1k$group[chr4_1k$value >=-0.1 & chr4_1k$value < -0] <- "-0.1 to 0"
chr4_1k$group[chr4_1k$value >=0 & chr4_1k$value < 0.1] <- "0 to 0.1"
chr4_1k$group[chr4_1k$value >=0.1 & chr4_1k$value < 0.2] <- "0.1 to 0.2"
chr4_1k$group[chr4_1k$value >=0.2 & chr4_1k$value < 0.3] <- "0.2 to 0.3"
chr4_1k$group[chr4_1k$value >=0.3 & chr4_1k$value < 0.4] <- "0.3 to 0.4"
chr4_1k$group[chr4_1k$value >=0.4 & chr4_1k$value < 0.5] <- "0.4 to 0.5"
chr4_1k$group[chr4_1k$value >=0.5 & chr4_1k$value < 0.6] <- "0.5 to 0.6"
chr4_1k$group[chr4_1k$value >=0.6 & chr4_1k$value < 0.7] <- "0.6 to 0.7"
chr4_1k$group[chr4_1k$value >=0.7 & chr4_1k$value < 0.8] <- "0.7 to 0.8"
chr4_1k$group[chr4_1k$value >=0.8 & chr4_1k$value < 0.9] <- "0.8 to 0.9"
chr4_1k$group[chr4_1k$value >=0.9 & chr4_1k$value < 1] <- "0.9 to 1"

chr4_1k$group <- factor(chr4_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr4_1k$group)
cor.group <- as.data.frame(table(chr4_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr4_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 4: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()



###### chr5 #####

load("/path/to/cis_correlations/chr5_1.Rdata")
load("/path/to/cis_correlations/chr5_2.Rdata")
load("/path/to/cis_correlations/chr5_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr5_1_dim <- lapply(chr5, FUN = getdim)
chr5_1_dim_total <- lapply(chr5_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr5_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr5_2_dim <- lapply(chr5_2, FUN = getdim)
chr5_2_dim_total <- lapply(chr5_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr5_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr5_3_dim <- lapply(chr5_3, FUN = getdim)
chr5_3_dim_total <- lapply(chr5_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr5_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr5_1 <- rbindlist(chr5)
dim(chr5_1)
chr5_2 <- rbindlist(chr5_2)
dim(chr5_2)
chr5_3 <- rbindlist(chr5_3)
dim(chr5_3)

chr5 <- rbind(chr5_1, chr5_2)
dim(chr5)
chr5 <- rbind(chr5, chr5_3)
dim(chr5)
F7_cis_dim <- list()
F7_cis_dim$chr5_F7 <- dim(chr5)

chr5_10k <- chr5[chr5$absdistance <= 10000]
dim(chr5_10k)
head(chr5_10k)
chr5_10k_sort <- chr5_10k[order(absdistance)]
head(chr5_10k_sort)

# split to positive and negative correlations
chr5_10k_pos <- chr5_10k_sort[chr5_10k_sort$value >= 0,]
dim(chr5_10k_pos)
chr5_10k_neg <- chr5_10k_sort[chr5_10k_sort$value < 0,]
dim(chr5_10k_neg)

# bin the positive correlations on distance
chr5_10k_pos$group <- as.numeric(cut2(chr5_10k_pos$absdistance, m=100))
summary(chr5_10k_pos$group)
chr5_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr5_10k_pos_valuemean <- aggregate(chr5_10k_pos[, 3], list(chr5_10k_pos$group), mean)
class(chr5_10k_pos_valuemean)
head(chr5_10k_pos_valuemean)
# calculate the median distance of each bin
chr5_10k_pos_mediandistance <- aggregate(chr5_10k_pos[, 4], list(chr5_10k_pos$group), median)
class(chr5_10k_pos_mediandistance)
head(chr5_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr5_10k_pos_df <- merge(chr5_10k_pos_valuemean, chr5_10k_pos_mediandistance, by=c("Group.1"))
head(chr5_10k_pos_df)

# bin the negative correlations on distance
chr5_10k_neg$group <- as.numeric(cut2(chr5_10k_neg$absdistance, m=100))
summary(chr5_10k_neg$group)
chr5_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr5_10k_neg_valuemean <- aggregate(chr5_10k_neg[, 3], list(chr5_10k_neg$group), mean)
class(chr5_10k_neg_valuemean)
head(chr5_10k_neg_valuemean)
# calculate the median distance of each bin
chr5_10k_neg_mediandistance <- aggregate(chr5_10k_neg[, 4], list(chr5_10k_neg$group), median)
class(chr5_10k_neg_mediandistance)
head(chr5_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr5_10k_neg_df <- merge(chr5_10k_neg_valuemean, chr5_10k_neg_mediandistance, by=c("Group.1"))
head(chr5_10k_neg_df)

# calculate the SD in each positive bin
chr5_10k_pos_sd <- aggregate(chr5_10k_pos[, 3], list(chr5_10k_pos$group), sd)
class(chr5_10k_pos_sd)
head(chr5_10k_pos_sd)
chr5_10k_pos_df <- merge(chr5_10k_pos_df, chr5_10k_pos_sd, by=c("Group.1"))
head(chr5_10k_pos_df)
chr5_10k_pos_df$n_in_group <- as.numeric(table(chr5_10k_pos$group))
head(chr5_10k_pos_df)
table(chr5_10k_pos$group)

# calculate the SD in each negative bin
chr5_10k_neg_sd <- aggregate(chr5_10k_neg[, 3], list(chr5_10k_neg$group), sd)
class(chr5_10k_neg_sd)
head(chr5_10k_neg_sd)
chr5_10k_neg_df <- merge(chr5_10k_neg_df, chr5_10k_neg_sd, by=c("Group.1"))
head(chr5_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr5_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr5_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr5_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr5_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr5_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 5: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr5_1k <- chr5_10k_sort[chr5_10k_sort$absdistance <= 1000,]
chr5_1k$group <- "group"
chr5_1k$group[chr5_1k$value >=-1 & chr5_1k$value < -0.9] <- "-1 to -0.9"
chr5_1k$group[chr5_1k$value >=-0.9 & chr5_1k$value < -0.8] <- "-0.9 to -0.8"
chr5_1k$group[chr5_1k$value >=-0.8 & chr5_1k$value < -0.7] <- "-0.8 to -0.7"
chr5_1k$group[chr5_1k$value >=-0.7 & chr5_1k$value < -0.6] <- "-0.7 to -0.6"
chr5_1k$group[chr5_1k$value >=-0.6 & chr5_1k$value < -0.5] <- "-0.6 to -0.5"
chr5_1k$group[chr5_1k$value >=-0.5 & chr5_1k$value < -0.4] <- "-0.5 to -0.4"
chr5_1k$group[chr5_1k$value >=-0.4 & chr5_1k$value < -0.3] <- "-0.4 to -0.3"
chr5_1k$group[chr5_1k$value >=-0.3 & chr5_1k$value < -0.2] <- "-0.3 to -0.2"
chr5_1k$group[chr5_1k$value >=-0.2 & chr5_1k$value < -0.1] <- "-0.2 to -0.1"
chr5_1k$group[chr5_1k$value >=-0.1 & chr5_1k$value < -0] <- "-0.1 to 0"
chr5_1k$group[chr5_1k$value >=0 & chr5_1k$value < 0.1] <- "0 to 0.1"
chr5_1k$group[chr5_1k$value >=0.1 & chr5_1k$value < 0.2] <- "0.1 to 0.2"
chr5_1k$group[chr5_1k$value >=0.2 & chr5_1k$value < 0.3] <- "0.2 to 0.3"
chr5_1k$group[chr5_1k$value >=0.3 & chr5_1k$value < 0.4] <- "0.3 to 0.4"
chr5_1k$group[chr5_1k$value >=0.4 & chr5_1k$value < 0.5] <- "0.4 to 0.5"
chr5_1k$group[chr5_1k$value >=0.5 & chr5_1k$value < 0.6] <- "0.5 to 0.6"
chr5_1k$group[chr5_1k$value >=0.6 & chr5_1k$value < 0.7] <- "0.6 to 0.7"
chr5_1k$group[chr5_1k$value >=0.7 & chr5_1k$value < 0.8] <- "0.7 to 0.8"
chr5_1k$group[chr5_1k$value >=0.8 & chr5_1k$value < 0.9] <- "0.8 to 0.9"
chr5_1k$group[chr5_1k$value >=0.9 & chr5_1k$value < 1] <- "0.9 to 1"

chr5_1k$group <- factor(chr5_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr5_1k$group)
cor.group <- as.data.frame(table(chr5_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr5_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 5: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()



###### chr6 #####

load("/path/to/cis_correlations/chr6_1.Rdata")
load("/path/to/cis_correlations/chr6_2.Rdata")
load("/path/to/cis_correlations/chr6_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr6_1_dim <- lapply(chr6, FUN = getdim)
chr6_1_dim_total <- lapply(chr6_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr6_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr6_2_dim <- lapply(chr6_2, FUN = getdim)
chr6_2_dim_total <- lapply(chr6_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr6_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr6_3_dim <- lapply(chr6_3, FUN = getdim)
chr6_3_dim_total <- lapply(chr6_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr6_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr6_1 <- rbindlist(chr6)
dim(chr6_1)
chr6_2 <- rbindlist(chr6_2)
dim(chr6_2)
chr6_3 <- rbindlist(chr6_3)
dim(chr6_3)

chr6 <- rbind(chr6_1, chr6_2)
dim(chr6)
chr6 <- rbind(chr6, chr6_3)
dim(chr6)
F7_cis_dim <- list()
F7_cis_dim$chr6_F7 <- dim(chr6)

chr6_10k <- chr6[chr6$absdistance <= 10000]
dim(chr6_10k)
head(chr6_10k)
chr6_10k_sort <- chr6_10k[order(absdistance)]
head(chr6_10k_sort)

# split to positive and negative correlations
chr6_10k_pos <- chr6_10k_sort[chr6_10k_sort$value >= 0,]
dim(chr6_10k_pos)
chr6_10k_neg <- chr6_10k_sort[chr6_10k_sort$value < 0,]
dim(chr6_10k_neg)

# bin the positive correlations on distance
chr6_10k_pos$group <- as.numeric(cut2(chr6_10k_pos$absdistance, m=100))
summary(chr6_10k_pos$group)
chr6_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr6_10k_pos_valuemean <- aggregate(chr6_10k_pos[, 3], list(chr6_10k_pos$group), mean)
class(chr6_10k_pos_valuemean)
head(chr6_10k_pos_valuemean)
# calculate the median distance of each bin
chr6_10k_pos_mediandistance <- aggregate(chr6_10k_pos[, 4], list(chr6_10k_pos$group), median)
class(chr6_10k_pos_mediandistance)
head(chr6_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr6_10k_pos_df <- merge(chr6_10k_pos_valuemean, chr6_10k_pos_mediandistance, by=c("Group.1"))
head(chr6_10k_pos_df)

# bin the negative correlations on distance
chr6_10k_neg$group <- as.numeric(cut2(chr6_10k_neg$absdistance, m=100))
summary(chr6_10k_neg$group)
chr6_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr6_10k_neg_valuemean <- aggregate(chr6_10k_neg[, 3], list(chr6_10k_neg$group), mean)
class(chr6_10k_neg_valuemean)
head(chr6_10k_neg_valuemean)
# calculate the median distance of each bin
chr6_10k_neg_mediandistance <- aggregate(chr6_10k_neg[, 4], list(chr6_10k_neg$group), median)
class(chr6_10k_neg_mediandistance)
head(chr6_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr6_10k_neg_df <- merge(chr6_10k_neg_valuemean, chr6_10k_neg_mediandistance, by=c("Group.1"))
head(chr6_10k_neg_df)

# calculate the SD in each positive bin
chr6_10k_pos_sd <- aggregate(chr6_10k_pos[, 3], list(chr6_10k_pos$group), sd)
class(chr6_10k_pos_sd)
head(chr6_10k_pos_sd)
chr6_10k_pos_df <- merge(chr6_10k_pos_df, chr6_10k_pos_sd, by=c("Group.1"))
head(chr6_10k_pos_df)
chr6_10k_pos_df$n_in_group <- as.numeric(table(chr6_10k_pos$group))
head(chr6_10k_pos_df)
table(chr6_10k_pos$group)

# calculate the SD in each negative bin
chr6_10k_neg_sd <- aggregate(chr6_10k_neg[, 3], list(chr6_10k_neg$group), sd)
class(chr6_10k_neg_sd)
head(chr6_10k_neg_sd)
chr6_10k_neg_df <- merge(chr6_10k_neg_df, chr6_10k_neg_sd, by=c("Group.1"))
head(chr6_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr6_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr6_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr6_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr6_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr6_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 6: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr6_1k <- chr6_10k_sort[chr6_10k_sort$absdistance <= 1000,]
chr6_1k$group <- "group"
chr6_1k$group[chr6_1k$value >=-1 & chr6_1k$value < -0.9] <- "-1 to -0.9"
chr6_1k$group[chr6_1k$value >=-0.9 & chr6_1k$value < -0.8] <- "-0.9 to -0.8"
chr6_1k$group[chr6_1k$value >=-0.8 & chr6_1k$value < -0.7] <- "-0.8 to -0.7"
chr6_1k$group[chr6_1k$value >=-0.7 & chr6_1k$value < -0.6] <- "-0.7 to -0.6"
chr6_1k$group[chr6_1k$value >=-0.6 & chr6_1k$value < -0.5] <- "-0.6 to -0.5"
chr6_1k$group[chr6_1k$value >=-0.5 & chr6_1k$value < -0.4] <- "-0.5 to -0.4"
chr6_1k$group[chr6_1k$value >=-0.4 & chr6_1k$value < -0.3] <- "-0.4 to -0.3"
chr6_1k$group[chr6_1k$value >=-0.3 & chr6_1k$value < -0.2] <- "-0.3 to -0.2"
chr6_1k$group[chr6_1k$value >=-0.2 & chr6_1k$value < -0.1] <- "-0.2 to -0.1"
chr6_1k$group[chr6_1k$value >=-0.1 & chr6_1k$value < -0] <- "-0.1 to 0"
chr6_1k$group[chr6_1k$value >=0 & chr6_1k$value < 0.1] <- "0 to 0.1"
chr6_1k$group[chr6_1k$value >=0.1 & chr6_1k$value < 0.2] <- "0.1 to 0.2"
chr6_1k$group[chr6_1k$value >=0.2 & chr6_1k$value < 0.3] <- "0.2 to 0.3"
chr6_1k$group[chr6_1k$value >=0.3 & chr6_1k$value < 0.4] <- "0.3 to 0.4"
chr6_1k$group[chr6_1k$value >=0.4 & chr6_1k$value < 0.5] <- "0.4 to 0.5"
chr6_1k$group[chr6_1k$value >=0.5 & chr6_1k$value < 0.6] <- "0.5 to 0.6"
chr6_1k$group[chr6_1k$value >=0.6 & chr6_1k$value < 0.7] <- "0.6 to 0.7"
chr6_1k$group[chr6_1k$value >=0.7 & chr6_1k$value < 0.8] <- "0.7 to 0.8"
chr6_1k$group[chr6_1k$value >=0.8 & chr6_1k$value < 0.9] <- "0.8 to 0.9"
chr6_1k$group[chr6_1k$value >=0.9 & chr6_1k$value < 1] <- "0.9 to 1"

chr6_1k$group <- factor(chr6_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr6_1k$group)
cor.group <- as.data.frame(table(chr6_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr6_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 6: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr7 #####

load("/path/to/cis_correlations/chr7_1.Rdata")
load("/path/to/cis_correlations/chr7_2.Rdata")
load("/path/to/cis_correlations/chr7_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr7_1_dim <- lapply(chr7, FUN = getdim)
chr7_1_dim_total <- lapply(chr7_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr7_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr7_2_dim <- lapply(chr7_2, FUN = getdim)
chr7_2_dim_total <- lapply(chr7_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr7_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr7_3_dim <- lapply(chr7_3, FUN = getdim)
chr7_3_dim_total <- lapply(chr7_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr7_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr7_1 <- rbindlist(chr7)
dim(chr7_1)
chr7_2 <- rbindlist(chr7_2)
dim(chr7_2)
chr7_3 <- rbindlist(chr7_3)
dim(chr7_3)

chr7 <- rbind(chr7_1, chr7_2)
dim(chr7)
chr7 <- rbind(chr7, chr7_3)
dim(chr7)
F7_cis_dim <- list()
F7_cis_dim$chr7_F7 <- dim(chr7)

chr7_10k <- chr7[chr7$absdistance <= 10000]
dim(chr7_10k)
head(chr7_10k)
chr7_10k_sort <- chr7_10k[order(absdistance)]
head(chr7_10k_sort)

# split to positive and negative correlations
chr7_10k_pos <- chr7_10k_sort[chr7_10k_sort$value >= 0,]
dim(chr7_10k_pos)
chr7_10k_neg <- chr7_10k_sort[chr7_10k_sort$value < 0,]
dim(chr7_10k_neg)

# bin the positive correlations on distance
chr7_10k_pos$group <- as.numeric(cut2(chr7_10k_pos$absdistance, m=100))
summary(chr7_10k_pos$group)
chr7_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr7_10k_pos_valuemean <- aggregate(chr7_10k_pos[, 3], list(chr7_10k_pos$group), mean)
class(chr7_10k_pos_valuemean)
head(chr7_10k_pos_valuemean)
# calculate the median distance of each bin
chr7_10k_pos_mediandistance <- aggregate(chr7_10k_pos[, 4], list(chr7_10k_pos$group), median)
class(chr7_10k_pos_mediandistance)
head(chr7_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr7_10k_pos_df <- merge(chr7_10k_pos_valuemean, chr7_10k_pos_mediandistance, by=c("Group.1"))
head(chr7_10k_pos_df)

# bin the negative correlations on distance
chr7_10k_neg$group <- as.numeric(cut2(chr7_10k_neg$absdistance, m=100))
summary(chr7_10k_neg$group)
chr7_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr7_10k_neg_valuemean <- aggregate(chr7_10k_neg[, 3], list(chr7_10k_neg$group), mean)
class(chr7_10k_neg_valuemean)
head(chr7_10k_neg_valuemean)
# calculate the median distance of each bin
chr7_10k_neg_mediandistance <- aggregate(chr7_10k_neg[, 4], list(chr7_10k_neg$group), median)
class(chr7_10k_neg_mediandistance)
head(chr7_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr7_10k_neg_df <- merge(chr7_10k_neg_valuemean, chr7_10k_neg_mediandistance, by=c("Group.1"))
head(chr7_10k_neg_df)

# calculate the SD in each positive bin
chr7_10k_pos_sd <- aggregate(chr7_10k_pos[, 3], list(chr7_10k_pos$group), sd)
class(chr7_10k_pos_sd)
head(chr7_10k_pos_sd)
chr7_10k_pos_df <- merge(chr7_10k_pos_df, chr7_10k_pos_sd, by=c("Group.1"))
head(chr7_10k_pos_df)
chr7_10k_pos_df$n_in_group <- as.numeric(table(chr7_10k_pos$group))
head(chr7_10k_pos_df)
table(chr7_10k_pos$group)

# calculate the SD in each negative bin
chr7_10k_neg_sd <- aggregate(chr7_10k_neg[, 3], list(chr7_10k_neg$group), sd)
class(chr7_10k_neg_sd)
head(chr7_10k_neg_sd)
chr7_10k_neg_df <- merge(chr7_10k_neg_df, chr7_10k_neg_sd, by=c("Group.1"))
head(chr7_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr7_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr7_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr7_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr7_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr7_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 7: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr7_1k <- chr7_10k_sort[chr7_10k_sort$absdistance <= 1000,]
chr7_1k$group <- "group"
chr7_1k$group[chr7_1k$value >=-1 & chr7_1k$value < -0.9] <- "-1 to -0.9"
chr7_1k$group[chr7_1k$value >=-0.9 & chr7_1k$value < -0.8] <- "-0.9 to -0.8"
chr7_1k$group[chr7_1k$value >=-0.8 & chr7_1k$value < -0.7] <- "-0.8 to -0.7"
chr7_1k$group[chr7_1k$value >=-0.7 & chr7_1k$value < -0.6] <- "-0.7 to -0.6"
chr7_1k$group[chr7_1k$value >=-0.6 & chr7_1k$value < -0.5] <- "-0.6 to -0.5"
chr7_1k$group[chr7_1k$value >=-0.5 & chr7_1k$value < -0.4] <- "-0.5 to -0.4"
chr7_1k$group[chr7_1k$value >=-0.4 & chr7_1k$value < -0.3] <- "-0.4 to -0.3"
chr7_1k$group[chr7_1k$value >=-0.3 & chr7_1k$value < -0.2] <- "-0.3 to -0.2"
chr7_1k$group[chr7_1k$value >=-0.2 & chr7_1k$value < -0.1] <- "-0.2 to -0.1"
chr7_1k$group[chr7_1k$value >=-0.1 & chr7_1k$value < -0] <- "-0.1 to 0"
chr7_1k$group[chr7_1k$value >=0 & chr7_1k$value < 0.1] <- "0 to 0.1"
chr7_1k$group[chr7_1k$value >=0.1 & chr7_1k$value < 0.2] <- "0.1 to 0.2"
chr7_1k$group[chr7_1k$value >=0.2 & chr7_1k$value < 0.3] <- "0.2 to 0.3"
chr7_1k$group[chr7_1k$value >=0.3 & chr7_1k$value < 0.4] <- "0.3 to 0.4"
chr7_1k$group[chr7_1k$value >=0.4 & chr7_1k$value < 0.5] <- "0.4 to 0.5"
chr7_1k$group[chr7_1k$value >=0.5 & chr7_1k$value < 0.6] <- "0.5 to 0.6"
chr7_1k$group[chr7_1k$value >=0.6 & chr7_1k$value < 0.7] <- "0.6 to 0.7"
chr7_1k$group[chr7_1k$value >=0.7 & chr7_1k$value < 0.8] <- "0.7 to 0.8"
chr7_1k$group[chr7_1k$value >=0.8 & chr7_1k$value < 0.9] <- "0.8 to 0.9"
chr7_1k$group[chr7_1k$value >=0.9 & chr7_1k$value < 1] <- "0.9 to 1"

chr7_1k$group <- factor(chr7_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr7_1k$group)
cor.group <- as.data.frame(table(chr7_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr7_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 7: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr8 #####

load("/path/to/cis_correlations/chr8_1.Rdata")
load("/path/to/cis_correlations/chr8_2.Rdata")
load("/path/to/cis_correlations/chr8_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr8_1_dim <- lapply(chr8, FUN = getdim)
chr8_1_dim_total <- lapply(chr8_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr8_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr8_2_dim <- lapply(chr8_2, FUN = getdim)
chr8_2_dim_total <- lapply(chr8_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr8_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr8_3_dim <- lapply(chr8_3, FUN = getdim)
chr8_3_dim_total <- lapply(chr8_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr8_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr8_1 <- rbindlist(chr8)
dim(chr8_1)
chr8_2 <- rbindlist(chr8_2)
dim(chr8_2)
chr8_3 <- rbindlist(chr8_3)
dim(chr8_3)

chr8 <- rbind(chr8_1, chr8_2)
dim(chr8)
chr8 <- rbind(chr8, chr8_3)
dim(chr8)
F7_cis_dim <- list()
F7_cis_dim$chr8_F7 <- dim(chr8)

chr8_10k <- chr8[chr8$absdistance <= 10000]
dim(chr8_10k)
head(chr8_10k)
chr8_10k_sort <- chr8_10k[order(absdistance)]
head(chr8_10k_sort)

# split to positive and negative correlations
chr8_10k_pos <- chr8_10k_sort[chr8_10k_sort$value >= 0,]
dim(chr8_10k_pos)
chr8_10k_neg <- chr8_10k_sort[chr8_10k_sort$value < 0,]
dim(chr8_10k_neg)

# bin the positive correlations on distance
chr8_10k_pos$group <- as.numeric(cut2(chr8_10k_pos$absdistance, m=100))
summary(chr8_10k_pos$group)
chr8_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr8_10k_pos_valuemean <- aggregate(chr8_10k_pos[, 3], list(chr8_10k_pos$group), mean)
class(chr8_10k_pos_valuemean)
head(chr8_10k_pos_valuemean)
# calculate the median distance of each bin
chr8_10k_pos_mediandistance <- aggregate(chr8_10k_pos[, 4], list(chr8_10k_pos$group), median)
class(chr8_10k_pos_mediandistance)
head(chr8_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr8_10k_pos_df <- merge(chr8_10k_pos_valuemean, chr8_10k_pos_mediandistance, by=c("Group.1"))
head(chr8_10k_pos_df)

# bin the negative correlations on distance
chr8_10k_neg$group <- as.numeric(cut2(chr8_10k_neg$absdistance, m=100))
summary(chr8_10k_neg$group)
chr8_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr8_10k_neg_valuemean <- aggregate(chr8_10k_neg[, 3], list(chr8_10k_neg$group), mean)
class(chr8_10k_neg_valuemean)
head(chr8_10k_neg_valuemean)
# calculate the median distance of each bin
chr8_10k_neg_mediandistance <- aggregate(chr8_10k_neg[, 4], list(chr8_10k_neg$group), median)
class(chr8_10k_neg_mediandistance)
head(chr8_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr8_10k_neg_df <- merge(chr8_10k_neg_valuemean, chr8_10k_neg_mediandistance, by=c("Group.1"))
head(chr8_10k_neg_df)

# calculate the SD in each positive bin
chr8_10k_pos_sd <- aggregate(chr8_10k_pos[, 3], list(chr8_10k_pos$group), sd)
class(chr8_10k_pos_sd)
head(chr8_10k_pos_sd)
chr8_10k_pos_df <- merge(chr8_10k_pos_df, chr8_10k_pos_sd, by=c("Group.1"))
head(chr8_10k_pos_df)
chr8_10k_pos_df$n_in_group <- as.numeric(table(chr8_10k_pos$group))
head(chr8_10k_pos_df)
table(chr8_10k_pos$group)

# calculate the SD in each negative bin
chr8_10k_neg_sd <- aggregate(chr8_10k_neg[, 3], list(chr8_10k_neg$group), sd)
class(chr8_10k_neg_sd)
head(chr8_10k_neg_sd)
chr8_10k_neg_df <- merge(chr8_10k_neg_df, chr8_10k_neg_sd, by=c("Group.1"))
head(chr8_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr8_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr8_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr8_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr8_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr8_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 8: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr8_1k <- chr8_10k_sort[chr8_10k_sort$absdistance <= 1000,]
chr8_1k$group <- "group"
chr8_1k$group[chr8_1k$value >=-1 & chr8_1k$value < -0.9] <- "-1 to -0.9"
chr8_1k$group[chr8_1k$value >=-0.9 & chr8_1k$value < -0.8] <- "-0.9 to -0.8"
chr8_1k$group[chr8_1k$value >=-0.8 & chr8_1k$value < -0.7] <- "-0.8 to -0.7"
chr8_1k$group[chr8_1k$value >=-0.7 & chr8_1k$value < -0.6] <- "-0.7 to -0.6"
chr8_1k$group[chr8_1k$value >=-0.6 & chr8_1k$value < -0.5] <- "-0.6 to -0.5"
chr8_1k$group[chr8_1k$value >=-0.5 & chr8_1k$value < -0.4] <- "-0.5 to -0.4"
chr8_1k$group[chr8_1k$value >=-0.4 & chr8_1k$value < -0.3] <- "-0.4 to -0.3"
chr8_1k$group[chr8_1k$value >=-0.3 & chr8_1k$value < -0.2] <- "-0.3 to -0.2"
chr8_1k$group[chr8_1k$value >=-0.2 & chr8_1k$value < -0.1] <- "-0.2 to -0.1"
chr8_1k$group[chr8_1k$value >=-0.1 & chr8_1k$value < -0] <- "-0.1 to 0"
chr8_1k$group[chr8_1k$value >=0 & chr8_1k$value < 0.1] <- "0 to 0.1"
chr8_1k$group[chr8_1k$value >=0.1 & chr8_1k$value < 0.2] <- "0.1 to 0.2"
chr8_1k$group[chr8_1k$value >=0.2 & chr8_1k$value < 0.3] <- "0.2 to 0.3"
chr8_1k$group[chr8_1k$value >=0.3 & chr8_1k$value < 0.4] <- "0.3 to 0.4"
chr8_1k$group[chr8_1k$value >=0.4 & chr8_1k$value < 0.5] <- "0.4 to 0.5"
chr8_1k$group[chr8_1k$value >=0.5 & chr8_1k$value < 0.6] <- "0.5 to 0.6"
chr8_1k$group[chr8_1k$value >=0.6 & chr8_1k$value < 0.7] <- "0.6 to 0.7"
chr8_1k$group[chr8_1k$value >=0.7 & chr8_1k$value < 0.8] <- "0.7 to 0.8"
chr8_1k$group[chr8_1k$value >=0.8 & chr8_1k$value < 0.9] <- "0.8 to 0.9"
chr8_1k$group[chr8_1k$value >=0.9 & chr8_1k$value < 1] <- "0.9 to 1"

chr8_1k$group <- factor(chr8_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr8_1k$group)
cor.group <- as.data.frame(table(chr8_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr8_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 8: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr9 #####

load("/path/to/cis_correlations/chr9_1.Rdata")
load("/path/to/cis_correlations/chr9_2.Rdata")
load("/path/to/cis_correlations/chr9_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr9_1_dim <- lapply(chr9, FUN = getdim)
chr9_1_dim_total <- lapply(chr9_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr9_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr9_2_dim <- lapply(chr9_2, FUN = getdim)
chr9_2_dim_total <- lapply(chr9_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr9_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr9_3_dim <- lapply(chr9_3, FUN = getdim)
chr9_3_dim_total <- lapply(chr9_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr9_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr9_1 <- rbindlist(chr9)
dim(chr9_1)
chr9_2 <- rbindlist(chr9_2)
dim(chr9_2)
chr9_3 <- rbindlist(chr9_3)
dim(chr9_3)

chr9 <- rbind(chr9_1, chr9_2)
dim(chr9)
chr9 <- rbind(chr9, chr9_3)
dim(chr9)
F7_cis_dim <- list()
F7_cis_dim$chr9_F7 <- dim(chr9)

chr9_10k <- chr9[chr9$absdistance <= 10000]
dim(chr9_10k)
head(chr9_10k)
chr9_10k_sort <- chr9_10k[order(absdistance)]
head(chr9_10k_sort)

# split to positive and negative correlations
chr9_10k_pos <- chr9_10k_sort[chr9_10k_sort$value >= 0,]
dim(chr9_10k_pos)
chr9_10k_neg <- chr9_10k_sort[chr9_10k_sort$value < 0,]
dim(chr9_10k_neg)

# bin the positive correlations on distance
chr9_10k_pos$group <- as.numeric(cut2(chr9_10k_pos$absdistance, m=100))
summary(chr9_10k_pos$group)
chr9_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr9_10k_pos_valuemean <- aggregate(chr9_10k_pos[, 3], list(chr9_10k_pos$group), mean)
class(chr9_10k_pos_valuemean)
head(chr9_10k_pos_valuemean)
# calculate the median distance of each bin
chr9_10k_pos_mediandistance <- aggregate(chr9_10k_pos[, 4], list(chr9_10k_pos$group), median)
class(chr9_10k_pos_mediandistance)
head(chr9_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr9_10k_pos_df <- merge(chr9_10k_pos_valuemean, chr9_10k_pos_mediandistance, by=c("Group.1"))
head(chr9_10k_pos_df)

# bin the negative correlations on distance
chr9_10k_neg$group <- as.numeric(cut2(chr9_10k_neg$absdistance, m=100))
summary(chr9_10k_neg$group)
chr9_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr9_10k_neg_valuemean <- aggregate(chr9_10k_neg[, 3], list(chr9_10k_neg$group), mean)
class(chr9_10k_neg_valuemean)
head(chr9_10k_neg_valuemean)
# calculate the median distance of each bin
chr9_10k_neg_mediandistance <- aggregate(chr9_10k_neg[, 4], list(chr9_10k_neg$group), median)
class(chr9_10k_neg_mediandistance)
head(chr9_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr9_10k_neg_df <- merge(chr9_10k_neg_valuemean, chr9_10k_neg_mediandistance, by=c("Group.1"))
head(chr9_10k_neg_df)

# calculate the SD in each positive bin
chr9_10k_pos_sd <- aggregate(chr9_10k_pos[, 3], list(chr9_10k_pos$group), sd)
class(chr9_10k_pos_sd)
head(chr9_10k_pos_sd)
chr9_10k_pos_df <- merge(chr9_10k_pos_df, chr9_10k_pos_sd, by=c("Group.1"))
head(chr9_10k_pos_df)
chr9_10k_pos_df$n_in_group <- as.numeric(table(chr9_10k_pos$group))
head(chr9_10k_pos_df)
table(chr9_10k_pos$group)

# calculate the SD in each negative bin
chr9_10k_neg_sd <- aggregate(chr9_10k_neg[, 3], list(chr9_10k_neg$group), sd)
class(chr9_10k_neg_sd)
head(chr9_10k_neg_sd)
chr9_10k_neg_df <- merge(chr9_10k_neg_df, chr9_10k_neg_sd, by=c("Group.1"))
head(chr9_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr9_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr9_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr9_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr9_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr9_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 9: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr9_1k <- chr9_10k_sort[chr9_10k_sort$absdistance <= 1000,]
chr9_1k$group <- "group"
chr9_1k$group[chr9_1k$value >=-1 & chr9_1k$value < -0.9] <- "-1 to -0.9"
chr9_1k$group[chr9_1k$value >=-0.9 & chr9_1k$value < -0.8] <- "-0.9 to -0.8"
chr9_1k$group[chr9_1k$value >=-0.8 & chr9_1k$value < -0.7] <- "-0.8 to -0.7"
chr9_1k$group[chr9_1k$value >=-0.7 & chr9_1k$value < -0.6] <- "-0.7 to -0.6"
chr9_1k$group[chr9_1k$value >=-0.6 & chr9_1k$value < -0.5] <- "-0.6 to -0.5"
chr9_1k$group[chr9_1k$value >=-0.5 & chr9_1k$value < -0.4] <- "-0.5 to -0.4"
chr9_1k$group[chr9_1k$value >=-0.4 & chr9_1k$value < -0.3] <- "-0.4 to -0.3"
chr9_1k$group[chr9_1k$value >=-0.3 & chr9_1k$value < -0.2] <- "-0.3 to -0.2"
chr9_1k$group[chr9_1k$value >=-0.2 & chr9_1k$value < -0.1] <- "-0.2 to -0.1"
chr9_1k$group[chr9_1k$value >=-0.1 & chr9_1k$value < -0] <- "-0.1 to 0"
chr9_1k$group[chr9_1k$value >=0 & chr9_1k$value < 0.1] <- "0 to 0.1"
chr9_1k$group[chr9_1k$value >=0.1 & chr9_1k$value < 0.2] <- "0.1 to 0.2"
chr9_1k$group[chr9_1k$value >=0.2 & chr9_1k$value < 0.3] <- "0.2 to 0.3"
chr9_1k$group[chr9_1k$value >=0.3 & chr9_1k$value < 0.4] <- "0.3 to 0.4"
chr9_1k$group[chr9_1k$value >=0.4 & chr9_1k$value < 0.5] <- "0.4 to 0.5"
chr9_1k$group[chr9_1k$value >=0.5 & chr9_1k$value < 0.6] <- "0.5 to 0.6"
chr9_1k$group[chr9_1k$value >=0.6 & chr9_1k$value < 0.7] <- "0.6 to 0.7"
chr9_1k$group[chr9_1k$value >=0.7 & chr9_1k$value < 0.8] <- "0.7 to 0.8"
chr9_1k$group[chr9_1k$value >=0.8 & chr9_1k$value < 0.9] <- "0.8 to 0.9"
chr9_1k$group[chr9_1k$value >=0.9 & chr9_1k$value < 1] <- "0.9 to 1"

chr9_1k$group <- factor(chr9_1k$group, 
                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                 "-0.2 to -0.1", "-0.1 to 0",
                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                 "0.9 to 1"))


table(chr9_1k$group)
cor.group <- as.data.frame(table(chr9_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr9_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 9: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()


###### chr10 #####

load("/path/to/cis_correlations/chr10_1.Rdata")
load("/path/to/cis_correlations/chr10_2.Rdata")
load("/path/to/cis_correlations/chr10_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr10_1_dim <- lapply(chr10, FUN = getdim)
chr10_1_dim_total <- lapply(chr10_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr10_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr10_2_dim <- lapply(chr10_2, FUN = getdim)
chr10_2_dim_total <- lapply(chr10_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr10_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr10_3_dim <- lapply(chr10_3, FUN = getdim)
chr10_3_dim_total <- lapply(chr10_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr10_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr10_1 <- rbindlist(chr10)
dim(chr10_1)
chr10_2 <- rbindlist(chr10_2)
dim(chr10_2)
chr10_3 <- rbindlist(chr10_3)
dim(chr10_3)

chr10 <- rbind(chr10_1, chr10_2)
dim(chr10)
chr10 <- rbind(chr10, chr10_3)
dim(chr10)
F7_cis_dim <- list()
F7_cis_dim$chr10_F7 <- dim(chr10)

chr10_10k <- chr10[chr10$absdistance <= 10000]
dim(chr10_10k)
head(chr10_10k)
chr10_10k_sort <- chr10_10k[order(absdistance)]
head(chr10_10k_sort)

# split to positive and negative correlations
chr10_10k_pos <- chr10_10k_sort[chr10_10k_sort$value >= 0,]
dim(chr10_10k_pos)
chr10_10k_neg <- chr10_10k_sort[chr10_10k_sort$value < 0,]
dim(chr10_10k_neg)

# bin the positive correlations on distance
chr10_10k_pos$group <- as.numeric(cut2(chr10_10k_pos$absdistance, m=100))
summary(chr10_10k_pos$group)
chr10_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr10_10k_pos_valuemean <- aggregate(chr10_10k_pos[, 3], list(chr10_10k_pos$group), mean)
class(chr10_10k_pos_valuemean)
head(chr10_10k_pos_valuemean)
# calculate the median distance of each bin
chr10_10k_pos_mediandistance <- aggregate(chr10_10k_pos[, 4], list(chr10_10k_pos$group), median)
class(chr10_10k_pos_mediandistance)
head(chr10_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr10_10k_pos_df <- merge(chr10_10k_pos_valuemean, chr10_10k_pos_mediandistance, by=c("Group.1"))
head(chr10_10k_pos_df)

# bin the negative correlations on distance
chr10_10k_neg$group <- as.numeric(cut2(chr10_10k_neg$absdistance, m=100))
summary(chr10_10k_neg$group)
chr10_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr10_10k_neg_valuemean <- aggregate(chr10_10k_neg[, 3], list(chr10_10k_neg$group), mean)
class(chr10_10k_neg_valuemean)
head(chr10_10k_neg_valuemean)
# calculate the median distance of each bin
chr10_10k_neg_mediandistance <- aggregate(chr10_10k_neg[, 4], list(chr10_10k_neg$group), median)
class(chr10_10k_neg_mediandistance)
head(chr10_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr10_10k_neg_df <- merge(chr10_10k_neg_valuemean, chr10_10k_neg_mediandistance, by=c("Group.1"))
head(chr10_10k_neg_df)

# calculate the SD in each positive bin
chr10_10k_pos_sd <- aggregate(chr10_10k_pos[, 3], list(chr10_10k_pos$group), sd)
class(chr10_10k_pos_sd)
head(chr10_10k_pos_sd)
chr10_10k_pos_df <- merge(chr10_10k_pos_df, chr10_10k_pos_sd, by=c("Group.1"))
head(chr10_10k_pos_df)
chr10_10k_pos_df$n_in_group <- as.numeric(table(chr10_10k_pos$group))
head(chr10_10k_pos_df)
table(chr10_10k_pos$group)

# calculate the SD in each negative bin
chr10_10k_neg_sd <- aggregate(chr10_10k_neg[, 3], list(chr10_10k_neg$group), sd)
class(chr10_10k_neg_sd)
head(chr10_10k_neg_sd)
chr10_10k_neg_df <- merge(chr10_10k_neg_df, chr10_10k_neg_sd, by=c("Group.1"))
head(chr10_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr10_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr10_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr10_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr10_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr10_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 10: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr10_1k <- chr10_10k_sort[chr10_10k_sort$absdistance <= 1000,]
chr10_1k$group <- "group"
chr10_1k$group[chr10_1k$value >=-1 & chr10_1k$value < -0.9] <- "-1 to -0.9"
chr10_1k$group[chr10_1k$value >=-0.9 & chr10_1k$value < -0.8] <- "-0.9 to -0.8"
chr10_1k$group[chr10_1k$value >=-0.8 & chr10_1k$value < -0.7] <- "-0.8 to -0.7"
chr10_1k$group[chr10_1k$value >=-0.7 & chr10_1k$value < -0.6] <- "-0.7 to -0.6"
chr10_1k$group[chr10_1k$value >=-0.6 & chr10_1k$value < -0.5] <- "-0.6 to -0.5"
chr10_1k$group[chr10_1k$value >=-0.5 & chr10_1k$value < -0.4] <- "-0.5 to -0.4"
chr10_1k$group[chr10_1k$value >=-0.4 & chr10_1k$value < -0.3] <- "-0.4 to -0.3"
chr10_1k$group[chr10_1k$value >=-0.3 & chr10_1k$value < -0.2] <- "-0.3 to -0.2"
chr10_1k$group[chr10_1k$value >=-0.2 & chr10_1k$value < -0.1] <- "-0.2 to -0.1"
chr10_1k$group[chr10_1k$value >=-0.1 & chr10_1k$value < -0] <- "-0.1 to 0"
chr10_1k$group[chr10_1k$value >=0 & chr10_1k$value < 0.1] <- "0 to 0.1"
chr10_1k$group[chr10_1k$value >=0.1 & chr10_1k$value < 0.2] <- "0.1 to 0.2"
chr10_1k$group[chr10_1k$value >=0.2 & chr10_1k$value < 0.3] <- "0.2 to 0.3"
chr10_1k$group[chr10_1k$value >=0.3 & chr10_1k$value < 0.4] <- "0.3 to 0.4"
chr10_1k$group[chr10_1k$value >=0.4 & chr10_1k$value < 0.5] <- "0.4 to 0.5"
chr10_1k$group[chr10_1k$value >=0.5 & chr10_1k$value < 0.6] <- "0.5 to 0.6"
chr10_1k$group[chr10_1k$value >=0.6 & chr10_1k$value < 0.7] <- "0.6 to 0.7"
chr10_1k$group[chr10_1k$value >=0.7 & chr10_1k$value < 0.8] <- "0.7 to 0.8"
chr10_1k$group[chr10_1k$value >=0.8 & chr10_1k$value < 0.9] <- "0.8 to 0.9"
chr10_1k$group[chr10_1k$value >=0.9 & chr10_1k$value < 1] <- "0.9 to 1"

chr10_1k$group <- factor(chr10_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr10_1k$group)
cor.group <- as.data.frame(table(chr10_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr10_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 10: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()





###### chr11 #####

load("/path/to/cis_correlations/chr11_1.Rdata")
load("/path/to/cis_correlations/chr11_2.Rdata")
load("/path/to/cis_correlations/chr11_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr11_1_dim <- lapply(chr11, FUN = getdim)
chr11_1_dim_total <- lapply(chr11_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr11_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr11_2_dim <- lapply(chr11_2, FUN = getdim)
chr11_2_dim_total <- lapply(chr11_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr11_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr11_3_dim <- lapply(chr11_3, FUN = getdim)
chr11_3_dim_total <- lapply(chr11_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr11_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr11_1 <- rbindlist(chr11)
dim(chr11_1)
chr11_2 <- rbindlist(chr11_2)
dim(chr11_2)
chr11_3 <- rbindlist(chr11_3)
dim(chr11_3)

chr11 <- rbind(chr11_1, chr11_2)
dim(chr11)
chr11 <- rbind(chr11, chr11_3)
dim(chr11)
F7_cis_dim <- list()
F7_cis_dim$chr11_F7 <- dim(chr11)

chr11_10k <- chr11[chr11$absdistance <= 10000]
dim(chr11_10k)
head(chr11_10k)
chr11_10k_sort <- chr11_10k[order(absdistance)]
head(chr11_10k_sort)

# split to positive and negative correlations
chr11_10k_pos <- chr11_10k_sort[chr11_10k_sort$value >= 0,]
dim(chr11_10k_pos)
chr11_10k_neg <- chr11_10k_sort[chr11_10k_sort$value < 0,]
dim(chr11_10k_neg)

# bin the positive correlations on distance
chr11_10k_pos$group <- as.numeric(cut2(chr11_10k_pos$absdistance, m=100))
summary(chr11_10k_pos$group)
chr11_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr11_10k_pos_valuemean <- aggregate(chr11_10k_pos[, 3], list(chr11_10k_pos$group), mean)
class(chr11_10k_pos_valuemean)
head(chr11_10k_pos_valuemean)
# calculate the median distance of each bin
chr11_10k_pos_mediandistance <- aggregate(chr11_10k_pos[, 4], list(chr11_10k_pos$group), median)
class(chr11_10k_pos_mediandistance)
head(chr11_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr11_10k_pos_df <- merge(chr11_10k_pos_valuemean, chr11_10k_pos_mediandistance, by=c("Group.1"))
head(chr11_10k_pos_df)

# bin the negative correlations on distance
chr11_10k_neg$group <- as.numeric(cut2(chr11_10k_neg$absdistance, m=100))
summary(chr11_10k_neg$group)
chr11_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr11_10k_neg_valuemean <- aggregate(chr11_10k_neg[, 3], list(chr11_10k_neg$group), mean)
class(chr11_10k_neg_valuemean)
head(chr11_10k_neg_valuemean)
# calculate the median distance of each bin
chr11_10k_neg_mediandistance <- aggregate(chr11_10k_neg[, 4], list(chr11_10k_neg$group), median)
class(chr11_10k_neg_mediandistance)
head(chr11_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr11_10k_neg_df <- merge(chr11_10k_neg_valuemean, chr11_10k_neg_mediandistance, by=c("Group.1"))
head(chr11_10k_neg_df)

# calculate the SD in each positive bin
chr11_10k_pos_sd <- aggregate(chr11_10k_pos[, 3], list(chr11_10k_pos$group), sd)
class(chr11_10k_pos_sd)
head(chr11_10k_pos_sd)
chr11_10k_pos_df <- merge(chr11_10k_pos_df, chr11_10k_pos_sd, by=c("Group.1"))
head(chr11_10k_pos_df)
chr11_10k_pos_df$n_in_group <- as.numeric(table(chr11_10k_pos$group))
head(chr11_10k_pos_df)
table(chr11_10k_pos$group)

# calculate the SD in each negative bin
chr11_10k_neg_sd <- aggregate(chr11_10k_neg[, 3], list(chr11_10k_neg$group), sd)
class(chr11_10k_neg_sd)
head(chr11_10k_neg_sd)
chr11_10k_neg_df <- merge(chr11_10k_neg_df, chr11_10k_neg_sd, by=c("Group.1"))
head(chr11_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr11_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr11_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr11_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr11_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr11_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 11: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr11_1k <- chr11_10k_sort[chr11_10k_sort$absdistance <= 1000,]
chr11_1k$group <- "group"
chr11_1k$group[chr11_1k$value >=-1 & chr11_1k$value < -0.9] <- "-1 to -0.9"
chr11_1k$group[chr11_1k$value >=-0.9 & chr11_1k$value < -0.8] <- "-0.9 to -0.8"
chr11_1k$group[chr11_1k$value >=-0.8 & chr11_1k$value < -0.7] <- "-0.8 to -0.7"
chr11_1k$group[chr11_1k$value >=-0.7 & chr11_1k$value < -0.6] <- "-0.7 to -0.6"
chr11_1k$group[chr11_1k$value >=-0.6 & chr11_1k$value < -0.5] <- "-0.6 to -0.5"
chr11_1k$group[chr11_1k$value >=-0.5 & chr11_1k$value < -0.4] <- "-0.5 to -0.4"
chr11_1k$group[chr11_1k$value >=-0.4 & chr11_1k$value < -0.3] <- "-0.4 to -0.3"
chr11_1k$group[chr11_1k$value >=-0.3 & chr11_1k$value < -0.2] <- "-0.3 to -0.2"
chr11_1k$group[chr11_1k$value >=-0.2 & chr11_1k$value < -0.1] <- "-0.2 to -0.1"
chr11_1k$group[chr11_1k$value >=-0.1 & chr11_1k$value < -0] <- "-0.1 to 0"
chr11_1k$group[chr11_1k$value >=0 & chr11_1k$value < 0.1] <- "0 to 0.1"
chr11_1k$group[chr11_1k$value >=0.1 & chr11_1k$value < 0.2] <- "0.1 to 0.2"
chr11_1k$group[chr11_1k$value >=0.2 & chr11_1k$value < 0.3] <- "0.2 to 0.3"
chr11_1k$group[chr11_1k$value >=0.3 & chr11_1k$value < 0.4] <- "0.3 to 0.4"
chr11_1k$group[chr11_1k$value >=0.4 & chr11_1k$value < 0.5] <- "0.4 to 0.5"
chr11_1k$group[chr11_1k$value >=0.5 & chr11_1k$value < 0.6] <- "0.5 to 0.6"
chr11_1k$group[chr11_1k$value >=0.6 & chr11_1k$value < 0.7] <- "0.6 to 0.7"
chr11_1k$group[chr11_1k$value >=0.7 & chr11_1k$value < 0.8] <- "0.7 to 0.8"
chr11_1k$group[chr11_1k$value >=0.8 & chr11_1k$value < 0.9] <- "0.8 to 0.9"
chr11_1k$group[chr11_1k$value >=0.9 & chr11_1k$value < 1] <- "0.9 to 1"

chr11_1k$group <- factor(chr11_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr11_1k$group)
cor.group <- as.data.frame(table(chr11_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr11_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 11: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()





###### chr12 #####

load("/path/to/cis_correlations/chr12_1.Rdata")
load("/path/to/cis_correlations/chr12_2.Rdata")
load("/path/to/cis_correlations/chr12_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr12_1_dim <- lapply(chr12, FUN = getdim)
chr12_1_dim_total <- lapply(chr12_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr12_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr12_2_dim <- lapply(chr12_2, FUN = getdim)
chr12_2_dim_total <- lapply(chr12_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr12_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr12_3_dim <- lapply(chr12_3, FUN = getdim)
chr12_3_dim_total <- lapply(chr12_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr12_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr12_1 <- rbindlist(chr12)
dim(chr12_1)
chr12_2 <- rbindlist(chr12_2)
dim(chr12_2)
chr12_3 <- rbindlist(chr12_3)
dim(chr12_3)

chr12 <- rbind(chr12_1, chr12_2)
dim(chr12)
chr12 <- rbind(chr12, chr12_3)
dim(chr12)
F7_cis_dim <- list()
F7_cis_dim$chr12_F7 <- dim(chr12)

chr12_10k <- chr12[chr12$absdistance <= 10000]
dim(chr12_10k)
head(chr12_10k)
chr12_10k_sort <- chr12_10k[order(absdistance)]
head(chr12_10k_sort)

# split to positive and negative correlations
chr12_10k_pos <- chr12_10k_sort[chr12_10k_sort$value >= 0,]
dim(chr12_10k_pos)
chr12_10k_neg <- chr12_10k_sort[chr12_10k_sort$value < 0,]
dim(chr12_10k_neg)

# bin the positive correlations on distance
chr12_10k_pos$group <- as.numeric(cut2(chr12_10k_pos$absdistance, m=100))
summary(chr12_10k_pos$group)
chr12_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr12_10k_pos_valuemean <- aggregate(chr12_10k_pos[, 3], list(chr12_10k_pos$group), mean)
class(chr12_10k_pos_valuemean)
head(chr12_10k_pos_valuemean)
# calculate the median distance of each bin
chr12_10k_pos_mediandistance <- aggregate(chr12_10k_pos[, 4], list(chr12_10k_pos$group), median)
class(chr12_10k_pos_mediandistance)
head(chr12_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr12_10k_pos_df <- merge(chr12_10k_pos_valuemean, chr12_10k_pos_mediandistance, by=c("Group.1"))
head(chr12_10k_pos_df)

# bin the negative correlations on distance
chr12_10k_neg$group <- as.numeric(cut2(chr12_10k_neg$absdistance, m=100))
summary(chr12_10k_neg$group)
chr12_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr12_10k_neg_valuemean <- aggregate(chr12_10k_neg[, 3], list(chr12_10k_neg$group), mean)
class(chr12_10k_neg_valuemean)
head(chr12_10k_neg_valuemean)
# calculate the median distance of each bin
chr12_10k_neg_mediandistance <- aggregate(chr12_10k_neg[, 4], list(chr12_10k_neg$group), median)
class(chr12_10k_neg_mediandistance)
head(chr12_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr12_10k_neg_df <- merge(chr12_10k_neg_valuemean, chr12_10k_neg_mediandistance, by=c("Group.1"))
head(chr12_10k_neg_df)

# calculate the SD in each positive bin
chr12_10k_pos_sd <- aggregate(chr12_10k_pos[, 3], list(chr12_10k_pos$group), sd)
class(chr12_10k_pos_sd)
head(chr12_10k_pos_sd)
chr12_10k_pos_df <- merge(chr12_10k_pos_df, chr12_10k_pos_sd, by=c("Group.1"))
head(chr12_10k_pos_df)
chr12_10k_pos_df$n_in_group <- as.numeric(table(chr12_10k_pos$group))
head(chr12_10k_pos_df)
table(chr12_10k_pos$group)

# calculate the SD in each negative bin
chr12_10k_neg_sd <- aggregate(chr12_10k_neg[, 3], list(chr12_10k_neg$group), sd)
class(chr12_10k_neg_sd)
head(chr12_10k_neg_sd)
chr12_10k_neg_df <- merge(chr12_10k_neg_df, chr12_10k_neg_sd, by=c("Group.1"))
head(chr12_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr12_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr12_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr12_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr12_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr12_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 12: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr12_1k <- chr12_10k_sort[chr12_10k_sort$absdistance <= 1000,]
chr12_1k$group <- "group"
chr12_1k$group[chr12_1k$value >=-1 & chr12_1k$value < -0.9] <- "-1 to -0.9"
chr12_1k$group[chr12_1k$value >=-0.9 & chr12_1k$value < -0.8] <- "-0.9 to -0.8"
chr12_1k$group[chr12_1k$value >=-0.8 & chr12_1k$value < -0.7] <- "-0.8 to -0.7"
chr12_1k$group[chr12_1k$value >=-0.7 & chr12_1k$value < -0.6] <- "-0.7 to -0.6"
chr12_1k$group[chr12_1k$value >=-0.6 & chr12_1k$value < -0.5] <- "-0.6 to -0.5"
chr12_1k$group[chr12_1k$value >=-0.5 & chr12_1k$value < -0.4] <- "-0.5 to -0.4"
chr12_1k$group[chr12_1k$value >=-0.4 & chr12_1k$value < -0.3] <- "-0.4 to -0.3"
chr12_1k$group[chr12_1k$value >=-0.3 & chr12_1k$value < -0.2] <- "-0.3 to -0.2"
chr12_1k$group[chr12_1k$value >=-0.2 & chr12_1k$value < -0.1] <- "-0.2 to -0.1"
chr12_1k$group[chr12_1k$value >=-0.1 & chr12_1k$value < -0] <- "-0.1 to 0"
chr12_1k$group[chr12_1k$value >=0 & chr12_1k$value < 0.1] <- "0 to 0.1"
chr12_1k$group[chr12_1k$value >=0.1 & chr12_1k$value < 0.2] <- "0.1 to 0.2"
chr12_1k$group[chr12_1k$value >=0.2 & chr12_1k$value < 0.3] <- "0.2 to 0.3"
chr12_1k$group[chr12_1k$value >=0.3 & chr12_1k$value < 0.4] <- "0.3 to 0.4"
chr12_1k$group[chr12_1k$value >=0.4 & chr12_1k$value < 0.5] <- "0.4 to 0.5"
chr12_1k$group[chr12_1k$value >=0.5 & chr12_1k$value < 0.6] <- "0.5 to 0.6"
chr12_1k$group[chr12_1k$value >=0.6 & chr12_1k$value < 0.7] <- "0.6 to 0.7"
chr12_1k$group[chr12_1k$value >=0.7 & chr12_1k$value < 0.8] <- "0.7 to 0.8"
chr12_1k$group[chr12_1k$value >=0.8 & chr12_1k$value < 0.9] <- "0.8 to 0.9"
chr12_1k$group[chr12_1k$value >=0.9 & chr12_1k$value < 1] <- "0.9 to 1"

chr12_1k$group <- factor(chr12_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr12_1k$group)
cor.group <- as.data.frame(table(chr12_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr12_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 12: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr13 #####

load("/path/to/cis_correlations/chr13_1.Rdata")
load("/path/to/cis_correlations/chr13_2.Rdata")
load("/path/to/cis_correlations/chr13_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr13_1_dim <- lapply(chr13, FUN = getdim)
chr13_1_dim_total <- lapply(chr13_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr13_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr13_2_dim <- lapply(chr13_2, FUN = getdim)
chr13_2_dim_total <- lapply(chr13_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr13_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr13_3_dim <- lapply(chr13_3, FUN = getdim)
chr13_3_dim_total <- lapply(chr13_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr13_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr13_1 <- rbindlist(chr13)
dim(chr13_1)
chr13_2 <- rbindlist(chr13_2)
dim(chr13_2)
chr13_3 <- rbindlist(chr13_3)
dim(chr13_3)

chr13 <- rbind(chr13_1, chr13_2)
dim(chr13)
chr13 <- rbind(chr13, chr13_3)
dim(chr13)
F7_cis_dim <- list()
F7_cis_dim$chr13_F7 <- dim(chr13)

chr13_10k <- chr13[chr13$absdistance <= 10000]
dim(chr13_10k)
head(chr13_10k)
chr13_10k_sort <- chr13_10k[order(absdistance)]
head(chr13_10k_sort)

# split to positive and negative correlations
chr13_10k_pos <- chr13_10k_sort[chr13_10k_sort$value >= 0,]
dim(chr13_10k_pos)
chr13_10k_neg <- chr13_10k_sort[chr13_10k_sort$value < 0,]
dim(chr13_10k_neg)

# bin the positive correlations on distance
chr13_10k_pos$group <- as.numeric(cut2(chr13_10k_pos$absdistance, m=100))
summary(chr13_10k_pos$group)
chr13_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr13_10k_pos_valuemean <- aggregate(chr13_10k_pos[, 3], list(chr13_10k_pos$group), mean)
class(chr13_10k_pos_valuemean)
head(chr13_10k_pos_valuemean)
# calculate the median distance of each bin
chr13_10k_pos_mediandistance <- aggregate(chr13_10k_pos[, 4], list(chr13_10k_pos$group), median)
class(chr13_10k_pos_mediandistance)
head(chr13_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr13_10k_pos_df <- merge(chr13_10k_pos_valuemean, chr13_10k_pos_mediandistance, by=c("Group.1"))
head(chr13_10k_pos_df)

# bin the negative correlations on distance
chr13_10k_neg$group <- as.numeric(cut2(chr13_10k_neg$absdistance, m=100))
summary(chr13_10k_neg$group)
chr13_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr13_10k_neg_valuemean <- aggregate(chr13_10k_neg[, 3], list(chr13_10k_neg$group), mean)
class(chr13_10k_neg_valuemean)
head(chr13_10k_neg_valuemean)
# calculate the median distance of each bin
chr13_10k_neg_mediandistance <- aggregate(chr13_10k_neg[, 4], list(chr13_10k_neg$group), median)
class(chr13_10k_neg_mediandistance)
head(chr13_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr13_10k_neg_df <- merge(chr13_10k_neg_valuemean, chr13_10k_neg_mediandistance, by=c("Group.1"))
head(chr13_10k_neg_df)

# calculate the SD in each positive bin
chr13_10k_pos_sd <- aggregate(chr13_10k_pos[, 3], list(chr13_10k_pos$group), sd)
class(chr13_10k_pos_sd)
head(chr13_10k_pos_sd)
chr13_10k_pos_df <- merge(chr13_10k_pos_df, chr13_10k_pos_sd, by=c("Group.1"))
head(chr13_10k_pos_df)
chr13_10k_pos_df$n_in_group <- as.numeric(table(chr13_10k_pos$group))
head(chr13_10k_pos_df)
table(chr13_10k_pos$group)

# calculate the SD in each negative bin
chr13_10k_neg_sd <- aggregate(chr13_10k_neg[, 3], list(chr13_10k_neg$group), sd)
class(chr13_10k_neg_sd)
head(chr13_10k_neg_sd)
chr13_10k_neg_df <- merge(chr13_10k_neg_df, chr13_10k_neg_sd, by=c("Group.1"))
head(chr13_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr13_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr13_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr13_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr13_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr13_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 13: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr13_1k <- chr13_10k_sort[chr13_10k_sort$absdistance <= 1000,]
chr13_1k$group <- "group"
chr13_1k$group[chr13_1k$value >=-1 & chr13_1k$value < -0.9] <- "-1 to -0.9"
chr13_1k$group[chr13_1k$value >=-0.9 & chr13_1k$value < -0.8] <- "-0.9 to -0.8"
chr13_1k$group[chr13_1k$value >=-0.8 & chr13_1k$value < -0.7] <- "-0.8 to -0.7"
chr13_1k$group[chr13_1k$value >=-0.7 & chr13_1k$value < -0.6] <- "-0.7 to -0.6"
chr13_1k$group[chr13_1k$value >=-0.6 & chr13_1k$value < -0.5] <- "-0.6 to -0.5"
chr13_1k$group[chr13_1k$value >=-0.5 & chr13_1k$value < -0.4] <- "-0.5 to -0.4"
chr13_1k$group[chr13_1k$value >=-0.4 & chr13_1k$value < -0.3] <- "-0.4 to -0.3"
chr13_1k$group[chr13_1k$value >=-0.3 & chr13_1k$value < -0.2] <- "-0.3 to -0.2"
chr13_1k$group[chr13_1k$value >=-0.2 & chr13_1k$value < -0.1] <- "-0.2 to -0.1"
chr13_1k$group[chr13_1k$value >=-0.1 & chr13_1k$value < -0] <- "-0.1 to 0"
chr13_1k$group[chr13_1k$value >=0 & chr13_1k$value < 0.1] <- "0 to 0.1"
chr13_1k$group[chr13_1k$value >=0.1 & chr13_1k$value < 0.2] <- "0.1 to 0.2"
chr13_1k$group[chr13_1k$value >=0.2 & chr13_1k$value < 0.3] <- "0.2 to 0.3"
chr13_1k$group[chr13_1k$value >=0.3 & chr13_1k$value < 0.4] <- "0.3 to 0.4"
chr13_1k$group[chr13_1k$value >=0.4 & chr13_1k$value < 0.5] <- "0.4 to 0.5"
chr13_1k$group[chr13_1k$value >=0.5 & chr13_1k$value < 0.6] <- "0.5 to 0.6"
chr13_1k$group[chr13_1k$value >=0.6 & chr13_1k$value < 0.7] <- "0.6 to 0.7"
chr13_1k$group[chr13_1k$value >=0.7 & chr13_1k$value < 0.8] <- "0.7 to 0.8"
chr13_1k$group[chr13_1k$value >=0.8 & chr13_1k$value < 0.9] <- "0.8 to 0.9"
chr13_1k$group[chr13_1k$value >=0.9 & chr13_1k$value < 1] <- "0.9 to 1"

chr13_1k$group <- factor(chr13_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr13_1k$group)
cor.group <- as.data.frame(table(chr13_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr13_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 13: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()



###### chr14 #####

load("/path/to/cis_correlations/chr14_1.Rdata")
load("/path/to/cis_correlations/chr14_2.Rdata")
load("/path/to/cis_correlations/chr14_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr14_1_dim <- lapply(chr14, FUN = getdim)
chr14_1_dim_total <- lapply(chr14_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr14_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr14_2_dim <- lapply(chr14_2, FUN = getdim)
chr14_2_dim_total <- lapply(chr14_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr14_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr14_3_dim <- lapply(chr14_3, FUN = getdim)
chr14_3_dim_total <- lapply(chr14_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr14_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr14_1 <- rbindlist(chr14)
dim(chr14_1)
chr14_2 <- rbindlist(chr14_2)
dim(chr14_2)
chr14_3 <- rbindlist(chr14_3)
dim(chr14_3)

chr14 <- rbind(chr14_1, chr14_2)
dim(chr14)
chr14 <- rbind(chr14, chr14_3)
dim(chr14)
F7_cis_dim <- list()
F7_cis_dim$chr14_F7 <- dim(chr14)

chr14_10k <- chr14[chr14$absdistance <= 10000]
dim(chr14_10k)
head(chr14_10k)
chr14_10k_sort <- chr14_10k[order(absdistance)]
head(chr14_10k_sort)

# split to positive and negative correlations
chr14_10k_pos <- chr14_10k_sort[chr14_10k_sort$value >= 0,]
dim(chr14_10k_pos)
chr14_10k_neg <- chr14_10k_sort[chr14_10k_sort$value < 0,]
dim(chr14_10k_neg)

# bin the positive correlations on distance
chr14_10k_pos$group <- as.numeric(cut2(chr14_10k_pos$absdistance, m=100))
summary(chr14_10k_pos$group)
chr14_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr14_10k_pos_valuemean <- aggregate(chr14_10k_pos[, 3], list(chr14_10k_pos$group), mean)
class(chr14_10k_pos_valuemean)
head(chr14_10k_pos_valuemean)
# calculate the median distance of each bin
chr14_10k_pos_mediandistance <- aggregate(chr14_10k_pos[, 4], list(chr14_10k_pos$group), median)
class(chr14_10k_pos_mediandistance)
head(chr14_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr14_10k_pos_df <- merge(chr14_10k_pos_valuemean, chr14_10k_pos_mediandistance, by=c("Group.1"))
head(chr14_10k_pos_df)

# bin the negative correlations on distance
chr14_10k_neg$group <- as.numeric(cut2(chr14_10k_neg$absdistance, m=100))
summary(chr14_10k_neg$group)
chr14_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr14_10k_neg_valuemean <- aggregate(chr14_10k_neg[, 3], list(chr14_10k_neg$group), mean)
class(chr14_10k_neg_valuemean)
head(chr14_10k_neg_valuemean)
# calculate the median distance of each bin
chr14_10k_neg_mediandistance <- aggregate(chr14_10k_neg[, 4], list(chr14_10k_neg$group), median)
class(chr14_10k_neg_mediandistance)
head(chr14_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr14_10k_neg_df <- merge(chr14_10k_neg_valuemean, chr14_10k_neg_mediandistance, by=c("Group.1"))
head(chr14_10k_neg_df)

# calculate the SD in each positive bin
chr14_10k_pos_sd <- aggregate(chr14_10k_pos[, 3], list(chr14_10k_pos$group), sd)
class(chr14_10k_pos_sd)
head(chr14_10k_pos_sd)
chr14_10k_pos_df <- merge(chr14_10k_pos_df, chr14_10k_pos_sd, by=c("Group.1"))
head(chr14_10k_pos_df)
chr14_10k_pos_df$n_in_group <- as.numeric(table(chr14_10k_pos$group))
head(chr14_10k_pos_df)
table(chr14_10k_pos$group)

# calculate the SD in each negative bin
chr14_10k_neg_sd <- aggregate(chr14_10k_neg[, 3], list(chr14_10k_neg$group), sd)
class(chr14_10k_neg_sd)
head(chr14_10k_neg_sd)
chr14_10k_neg_df <- merge(chr14_10k_neg_df, chr14_10k_neg_sd, by=c("Group.1"))
head(chr14_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr14_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr14_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr14_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr14_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr14_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 14: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr14_1k <- chr14_10k_sort[chr14_10k_sort$absdistance <= 1000,]
chr14_1k$group <- "group"
chr14_1k$group[chr14_1k$value >=-1 & chr14_1k$value < -0.9] <- "-1 to -0.9"
chr14_1k$group[chr14_1k$value >=-0.9 & chr14_1k$value < -0.8] <- "-0.9 to -0.8"
chr14_1k$group[chr14_1k$value >=-0.8 & chr14_1k$value < -0.7] <- "-0.8 to -0.7"
chr14_1k$group[chr14_1k$value >=-0.7 & chr14_1k$value < -0.6] <- "-0.7 to -0.6"
chr14_1k$group[chr14_1k$value >=-0.6 & chr14_1k$value < -0.5] <- "-0.6 to -0.5"
chr14_1k$group[chr14_1k$value >=-0.5 & chr14_1k$value < -0.4] <- "-0.5 to -0.4"
chr14_1k$group[chr14_1k$value >=-0.4 & chr14_1k$value < -0.3] <- "-0.4 to -0.3"
chr14_1k$group[chr14_1k$value >=-0.3 & chr14_1k$value < -0.2] <- "-0.3 to -0.2"
chr14_1k$group[chr14_1k$value >=-0.2 & chr14_1k$value < -0.1] <- "-0.2 to -0.1"
chr14_1k$group[chr14_1k$value >=-0.1 & chr14_1k$value < -0] <- "-0.1 to 0"
chr14_1k$group[chr14_1k$value >=0 & chr14_1k$value < 0.1] <- "0 to 0.1"
chr14_1k$group[chr14_1k$value >=0.1 & chr14_1k$value < 0.2] <- "0.1 to 0.2"
chr14_1k$group[chr14_1k$value >=0.2 & chr14_1k$value < 0.3] <- "0.2 to 0.3"
chr14_1k$group[chr14_1k$value >=0.3 & chr14_1k$value < 0.4] <- "0.3 to 0.4"
chr14_1k$group[chr14_1k$value >=0.4 & chr14_1k$value < 0.5] <- "0.4 to 0.5"
chr14_1k$group[chr14_1k$value >=0.5 & chr14_1k$value < 0.6] <- "0.5 to 0.6"
chr14_1k$group[chr14_1k$value >=0.6 & chr14_1k$value < 0.7] <- "0.6 to 0.7"
chr14_1k$group[chr14_1k$value >=0.7 & chr14_1k$value < 0.8] <- "0.7 to 0.8"
chr14_1k$group[chr14_1k$value >=0.8 & chr14_1k$value < 0.9] <- "0.8 to 0.9"
chr14_1k$group[chr14_1k$value >=0.9 & chr14_1k$value < 1] <- "0.9 to 1"

chr14_1k$group <- factor(chr14_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr14_1k$group)
cor.group <- as.data.frame(table(chr14_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr14_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 14: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr15 #####

load("/path/to/cis_correlations/chr15_1.Rdata")
load("/path/to/cis_correlations/chr15_2.Rdata")
load("/path/to/cis_correlations/chr15_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr15_1_dim <- lapply(chr15, FUN = getdim)
chr15_1_dim_total <- lapply(chr15_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr15_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr15_2_dim <- lapply(chr15_2, FUN = getdim)
chr15_2_dim_total <- lapply(chr15_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr15_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr15_3_dim <- lapply(chr15_3, FUN = getdim)
chr15_3_dim_total <- lapply(chr15_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr15_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr15_1 <- rbindlist(chr15)
dim(chr15_1)
chr15_2 <- rbindlist(chr15_2)
dim(chr15_2)
chr15_3 <- rbindlist(chr15_3)
dim(chr15_3)

chr15 <- rbind(chr15_1, chr15_2)
dim(chr15)
chr15 <- rbind(chr15, chr15_3)
dim(chr15)
F7_cis_dim <- list()
F7_cis_dim$chr15_F7 <- dim(chr15)

chr15_10k <- chr15[chr15$absdistance <= 10000]
dim(chr15_10k)
head(chr15_10k)
chr15_10k_sort <- chr15_10k[order(absdistance)]
head(chr15_10k_sort)

# split to positive and negative correlations
chr15_10k_pos <- chr15_10k_sort[chr15_10k_sort$value >= 0,]
dim(chr15_10k_pos)
chr15_10k_neg <- chr15_10k_sort[chr15_10k_sort$value < 0,]
dim(chr15_10k_neg)

# bin the positive correlations on distance
chr15_10k_pos$group <- as.numeric(cut2(chr15_10k_pos$absdistance, m=100))
summary(chr15_10k_pos$group)
chr15_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr15_10k_pos_valuemean <- aggregate(chr15_10k_pos[, 3], list(chr15_10k_pos$group), mean)
class(chr15_10k_pos_valuemean)
head(chr15_10k_pos_valuemean)
# calculate the median distance of each bin
chr15_10k_pos_mediandistance <- aggregate(chr15_10k_pos[, 4], list(chr15_10k_pos$group), median)
class(chr15_10k_pos_mediandistance)
head(chr15_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr15_10k_pos_df <- merge(chr15_10k_pos_valuemean, chr15_10k_pos_mediandistance, by=c("Group.1"))
head(chr15_10k_pos_df)

# bin the negative correlations on distance
chr15_10k_neg$group <- as.numeric(cut2(chr15_10k_neg$absdistance, m=100))
summary(chr15_10k_neg$group)
chr15_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr15_10k_neg_valuemean <- aggregate(chr15_10k_neg[, 3], list(chr15_10k_neg$group), mean)
class(chr15_10k_neg_valuemean)
head(chr15_10k_neg_valuemean)
# calculate the median distance of each bin
chr15_10k_neg_mediandistance <- aggregate(chr15_10k_neg[, 4], list(chr15_10k_neg$group), median)
class(chr15_10k_neg_mediandistance)
head(chr15_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr15_10k_neg_df <- merge(chr15_10k_neg_valuemean, chr15_10k_neg_mediandistance, by=c("Group.1"))
head(chr15_10k_neg_df)

# calculate the SD in each positive bin
chr15_10k_pos_sd <- aggregate(chr15_10k_pos[, 3], list(chr15_10k_pos$group), sd)
class(chr15_10k_pos_sd)
head(chr15_10k_pos_sd)
chr15_10k_pos_df <- merge(chr15_10k_pos_df, chr15_10k_pos_sd, by=c("Group.1"))
head(chr15_10k_pos_df)
chr15_10k_pos_df$n_in_group <- as.numeric(table(chr15_10k_pos$group))
head(chr15_10k_pos_df)
table(chr15_10k_pos$group)

# calculate the SD in each negative bin
chr15_10k_neg_sd <- aggregate(chr15_10k_neg[, 3], list(chr15_10k_neg$group), sd)
class(chr15_10k_neg_sd)
head(chr15_10k_neg_sd)
chr15_10k_neg_df <- merge(chr15_10k_neg_df, chr15_10k_neg_sd, by=c("Group.1"))
head(chr15_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr15_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr15_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr15_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr15_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr15_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 15: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr15_1k <- chr15_10k_sort[chr15_10k_sort$absdistance <= 1000,]
chr15_1k$group <- "group"
chr15_1k$group[chr15_1k$value >=-1 & chr15_1k$value < -0.9] <- "-1 to -0.9"
chr15_1k$group[chr15_1k$value >=-0.9 & chr15_1k$value < -0.8] <- "-0.9 to -0.8"
chr15_1k$group[chr15_1k$value >=-0.8 & chr15_1k$value < -0.7] <- "-0.8 to -0.7"
chr15_1k$group[chr15_1k$value >=-0.7 & chr15_1k$value < -0.6] <- "-0.7 to -0.6"
chr15_1k$group[chr15_1k$value >=-0.6 & chr15_1k$value < -0.5] <- "-0.6 to -0.5"
chr15_1k$group[chr15_1k$value >=-0.5 & chr15_1k$value < -0.4] <- "-0.5 to -0.4"
chr15_1k$group[chr15_1k$value >=-0.4 & chr15_1k$value < -0.3] <- "-0.4 to -0.3"
chr15_1k$group[chr15_1k$value >=-0.3 & chr15_1k$value < -0.2] <- "-0.3 to -0.2"
chr15_1k$group[chr15_1k$value >=-0.2 & chr15_1k$value < -0.1] <- "-0.2 to -0.1"
chr15_1k$group[chr15_1k$value >=-0.1 & chr15_1k$value < -0] <- "-0.1 to 0"
chr15_1k$group[chr15_1k$value >=0 & chr15_1k$value < 0.1] <- "0 to 0.1"
chr15_1k$group[chr15_1k$value >=0.1 & chr15_1k$value < 0.2] <- "0.1 to 0.2"
chr15_1k$group[chr15_1k$value >=0.2 & chr15_1k$value < 0.3] <- "0.2 to 0.3"
chr15_1k$group[chr15_1k$value >=0.3 & chr15_1k$value < 0.4] <- "0.3 to 0.4"
chr15_1k$group[chr15_1k$value >=0.4 & chr15_1k$value < 0.5] <- "0.4 to 0.5"
chr15_1k$group[chr15_1k$value >=0.5 & chr15_1k$value < 0.6] <- "0.5 to 0.6"
chr15_1k$group[chr15_1k$value >=0.6 & chr15_1k$value < 0.7] <- "0.6 to 0.7"
chr15_1k$group[chr15_1k$value >=0.7 & chr15_1k$value < 0.8] <- "0.7 to 0.8"
chr15_1k$group[chr15_1k$value >=0.8 & chr15_1k$value < 0.9] <- "0.8 to 0.9"
chr15_1k$group[chr15_1k$value >=0.9 & chr15_1k$value < 1] <- "0.9 to 1"

chr15_1k$group <- factor(chr15_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr15_1k$group)
cor.group <- as.data.frame(table(chr15_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr15_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 15: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr16 #####

load("/path/to/cis_correlations/chr16_1.Rdata")
load("/path/to/cis_correlations/chr16_2.Rdata")
load("/path/to/cis_correlations/chr16_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr16_1_dim <- lapply(chr16, FUN = getdim)
chr16_1_dim_total <- lapply(chr16_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr16_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr16_2_dim <- lapply(chr16_2, FUN = getdim)
chr16_2_dim_total <- lapply(chr16_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr16_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr16_3_dim <- lapply(chr16_3, FUN = getdim)
chr16_3_dim_total <- lapply(chr16_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr16_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr16_1 <- rbindlist(chr16)
dim(chr16_1)
chr16_2 <- rbindlist(chr16_2)
dim(chr16_2)
chr16_3 <- rbindlist(chr16_3)
dim(chr16_3)

chr16 <- rbind(chr16_1, chr16_2)
dim(chr16)
chr16 <- rbind(chr16, chr16_3)
dim(chr16)
F7_cis_dim <- list()
F7_cis_dim$chr16_F7 <- dim(chr16)

chr16_10k <- chr16[chr16$absdistance <= 10000]
dim(chr16_10k)
head(chr16_10k)
chr16_10k_sort <- chr16_10k[order(absdistance)]
head(chr16_10k_sort)

# split to positive and negative correlations
chr16_10k_pos <- chr16_10k_sort[chr16_10k_sort$value >= 0,]
dim(chr16_10k_pos)
chr16_10k_neg <- chr16_10k_sort[chr16_10k_sort$value < 0,]
dim(chr16_10k_neg)

# bin the positive correlations on distance
chr16_10k_pos$group <- as.numeric(cut2(chr16_10k_pos$absdistance, m=100))
summary(chr16_10k_pos$group)
chr16_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr16_10k_pos_valuemean <- aggregate(chr16_10k_pos[, 3], list(chr16_10k_pos$group), mean)
class(chr16_10k_pos_valuemean)
head(chr16_10k_pos_valuemean)
# calculate the median distance of each bin
chr16_10k_pos_mediandistance <- aggregate(chr16_10k_pos[, 4], list(chr16_10k_pos$group), median)
class(chr16_10k_pos_mediandistance)
head(chr16_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr16_10k_pos_df <- merge(chr16_10k_pos_valuemean, chr16_10k_pos_mediandistance, by=c("Group.1"))
head(chr16_10k_pos_df)

# bin the negative correlations on distance
chr16_10k_neg$group <- as.numeric(cut2(chr16_10k_neg$absdistance, m=100))
summary(chr16_10k_neg$group)
chr16_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr16_10k_neg_valuemean <- aggregate(chr16_10k_neg[, 3], list(chr16_10k_neg$group), mean)
class(chr16_10k_neg_valuemean)
head(chr16_10k_neg_valuemean)
# calculate the median distance of each bin
chr16_10k_neg_mediandistance <- aggregate(chr16_10k_neg[, 4], list(chr16_10k_neg$group), median)
class(chr16_10k_neg_mediandistance)
head(chr16_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr16_10k_neg_df <- merge(chr16_10k_neg_valuemean, chr16_10k_neg_mediandistance, by=c("Group.1"))
head(chr16_10k_neg_df)

# calculate the SD in each positive bin
chr16_10k_pos_sd <- aggregate(chr16_10k_pos[, 3], list(chr16_10k_pos$group), sd)
class(chr16_10k_pos_sd)
head(chr16_10k_pos_sd)
chr16_10k_pos_df <- merge(chr16_10k_pos_df, chr16_10k_pos_sd, by=c("Group.1"))
head(chr16_10k_pos_df)
chr16_10k_pos_df$n_in_group <- as.numeric(table(chr16_10k_pos$group))
head(chr16_10k_pos_df)
table(chr16_10k_pos$group)

# calculate the SD in each negative bin
chr16_10k_neg_sd <- aggregate(chr16_10k_neg[, 3], list(chr16_10k_neg$group), sd)
class(chr16_10k_neg_sd)
head(chr16_10k_neg_sd)
chr16_10k_neg_df <- merge(chr16_10k_neg_df, chr16_10k_neg_sd, by=c("Group.1"))
head(chr16_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr16_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr16_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr16_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr16_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr16_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 16: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr16_1k <- chr16_10k_sort[chr16_10k_sort$absdistance <= 1000,]
chr16_1k$group <- "group"
chr16_1k$group[chr16_1k$value >=-1 & chr16_1k$value < -0.9] <- "-1 to -0.9"
chr16_1k$group[chr16_1k$value >=-0.9 & chr16_1k$value < -0.8] <- "-0.9 to -0.8"
chr16_1k$group[chr16_1k$value >=-0.8 & chr16_1k$value < -0.7] <- "-0.8 to -0.7"
chr16_1k$group[chr16_1k$value >=-0.7 & chr16_1k$value < -0.6] <- "-0.7 to -0.6"
chr16_1k$group[chr16_1k$value >=-0.6 & chr16_1k$value < -0.5] <- "-0.6 to -0.5"
chr16_1k$group[chr16_1k$value >=-0.5 & chr16_1k$value < -0.4] <- "-0.5 to -0.4"
chr16_1k$group[chr16_1k$value >=-0.4 & chr16_1k$value < -0.3] <- "-0.4 to -0.3"
chr16_1k$group[chr16_1k$value >=-0.3 & chr16_1k$value < -0.2] <- "-0.3 to -0.2"
chr16_1k$group[chr16_1k$value >=-0.2 & chr16_1k$value < -0.1] <- "-0.2 to -0.1"
chr16_1k$group[chr16_1k$value >=-0.1 & chr16_1k$value < -0] <- "-0.1 to 0"
chr16_1k$group[chr16_1k$value >=0 & chr16_1k$value < 0.1] <- "0 to 0.1"
chr16_1k$group[chr16_1k$value >=0.1 & chr16_1k$value < 0.2] <- "0.1 to 0.2"
chr16_1k$group[chr16_1k$value >=0.2 & chr16_1k$value < 0.3] <- "0.2 to 0.3"
chr16_1k$group[chr16_1k$value >=0.3 & chr16_1k$value < 0.4] <- "0.3 to 0.4"
chr16_1k$group[chr16_1k$value >=0.4 & chr16_1k$value < 0.5] <- "0.4 to 0.5"
chr16_1k$group[chr16_1k$value >=0.5 & chr16_1k$value < 0.6] <- "0.5 to 0.6"
chr16_1k$group[chr16_1k$value >=0.6 & chr16_1k$value < 0.7] <- "0.6 to 0.7"
chr16_1k$group[chr16_1k$value >=0.7 & chr16_1k$value < 0.8] <- "0.7 to 0.8"
chr16_1k$group[chr16_1k$value >=0.8 & chr16_1k$value < 0.9] <- "0.8 to 0.9"
chr16_1k$group[chr16_1k$value >=0.9 & chr16_1k$value < 1] <- "0.9 to 1"

chr16_1k$group <- factor(chr16_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr16_1k$group)
cor.group <- as.data.frame(table(chr16_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr16_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 16: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr17 #####

load("/path/to/cis_correlations/chr17_1.Rdata")
load("/path/to/cis_correlations/chr17_2.Rdata")
load("/path/to/cis_correlations/chr17_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr17_1_dim <- lapply(chr17, FUN = getdim)
chr17_1_dim_total <- lapply(chr17_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr17_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr17_2_dim <- lapply(chr17_2, FUN = getdim)
chr17_2_dim_total <- lapply(chr17_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr17_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr17_3_dim <- lapply(chr17_3, FUN = getdim)
chr17_3_dim_total <- lapply(chr17_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr17_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr17_1 <- rbindlist(chr17)
dim(chr17_1)
chr17_2 <- rbindlist(chr17_2)
dim(chr17_2)
chr17_3 <- rbindlist(chr17_3)
dim(chr17_3)

chr17 <- rbind(chr17_1, chr17_2)
dim(chr17)
chr17 <- rbind(chr17, chr17_3)
dim(chr17)
F7_cis_dim <- list()
F7_cis_dim$chr17_F7 <- dim(chr17)

chr17_10k <- chr17[chr17$absdistance <= 10000]
dim(chr17_10k)
head(chr17_10k)
chr17_10k_sort <- chr17_10k[order(absdistance)]
head(chr17_10k_sort)

# split to positive and negative correlations
chr17_10k_pos <- chr17_10k_sort[chr17_10k_sort$value >= 0,]
dim(chr17_10k_pos)
chr17_10k_neg <- chr17_10k_sort[chr17_10k_sort$value < 0,]
dim(chr17_10k_neg)

# bin the positive correlations on distance
chr17_10k_pos$group <- as.numeric(cut2(chr17_10k_pos$absdistance, m=100))
summary(chr17_10k_pos$group)
chr17_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr17_10k_pos_valuemean <- aggregate(chr17_10k_pos[, 3], list(chr17_10k_pos$group), mean)
class(chr17_10k_pos_valuemean)
head(chr17_10k_pos_valuemean)
# calculate the median distance of each bin
chr17_10k_pos_mediandistance <- aggregate(chr17_10k_pos[, 4], list(chr17_10k_pos$group), median)
class(chr17_10k_pos_mediandistance)
head(chr17_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr17_10k_pos_df <- merge(chr17_10k_pos_valuemean, chr17_10k_pos_mediandistance, by=c("Group.1"))
head(chr17_10k_pos_df)

# bin the negative correlations on distance
chr17_10k_neg$group <- as.numeric(cut2(chr17_10k_neg$absdistance, m=100))
summary(chr17_10k_neg$group)
chr17_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr17_10k_neg_valuemean <- aggregate(chr17_10k_neg[, 3], list(chr17_10k_neg$group), mean)
class(chr17_10k_neg_valuemean)
head(chr17_10k_neg_valuemean)
# calculate the median distance of each bin
chr17_10k_neg_mediandistance <- aggregate(chr17_10k_neg[, 4], list(chr17_10k_neg$group), median)
class(chr17_10k_neg_mediandistance)
head(chr17_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr17_10k_neg_df <- merge(chr17_10k_neg_valuemean, chr17_10k_neg_mediandistance, by=c("Group.1"))
head(chr17_10k_neg_df)

# calculate the SD in each positive bin
chr17_10k_pos_sd <- aggregate(chr17_10k_pos[, 3], list(chr17_10k_pos$group), sd)
class(chr17_10k_pos_sd)
head(chr17_10k_pos_sd)
chr17_10k_pos_df <- merge(chr17_10k_pos_df, chr17_10k_pos_sd, by=c("Group.1"))
head(chr17_10k_pos_df)
chr17_10k_pos_df$n_in_group <- as.numeric(table(chr17_10k_pos$group))
head(chr17_10k_pos_df)
table(chr17_10k_pos$group)

# calculate the SD in each negative bin
chr17_10k_neg_sd <- aggregate(chr17_10k_neg[, 3], list(chr17_10k_neg$group), sd)
class(chr17_10k_neg_sd)
head(chr17_10k_neg_sd)
chr17_10k_neg_df <- merge(chr17_10k_neg_df, chr17_10k_neg_sd, by=c("Group.1"))
head(chr17_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr17_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr17_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr17_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr17_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr17_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 17: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr17_1k <- chr17_10k_sort[chr17_10k_sort$absdistance <= 1000,]
chr17_1k$group <- "group"
chr17_1k$group[chr17_1k$value >=-1 & chr17_1k$value < -0.9] <- "-1 to -0.9"
chr17_1k$group[chr17_1k$value >=-0.9 & chr17_1k$value < -0.8] <- "-0.9 to -0.8"
chr17_1k$group[chr17_1k$value >=-0.8 & chr17_1k$value < -0.7] <- "-0.8 to -0.7"
chr17_1k$group[chr17_1k$value >=-0.7 & chr17_1k$value < -0.6] <- "-0.7 to -0.6"
chr17_1k$group[chr17_1k$value >=-0.6 & chr17_1k$value < -0.5] <- "-0.6 to -0.5"
chr17_1k$group[chr17_1k$value >=-0.5 & chr17_1k$value < -0.4] <- "-0.5 to -0.4"
chr17_1k$group[chr17_1k$value >=-0.4 & chr17_1k$value < -0.3] <- "-0.4 to -0.3"
chr17_1k$group[chr17_1k$value >=-0.3 & chr17_1k$value < -0.2] <- "-0.3 to -0.2"
chr17_1k$group[chr17_1k$value >=-0.2 & chr17_1k$value < -0.1] <- "-0.2 to -0.1"
chr17_1k$group[chr17_1k$value >=-0.1 & chr17_1k$value < -0] <- "-0.1 to 0"
chr17_1k$group[chr17_1k$value >=0 & chr17_1k$value < 0.1] <- "0 to 0.1"
chr17_1k$group[chr17_1k$value >=0.1 & chr17_1k$value < 0.2] <- "0.1 to 0.2"
chr17_1k$group[chr17_1k$value >=0.2 & chr17_1k$value < 0.3] <- "0.2 to 0.3"
chr17_1k$group[chr17_1k$value >=0.3 & chr17_1k$value < 0.4] <- "0.3 to 0.4"
chr17_1k$group[chr17_1k$value >=0.4 & chr17_1k$value < 0.5] <- "0.4 to 0.5"
chr17_1k$group[chr17_1k$value >=0.5 & chr17_1k$value < 0.6] <- "0.5 to 0.6"
chr17_1k$group[chr17_1k$value >=0.6 & chr17_1k$value < 0.7] <- "0.6 to 0.7"
chr17_1k$group[chr17_1k$value >=0.7 & chr17_1k$value < 0.8] <- "0.7 to 0.8"
chr17_1k$group[chr17_1k$value >=0.8 & chr17_1k$value < 0.9] <- "0.8 to 0.9"
chr17_1k$group[chr17_1k$value >=0.9 & chr17_1k$value < 1] <- "0.9 to 1"

chr17_1k$group <- factor(chr17_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr17_1k$group)
cor.group <- as.data.frame(table(chr17_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr17_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 17: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr18 #####

load("/path/to/cis_correlations/chr18_1.Rdata")
load("/path/to/cis_correlations/chr18_2.Rdata")
load("/path/to/cis_correlations/chr18_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr18_1_dim <- lapply(chr18, FUN = getdim)
chr18_1_dim_total <- lapply(chr18_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr18_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr18_2_dim <- lapply(chr18_2, FUN = getdim)
chr18_2_dim_total <- lapply(chr18_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr18_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr18_3_dim <- lapply(chr18_3, FUN = getdim)
chr18_3_dim_total <- lapply(chr18_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr18_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr18_1 <- rbindlist(chr18)
dim(chr18_1)
chr18_2 <- rbindlist(chr18_2)
dim(chr18_2)
chr18_3 <- rbindlist(chr18_3)
dim(chr18_3)

chr18 <- rbind(chr18_1, chr18_2)
dim(chr18)
chr18 <- rbind(chr18, chr18_3)
dim(chr18)
F7_cis_dim <- list()
F7_cis_dim$chr18_F7 <- dim(chr18)

chr18_10k <- chr18[chr18$absdistance <= 10000]
dim(chr18_10k)
head(chr18_10k)
chr18_10k_sort <- chr18_10k[order(absdistance)]
head(chr18_10k_sort)

# split to positive and negative correlations
chr18_10k_pos <- chr18_10k_sort[chr18_10k_sort$value >= 0,]
dim(chr18_10k_pos)
chr18_10k_neg <- chr18_10k_sort[chr18_10k_sort$value < 0,]
dim(chr18_10k_neg)

# bin the positive correlations on distance
chr18_10k_pos$group <- as.numeric(cut2(chr18_10k_pos$absdistance, m=100))
summary(chr18_10k_pos$group)
chr18_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr18_10k_pos_valuemean <- aggregate(chr18_10k_pos[, 3], list(chr18_10k_pos$group), mean)
class(chr18_10k_pos_valuemean)
head(chr18_10k_pos_valuemean)
# calculate the median distance of each bin
chr18_10k_pos_mediandistance <- aggregate(chr18_10k_pos[, 4], list(chr18_10k_pos$group), median)
class(chr18_10k_pos_mediandistance)
head(chr18_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr18_10k_pos_df <- merge(chr18_10k_pos_valuemean, chr18_10k_pos_mediandistance, by=c("Group.1"))
head(chr18_10k_pos_df)

# bin the negative correlations on distance
chr18_10k_neg$group <- as.numeric(cut2(chr18_10k_neg$absdistance, m=100))
summary(chr18_10k_neg$group)
chr18_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr18_10k_neg_valuemean <- aggregate(chr18_10k_neg[, 3], list(chr18_10k_neg$group), mean)
class(chr18_10k_neg_valuemean)
head(chr18_10k_neg_valuemean)
# calculate the median distance of each bin
chr18_10k_neg_mediandistance <- aggregate(chr18_10k_neg[, 4], list(chr18_10k_neg$group), median)
class(chr18_10k_neg_mediandistance)
head(chr18_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr18_10k_neg_df <- merge(chr18_10k_neg_valuemean, chr18_10k_neg_mediandistance, by=c("Group.1"))
head(chr18_10k_neg_df)

# calculate the SD in each positive bin
chr18_10k_pos_sd <- aggregate(chr18_10k_pos[, 3], list(chr18_10k_pos$group), sd)
class(chr18_10k_pos_sd)
head(chr18_10k_pos_sd)
chr18_10k_pos_df <- merge(chr18_10k_pos_df, chr18_10k_pos_sd, by=c("Group.1"))
head(chr18_10k_pos_df)
chr18_10k_pos_df$n_in_group <- as.numeric(table(chr18_10k_pos$group))
head(chr18_10k_pos_df)
table(chr18_10k_pos$group)

# calculate the SD in each negative bin
chr18_10k_neg_sd <- aggregate(chr18_10k_neg[, 3], list(chr18_10k_neg$group), sd)
class(chr18_10k_neg_sd)
head(chr18_10k_neg_sd)
chr18_10k_neg_df <- merge(chr18_10k_neg_df, chr18_10k_neg_sd, by=c("Group.1"))
head(chr18_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr18_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr18_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr18_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr18_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr18_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 18: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr18_1k <- chr18_10k_sort[chr18_10k_sort$absdistance <= 1000,]
chr18_1k$group <- "group"
chr18_1k$group[chr18_1k$value >=-1 & chr18_1k$value < -0.9] <- "-1 to -0.9"
chr18_1k$group[chr18_1k$value >=-0.9 & chr18_1k$value < -0.8] <- "-0.9 to -0.8"
chr18_1k$group[chr18_1k$value >=-0.8 & chr18_1k$value < -0.7] <- "-0.8 to -0.7"
chr18_1k$group[chr18_1k$value >=-0.7 & chr18_1k$value < -0.6] <- "-0.7 to -0.6"
chr18_1k$group[chr18_1k$value >=-0.6 & chr18_1k$value < -0.5] <- "-0.6 to -0.5"
chr18_1k$group[chr18_1k$value >=-0.5 & chr18_1k$value < -0.4] <- "-0.5 to -0.4"
chr18_1k$group[chr18_1k$value >=-0.4 & chr18_1k$value < -0.3] <- "-0.4 to -0.3"
chr18_1k$group[chr18_1k$value >=-0.3 & chr18_1k$value < -0.2] <- "-0.3 to -0.2"
chr18_1k$group[chr18_1k$value >=-0.2 & chr18_1k$value < -0.1] <- "-0.2 to -0.1"
chr18_1k$group[chr18_1k$value >=-0.1 & chr18_1k$value < -0] <- "-0.1 to 0"
chr18_1k$group[chr18_1k$value >=0 & chr18_1k$value < 0.1] <- "0 to 0.1"
chr18_1k$group[chr18_1k$value >=0.1 & chr18_1k$value < 0.2] <- "0.1 to 0.2"
chr18_1k$group[chr18_1k$value >=0.2 & chr18_1k$value < 0.3] <- "0.2 to 0.3"
chr18_1k$group[chr18_1k$value >=0.3 & chr18_1k$value < 0.4] <- "0.3 to 0.4"
chr18_1k$group[chr18_1k$value >=0.4 & chr18_1k$value < 0.5] <- "0.4 to 0.5"
chr18_1k$group[chr18_1k$value >=0.5 & chr18_1k$value < 0.6] <- "0.5 to 0.6"
chr18_1k$group[chr18_1k$value >=0.6 & chr18_1k$value < 0.7] <- "0.6 to 0.7"
chr18_1k$group[chr18_1k$value >=0.7 & chr18_1k$value < 0.8] <- "0.7 to 0.8"
chr18_1k$group[chr18_1k$value >=0.8 & chr18_1k$value < 0.9] <- "0.8 to 0.9"
chr18_1k$group[chr18_1k$value >=0.9 & chr18_1k$value < 1] <- "0.9 to 1"

chr18_1k$group <- factor(chr18_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr18_1k$group)
cor.group <- as.data.frame(table(chr18_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr18_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 18: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr19 #####

load("/path/to/cis_correlations/chr19_1.Rdata")
load("/path/to/cis_correlations/chr19_2.Rdata")
load("/path/to/cis_correlations/chr19_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr19_1_dim <- lapply(chr19, FUN = getdim)
chr19_1_dim_total <- lapply(chr19_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr19_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr19_2_dim <- lapply(chr19_2, FUN = getdim)
chr19_2_dim_total <- lapply(chr19_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr19_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr19_3_dim <- lapply(chr19_3, FUN = getdim)
chr19_3_dim_total <- lapply(chr19_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr19_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr19_1 <- rbindlist(chr19)
dim(chr19_1)
chr19_2 <- rbindlist(chr19_2)
dim(chr19_2)
chr19_3 <- rbindlist(chr19_3)
dim(chr19_3)

chr19 <- rbind(chr19_1, chr19_2)
dim(chr19)
chr19 <- rbind(chr19, chr19_3)
dim(chr19)
F7_cis_dim <- list()
F7_cis_dim$chr19_F7 <- dim(chr19)

chr19_10k <- chr19[chr19$absdistance <= 10000]
dim(chr19_10k)
head(chr19_10k)
chr19_10k_sort <- chr19_10k[order(absdistance)]
head(chr19_10k_sort)

# split to positive and negative correlations
chr19_10k_pos <- chr19_10k_sort[chr19_10k_sort$value >= 0,]
dim(chr19_10k_pos)
chr19_10k_neg <- chr19_10k_sort[chr19_10k_sort$value < 0,]
dim(chr19_10k_neg)

# bin the positive correlations on distance
chr19_10k_pos$group <- as.numeric(cut2(chr19_10k_pos$absdistance, m=100))
summary(chr19_10k_pos$group)
chr19_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr19_10k_pos_valuemean <- aggregate(chr19_10k_pos[, 3], list(chr19_10k_pos$group), mean)
class(chr19_10k_pos_valuemean)
head(chr19_10k_pos_valuemean)
# calculate the median distance of each bin
chr19_10k_pos_mediandistance <- aggregate(chr19_10k_pos[, 4], list(chr19_10k_pos$group), median)
class(chr19_10k_pos_mediandistance)
head(chr19_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr19_10k_pos_df <- merge(chr19_10k_pos_valuemean, chr19_10k_pos_mediandistance, by=c("Group.1"))
head(chr19_10k_pos_df)

# bin the negative correlations on distance
chr19_10k_neg$group <- as.numeric(cut2(chr19_10k_neg$absdistance, m=100))
summary(chr19_10k_neg$group)
chr19_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr19_10k_neg_valuemean <- aggregate(chr19_10k_neg[, 3], list(chr19_10k_neg$group), mean)
class(chr19_10k_neg_valuemean)
head(chr19_10k_neg_valuemean)
# calculate the median distance of each bin
chr19_10k_neg_mediandistance <- aggregate(chr19_10k_neg[, 4], list(chr19_10k_neg$group), median)
class(chr19_10k_neg_mediandistance)
head(chr19_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr19_10k_neg_df <- merge(chr19_10k_neg_valuemean, chr19_10k_neg_mediandistance, by=c("Group.1"))
head(chr19_10k_neg_df)

# calculate the SD in each positive bin
chr19_10k_pos_sd <- aggregate(chr19_10k_pos[, 3], list(chr19_10k_pos$group), sd)
class(chr19_10k_pos_sd)
head(chr19_10k_pos_sd)
chr19_10k_pos_df <- merge(chr19_10k_pos_df, chr19_10k_pos_sd, by=c("Group.1"))
head(chr19_10k_pos_df)
chr19_10k_pos_df$n_in_group <- as.numeric(table(chr19_10k_pos$group))
head(chr19_10k_pos_df)
table(chr19_10k_pos$group)

# calculate the SD in each negative bin
chr19_10k_neg_sd <- aggregate(chr19_10k_neg[, 3], list(chr19_10k_neg$group), sd)
class(chr19_10k_neg_sd)
head(chr19_10k_neg_sd)
chr19_10k_neg_df <- merge(chr19_10k_neg_df, chr19_10k_neg_sd, by=c("Group.1"))
head(chr19_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr19_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr19_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr19_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr19_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr19_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 19: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr19_1k <- chr19_10k_sort[chr19_10k_sort$absdistance <= 1000,]
chr19_1k$group <- "group"
chr19_1k$group[chr19_1k$value >=-1 & chr19_1k$value < -0.9] <- "-1 to -0.9"
chr19_1k$group[chr19_1k$value >=-0.9 & chr19_1k$value < -0.8] <- "-0.9 to -0.8"
chr19_1k$group[chr19_1k$value >=-0.8 & chr19_1k$value < -0.7] <- "-0.8 to -0.7"
chr19_1k$group[chr19_1k$value >=-0.7 & chr19_1k$value < -0.6] <- "-0.7 to -0.6"
chr19_1k$group[chr19_1k$value >=-0.6 & chr19_1k$value < -0.5] <- "-0.6 to -0.5"
chr19_1k$group[chr19_1k$value >=-0.5 & chr19_1k$value < -0.4] <- "-0.5 to -0.4"
chr19_1k$group[chr19_1k$value >=-0.4 & chr19_1k$value < -0.3] <- "-0.4 to -0.3"
chr19_1k$group[chr19_1k$value >=-0.3 & chr19_1k$value < -0.2] <- "-0.3 to -0.2"
chr19_1k$group[chr19_1k$value >=-0.2 & chr19_1k$value < -0.1] <- "-0.2 to -0.1"
chr19_1k$group[chr19_1k$value >=-0.1 & chr19_1k$value < -0] <- "-0.1 to 0"
chr19_1k$group[chr19_1k$value >=0 & chr19_1k$value < 0.1] <- "0 to 0.1"
chr19_1k$group[chr19_1k$value >=0.1 & chr19_1k$value < 0.2] <- "0.1 to 0.2"
chr19_1k$group[chr19_1k$value >=0.2 & chr19_1k$value < 0.3] <- "0.2 to 0.3"
chr19_1k$group[chr19_1k$value >=0.3 & chr19_1k$value < 0.4] <- "0.3 to 0.4"
chr19_1k$group[chr19_1k$value >=0.4 & chr19_1k$value < 0.5] <- "0.4 to 0.5"
chr19_1k$group[chr19_1k$value >=0.5 & chr19_1k$value < 0.6] <- "0.5 to 0.6"
chr19_1k$group[chr19_1k$value >=0.6 & chr19_1k$value < 0.7] <- "0.6 to 0.7"
chr19_1k$group[chr19_1k$value >=0.7 & chr19_1k$value < 0.8] <- "0.7 to 0.8"
chr19_1k$group[chr19_1k$value >=0.8 & chr19_1k$value < 0.9] <- "0.8 to 0.9"
chr19_1k$group[chr19_1k$value >=0.9 & chr19_1k$value < 1] <- "0.9 to 1"

chr19_1k$group <- factor(chr19_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr19_1k$group)
cor.group <- as.data.frame(table(chr19_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr19_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 19: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr20 #####

load("/path/to/cis_correlations/chr20_1.Rdata")
load("/path/to/cis_correlations/chr20_2.Rdata")
load("/path/to/cis_correlations/chr20_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr20_1_dim <- lapply(chr20, FUN = getdim)
chr20_1_dim_total <- lapply(chr20_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr20_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr20_2_dim <- lapply(chr20_2, FUN = getdim)
chr20_2_dim_total <- lapply(chr20_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr20_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr20_3_dim <- lapply(chr20_3, FUN = getdim)
chr20_3_dim_total <- lapply(chr20_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr20_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr20_1 <- rbindlist(chr20)
dim(chr20_1)
chr20_2 <- rbindlist(chr20_2)
dim(chr20_2)
chr20_3 <- rbindlist(chr20_3)
dim(chr20_3)

chr20 <- rbind(chr20_1, chr20_2)
dim(chr20)
chr20 <- rbind(chr20, chr20_3)
dim(chr20)
F7_cis_dim <- list()
F7_cis_dim$chr20_F7 <- dim(chr20)

chr20_10k <- chr20[chr20$absdistance <= 10000]
dim(chr20_10k)
head(chr20_10k)
chr20_10k_sort <- chr20_10k[order(absdistance)]
head(chr20_10k_sort)

# split to positive and negative correlations
chr20_10k_pos <- chr20_10k_sort[chr20_10k_sort$value >= 0,]
dim(chr20_10k_pos)
chr20_10k_neg <- chr20_10k_sort[chr20_10k_sort$value < 0,]
dim(chr20_10k_neg)

# bin the positive correlations on distance
chr20_10k_pos$group <- as.numeric(cut2(chr20_10k_pos$absdistance, m=100))
summary(chr20_10k_pos$group)
chr20_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr20_10k_pos_valuemean <- aggregate(chr20_10k_pos[, 3], list(chr20_10k_pos$group), mean)
class(chr20_10k_pos_valuemean)
head(chr20_10k_pos_valuemean)
# calculate the median distance of each bin
chr20_10k_pos_mediandistance <- aggregate(chr20_10k_pos[, 4], list(chr20_10k_pos$group), median)
class(chr20_10k_pos_mediandistance)
head(chr20_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr20_10k_pos_df <- merge(chr20_10k_pos_valuemean, chr20_10k_pos_mediandistance, by=c("Group.1"))
head(chr20_10k_pos_df)

# bin the negative correlations on distance
chr20_10k_neg$group <- as.numeric(cut2(chr20_10k_neg$absdistance, m=100))
summary(chr20_10k_neg$group)
chr20_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr20_10k_neg_valuemean <- aggregate(chr20_10k_neg[, 3], list(chr20_10k_neg$group), mean)
class(chr20_10k_neg_valuemean)
head(chr20_10k_neg_valuemean)
# calculate the median distance of each bin
chr20_10k_neg_mediandistance <- aggregate(chr20_10k_neg[, 4], list(chr20_10k_neg$group), median)
class(chr20_10k_neg_mediandistance)
head(chr20_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr20_10k_neg_df <- merge(chr20_10k_neg_valuemean, chr20_10k_neg_mediandistance, by=c("Group.1"))
head(chr20_10k_neg_df)

# calculate the SD in each positive bin
chr20_10k_pos_sd <- aggregate(chr20_10k_pos[, 3], list(chr20_10k_pos$group), sd)
class(chr20_10k_pos_sd)
head(chr20_10k_pos_sd)
chr20_10k_pos_df <- merge(chr20_10k_pos_df, chr20_10k_pos_sd, by=c("Group.1"))
head(chr20_10k_pos_df)
chr20_10k_pos_df$n_in_group <- as.numeric(table(chr20_10k_pos$group))
head(chr20_10k_pos_df)
table(chr20_10k_pos$group)

# calculate the SD in each negative bin
chr20_10k_neg_sd <- aggregate(chr20_10k_neg[, 3], list(chr20_10k_neg$group), sd)
class(chr20_10k_neg_sd)
head(chr20_10k_neg_sd)
chr20_10k_neg_df <- merge(chr20_10k_neg_df, chr20_10k_neg_sd, by=c("Group.1"))
head(chr20_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr20_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr20_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr20_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr20_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr20_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 20: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr20_1k <- chr20_10k_sort[chr20_10k_sort$absdistance <= 1000,]
chr20_1k$group <- "group"
chr20_1k$group[chr20_1k$value >=-1 & chr20_1k$value < -0.9] <- "-1 to -0.9"
chr20_1k$group[chr20_1k$value >=-0.9 & chr20_1k$value < -0.8] <- "-0.9 to -0.8"
chr20_1k$group[chr20_1k$value >=-0.8 & chr20_1k$value < -0.7] <- "-0.8 to -0.7"
chr20_1k$group[chr20_1k$value >=-0.7 & chr20_1k$value < -0.6] <- "-0.7 to -0.6"
chr20_1k$group[chr20_1k$value >=-0.6 & chr20_1k$value < -0.5] <- "-0.6 to -0.5"
chr20_1k$group[chr20_1k$value >=-0.5 & chr20_1k$value < -0.4] <- "-0.5 to -0.4"
chr20_1k$group[chr20_1k$value >=-0.4 & chr20_1k$value < -0.3] <- "-0.4 to -0.3"
chr20_1k$group[chr20_1k$value >=-0.3 & chr20_1k$value < -0.2] <- "-0.3 to -0.2"
chr20_1k$group[chr20_1k$value >=-0.2 & chr20_1k$value < -0.1] <- "-0.2 to -0.1"
chr20_1k$group[chr20_1k$value >=-0.1 & chr20_1k$value < -0] <- "-0.1 to 0"
chr20_1k$group[chr20_1k$value >=0 & chr20_1k$value < 0.1] <- "0 to 0.1"
chr20_1k$group[chr20_1k$value >=0.1 & chr20_1k$value < 0.2] <- "0.1 to 0.2"
chr20_1k$group[chr20_1k$value >=0.2 & chr20_1k$value < 0.3] <- "0.2 to 0.3"
chr20_1k$group[chr20_1k$value >=0.3 & chr20_1k$value < 0.4] <- "0.3 to 0.4"
chr20_1k$group[chr20_1k$value >=0.4 & chr20_1k$value < 0.5] <- "0.4 to 0.5"
chr20_1k$group[chr20_1k$value >=0.5 & chr20_1k$value < 0.6] <- "0.5 to 0.6"
chr20_1k$group[chr20_1k$value >=0.6 & chr20_1k$value < 0.7] <- "0.6 to 0.7"
chr20_1k$group[chr20_1k$value >=0.7 & chr20_1k$value < 0.8] <- "0.7 to 0.8"
chr20_1k$group[chr20_1k$value >=0.8 & chr20_1k$value < 0.9] <- "0.8 to 0.9"
chr20_1k$group[chr20_1k$value >=0.9 & chr20_1k$value < 1] <- "0.9 to 1"

chr20_1k$group <- factor(chr20_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr20_1k$group)
cor.group <- as.data.frame(table(chr20_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr20_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 20: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr21 #####

load("/path/to/cis_correlations/chr21_1.Rdata")
load("/path/to/cis_correlations/chr21_2.Rdata")
load("/path/to/cis_correlations/chr21_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr21_1_dim <- lapply(chr21, FUN = getdim)
chr21_1_dim_total <- lapply(chr21_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr21_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr21_2_dim <- lapply(chr21_2, FUN = getdim)
chr21_2_dim_total <- lapply(chr21_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr21_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr21_3_dim <- lapply(chr21_3, FUN = getdim)
chr21_3_dim_total <- lapply(chr21_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr21_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr21_1 <- rbindlist(chr21)
dim(chr21_1)
chr21_2 <- rbindlist(chr21_2)
dim(chr21_2)
chr21_3 <- rbindlist(chr21_3)
dim(chr21_3)

chr21 <- rbind(chr21_1, chr21_2)
dim(chr21)
chr21 <- rbind(chr21, chr21_3)
dim(chr21)
F7_cis_dim <- list()
F7_cis_dim$chr21_F7 <- dim(chr21)

chr21_10k <- chr21[chr21$absdistance <= 10000]
dim(chr21_10k)
head(chr21_10k)
chr21_10k_sort <- chr21_10k[order(absdistance)]
head(chr21_10k_sort)

# split to positive and negative correlations
chr21_10k_pos <- chr21_10k_sort[chr21_10k_sort$value >= 0,]
dim(chr21_10k_pos)
chr21_10k_neg <- chr21_10k_sort[chr21_10k_sort$value < 0,]
dim(chr21_10k_neg)

# bin the positive correlations on distance
chr21_10k_pos$group <- as.numeric(cut2(chr21_10k_pos$absdistance, m=100))
summary(chr21_10k_pos$group)
chr21_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr21_10k_pos_valuemean <- aggregate(chr21_10k_pos[, 3], list(chr21_10k_pos$group), mean)
class(chr21_10k_pos_valuemean)
head(chr21_10k_pos_valuemean)
# calculate the median distance of each bin
chr21_10k_pos_mediandistance <- aggregate(chr21_10k_pos[, 4], list(chr21_10k_pos$group), median)
class(chr21_10k_pos_mediandistance)
head(chr21_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr21_10k_pos_df <- merge(chr21_10k_pos_valuemean, chr21_10k_pos_mediandistance, by=c("Group.1"))
head(chr21_10k_pos_df)

# bin the negative correlations on distance
chr21_10k_neg$group <- as.numeric(cut2(chr21_10k_neg$absdistance, m=100))
summary(chr21_10k_neg$group)
chr21_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr21_10k_neg_valuemean <- aggregate(chr21_10k_neg[, 3], list(chr21_10k_neg$group), mean)
class(chr21_10k_neg_valuemean)
head(chr21_10k_neg_valuemean)
# calculate the median distance of each bin
chr21_10k_neg_mediandistance <- aggregate(chr21_10k_neg[, 4], list(chr21_10k_neg$group), median)
class(chr21_10k_neg_mediandistance)
head(chr21_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr21_10k_neg_df <- merge(chr21_10k_neg_valuemean, chr21_10k_neg_mediandistance, by=c("Group.1"))
head(chr21_10k_neg_df)

# calculate the SD in each positive bin
chr21_10k_pos_sd <- aggregate(chr21_10k_pos[, 3], list(chr21_10k_pos$group), sd)
class(chr21_10k_pos_sd)
head(chr21_10k_pos_sd)
chr21_10k_pos_df <- merge(chr21_10k_pos_df, chr21_10k_pos_sd, by=c("Group.1"))
head(chr21_10k_pos_df)
chr21_10k_pos_df$n_in_group <- as.numeric(table(chr21_10k_pos$group))
head(chr21_10k_pos_df)
table(chr21_10k_pos$group)

# calculate the SD in each negative bin
chr21_10k_neg_sd <- aggregate(chr21_10k_neg[, 3], list(chr21_10k_neg$group), sd)
class(chr21_10k_neg_sd)
head(chr21_10k_neg_sd)
chr21_10k_neg_df <- merge(chr21_10k_neg_df, chr21_10k_neg_sd, by=c("Group.1"))
head(chr21_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr21_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr21_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr21_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr21_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr21_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 21: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr21_1k <- chr21_10k_sort[chr21_10k_sort$absdistance <= 1000,]
chr21_1k$group <- "group"
chr21_1k$group[chr21_1k$value >=-1 & chr21_1k$value < -0.9] <- "-1 to -0.9"
chr21_1k$group[chr21_1k$value >=-0.9 & chr21_1k$value < -0.8] <- "-0.9 to -0.8"
chr21_1k$group[chr21_1k$value >=-0.8 & chr21_1k$value < -0.7] <- "-0.8 to -0.7"
chr21_1k$group[chr21_1k$value >=-0.7 & chr21_1k$value < -0.6] <- "-0.7 to -0.6"
chr21_1k$group[chr21_1k$value >=-0.6 & chr21_1k$value < -0.5] <- "-0.6 to -0.5"
chr21_1k$group[chr21_1k$value >=-0.5 & chr21_1k$value < -0.4] <- "-0.5 to -0.4"
chr21_1k$group[chr21_1k$value >=-0.4 & chr21_1k$value < -0.3] <- "-0.4 to -0.3"
chr21_1k$group[chr21_1k$value >=-0.3 & chr21_1k$value < -0.2] <- "-0.3 to -0.2"
chr21_1k$group[chr21_1k$value >=-0.2 & chr21_1k$value < -0.1] <- "-0.2 to -0.1"
chr21_1k$group[chr21_1k$value >=-0.1 & chr21_1k$value < -0] <- "-0.1 to 0"
chr21_1k$group[chr21_1k$value >=0 & chr21_1k$value < 0.1] <- "0 to 0.1"
chr21_1k$group[chr21_1k$value >=0.1 & chr21_1k$value < 0.2] <- "0.1 to 0.2"
chr21_1k$group[chr21_1k$value >=0.2 & chr21_1k$value < 0.3] <- "0.2 to 0.3"
chr21_1k$group[chr21_1k$value >=0.3 & chr21_1k$value < 0.4] <- "0.3 to 0.4"
chr21_1k$group[chr21_1k$value >=0.4 & chr21_1k$value < 0.5] <- "0.4 to 0.5"
chr21_1k$group[chr21_1k$value >=0.5 & chr21_1k$value < 0.6] <- "0.5 to 0.6"
chr21_1k$group[chr21_1k$value >=0.6 & chr21_1k$value < 0.7] <- "0.6 to 0.7"
chr21_1k$group[chr21_1k$value >=0.7 & chr21_1k$value < 0.8] <- "0.7 to 0.8"
chr21_1k$group[chr21_1k$value >=0.8 & chr21_1k$value < 0.9] <- "0.8 to 0.9"
chr21_1k$group[chr21_1k$value >=0.9 & chr21_1k$value < 1] <- "0.9 to 1"

chr21_1k$group <- factor(chr21_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr21_1k$group)
cor.group <- as.data.frame(table(chr21_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr21_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 21: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()




###### chr22 #####

load("/path/to/cis_correlations/chr22_1.Rdata")
load("/path/to/cis_correlations/chr22_2.Rdata")
load("/path/to/cis_correlations/chr22_3.Rdata")

ls()

getdim <- function(x){
  print(dim(x))
  return(dim(x))
}
chr22_1_dim <- lapply(chr22, FUN = getdim)
chr22_1_dim_total <- lapply(chr22_1_dim, FUN=function(x) x[1])
total_sum_1_1 <- 0
for (i in chr22_1_dim_total){
  total_sum_1_1 <- total_sum_1_1 + i
}

chr22_2_dim <- lapply(chr22_2, FUN = getdim)
chr22_2_dim_total <- lapply(chr22_2_dim, FUN=function(x) x[1])
total_sum_1_2 <- 0
for (i in chr22_2_dim_total){
  total_sum_1_2 <- total_sum_1_2 + i
}

chr22_3_dim <- lapply(chr22_3, FUN = getdim)
chr22_3_dim_total <- lapply(chr22_3_dim, FUN=function(x) x[1])
total_sum_1_3 <- 0
for (i in chr22_3_dim_total){
  total_sum_1_3 <- total_sum_1_3 + i
}

chr22_1 <- rbindlist(chr22)
dim(chr22_1)
chr22_2 <- rbindlist(chr22_2)
dim(chr22_2)
chr22_3 <- rbindlist(chr22_3)
dim(chr22_3)

chr22 <- rbind(chr22_1, chr22_2)
dim(chr22)
chr22 <- rbind(chr22, chr22_3)
dim(chr22)
F7_cis_dim <- list()
F7_cis_dim$chr22_F7 <- dim(chr22)

chr22_10k <- chr22[chr22$absdistance <= 10000]
dim(chr22_10k)
head(chr22_10k)
chr22_10k_sort <- chr22_10k[order(absdistance)]
head(chr22_10k_sort)

# split to positive and negative correlations
chr22_10k_pos <- chr22_10k_sort[chr22_10k_sort$value >= 0,]
dim(chr22_10k_pos)
chr22_10k_neg <- chr22_10k_sort[chr22_10k_sort$value < 0,]
dim(chr22_10k_neg)

# bin the positive correlations on distance
chr22_10k_pos$group <- as.numeric(cut2(chr22_10k_pos$absdistance, m=100))
summary(chr22_10k_pos$group)
chr22_10k_pos$group[1:900]
# calculate the mean correlation in each bin
chr22_10k_pos_valuemean <- aggregate(chr22_10k_pos[, 3], list(chr22_10k_pos$group), mean)
class(chr22_10k_pos_valuemean)
head(chr22_10k_pos_valuemean)
# calculate the median distance of each bin
chr22_10k_pos_mediandistance <- aggregate(chr22_10k_pos[, 4], list(chr22_10k_pos$group), median)
class(chr22_10k_pos_mediandistance)
head(chr22_10k_pos_mediandistance)
# merge together the mean correlation and median distance
chr22_10k_pos_df <- merge(chr22_10k_pos_valuemean, chr22_10k_pos_mediandistance, by=c("Group.1"))
head(chr22_10k_pos_df)

# bin the negative correlations on distance
chr22_10k_neg$group <- as.numeric(cut2(chr22_10k_neg$absdistance, m=100))
summary(chr22_10k_neg$group)
chr22_10k_neg$group[1:900]
# calculate the mean correlation in each bin
chr22_10k_neg_valuemean <- aggregate(chr22_10k_neg[, 3], list(chr22_10k_neg$group), mean)
class(chr22_10k_neg_valuemean)
head(chr22_10k_neg_valuemean)
# calculate the median distance of each bin
chr22_10k_neg_mediandistance <- aggregate(chr22_10k_neg[, 4], list(chr22_10k_neg$group), median)
class(chr22_10k_neg_mediandistance)
head(chr22_10k_neg_mediandistance)
# merge together the mean correlation and median distance
chr22_10k_neg_df <- merge(chr22_10k_neg_valuemean, chr22_10k_neg_mediandistance, by=c("Group.1"))
head(chr22_10k_neg_df)

# calculate the SD in each positive bin
chr22_10k_pos_sd <- aggregate(chr22_10k_pos[, 3], list(chr22_10k_pos$group), sd)
class(chr22_10k_pos_sd)
head(chr22_10k_pos_sd)
chr22_10k_pos_df <- merge(chr22_10k_pos_df, chr22_10k_pos_sd, by=c("Group.1"))
head(chr22_10k_pos_df)
chr22_10k_pos_df$n_in_group <- as.numeric(table(chr22_10k_pos$group))
head(chr22_10k_pos_df)
table(chr22_10k_pos$group)

# calculate the SD in each negative bin
chr22_10k_neg_sd <- aggregate(chr22_10k_neg[, 3], list(chr22_10k_neg$group), sd)
class(chr22_10k_neg_sd)
head(chr22_10k_neg_sd)
chr22_10k_neg_df <- merge(chr22_10k_neg_df, chr22_10k_neg_sd, by=c("Group.1"))
head(chr22_10k_neg_df)

pdf(file = "/path/to/F7_cisdecayplots/chr22_decayplot_F7_10k_SD.pdf", width = 7, height = 5)
ggplot()+
  geom_errorbar(data=chr22_10k_pos_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr22_10k_pos_df,aes(x=absdistance,y=value.x,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr22_10k_neg_df,aes(x=absdistance,ymin=value.x-value.y, ymax=value.x+value.y), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr22_10k_neg_df,aes(x=absdistance,y=value.x,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 22: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# histogram to show 1k correlation values

chr22_1k <- chr22_10k_sort[chr22_10k_sort$absdistance <= 1000,]
chr22_1k$group <- "group"
chr22_1k$group[chr22_1k$value >=-1 & chr22_1k$value < -0.9] <- "-1 to -0.9"
chr22_1k$group[chr22_1k$value >=-0.9 & chr22_1k$value < -0.8] <- "-0.9 to -0.8"
chr22_1k$group[chr22_1k$value >=-0.8 & chr22_1k$value < -0.7] <- "-0.8 to -0.7"
chr22_1k$group[chr22_1k$value >=-0.7 & chr22_1k$value < -0.6] <- "-0.7 to -0.6"
chr22_1k$group[chr22_1k$value >=-0.6 & chr22_1k$value < -0.5] <- "-0.6 to -0.5"
chr22_1k$group[chr22_1k$value >=-0.5 & chr22_1k$value < -0.4] <- "-0.5 to -0.4"
chr22_1k$group[chr22_1k$value >=-0.4 & chr22_1k$value < -0.3] <- "-0.4 to -0.3"
chr22_1k$group[chr22_1k$value >=-0.3 & chr22_1k$value < -0.2] <- "-0.3 to -0.2"
chr22_1k$group[chr22_1k$value >=-0.2 & chr22_1k$value < -0.1] <- "-0.2 to -0.1"
chr22_1k$group[chr22_1k$value >=-0.1 & chr22_1k$value < -0] <- "-0.1 to 0"
chr22_1k$group[chr22_1k$value >=0 & chr22_1k$value < 0.1] <- "0 to 0.1"
chr22_1k$group[chr22_1k$value >=0.1 & chr22_1k$value < 0.2] <- "0.1 to 0.2"
chr22_1k$group[chr22_1k$value >=0.2 & chr22_1k$value < 0.3] <- "0.2 to 0.3"
chr22_1k$group[chr22_1k$value >=0.3 & chr22_1k$value < 0.4] <- "0.3 to 0.4"
chr22_1k$group[chr22_1k$value >=0.4 & chr22_1k$value < 0.5] <- "0.4 to 0.5"
chr22_1k$group[chr22_1k$value >=0.5 & chr22_1k$value < 0.6] <- "0.5 to 0.6"
chr22_1k$group[chr22_1k$value >=0.6 & chr22_1k$value < 0.7] <- "0.6 to 0.7"
chr22_1k$group[chr22_1k$value >=0.7 & chr22_1k$value < 0.8] <- "0.7 to 0.8"
chr22_1k$group[chr22_1k$value >=0.8 & chr22_1k$value < 0.9] <- "0.8 to 0.9"
chr22_1k$group[chr22_1k$value >=0.9 & chr22_1k$value < 1] <- "0.9 to 1"

chr22_1k$group <- factor(chr22_1k$group, 
                         levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                  "-0.8 to -0.7", "-0.7 to -0.6",
                                  "-0.6 to -0.5", "-0.5 to -0.4",
                                  "-0.4 to -0.3", "-0.3 to -0.2",
                                  "-0.2 to -0.1", "-0.1 to 0",
                                  "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                  "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                  "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                  "0.9 to 1"))


table(chr22_1k$group)
cor.group <- as.data.frame(table(chr22_1k$group))
cor.group

pdf(file = "/path/to/F7_cisdecayplots/chr22_histogram_F7_1k.pdf", width = 7, height = 5)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 22: values of cis correlations\nwithin 1kb in ARIES 7 year olds"), x = "Correlation", y = "Frequency")
dev.off()