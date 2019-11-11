# kTotal for calculating DNAm connectivity

library(WGCNA)
library(ggplot2)
library(meffil)
library(data.table)
library(dplyr)

load("/path/to/F7readyfor450kcor.Rdata")
probeDetails <- meffil.featureset("450k")
rownames(probeDetails) <- probeDetails$name
head(probeDetails)
probeDetails <- probeDetails[probeDetails$chromosome == "chr20",]
dim(probeDetails)
chr20 <- as.character(probeDetails$name)
F7Data.sort <- F7Data.sort[,colnames(F7Data.sort) %in% chr20]
dim(F7Data.sort)
F7Data.sort[1:5,1:5]
F7_chr20_cormat <- bicor(F7Data.sort, maxPOutliers = 0.05)
F7_chr20_cormat[lower.tri(F7_chr20_cormat, diag = TRUE)] <- NA
F7_chr20_cormat[1:5,1:5]

kTotal = apply(F7_chr20_cormat, 2, sum, na.rm = TRUE)
class(kTotal)
length(kTotal)
kTotal[1:10]
summary(kTotal)

kTotal.df <- as.data.frame(kTotal)
jpeg(filename = "~/F7/chr20_kTotal_F7_hist.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot() + 
  geom_histogram(data=kTotal.df,aes(x=kTotal),color="black", fill="#238A8DFF")+
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  labs(title=paste("Chromosome 20: histogram of kTotal scores, in ARIES 7 year olds"), x = "kTotal score", y = "Frequency")
dev.off()

# getting alspac kids for chr20

library(data.table)
chr20.fam <- fread("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2017-05-04/data/plink/data_chr20.fam")
chr20.kids <- chr20.fam[substr(chr20.fam$V1,nchar(chr20.fam$V1),nchar(chr20.fam$V1)) == "A",]
chr20.kids <- chr20.kids[,c(1:2)]
write.table(chr20.kids, file = "~/F7/chr20kidsIDs.txt", row.names = F, col.names = F, quote = F)

# then run the LDscore_ALSPACchr20.txt file to get LD score for chr 20 SNPs

# and read in the LD scores:
chr20_LD = fread("~/F7/F7_chr20_LDscores.l2.ldscore")
head(chr20_LD)
dim(chr20_LD)

# read in the godmc assoc file and reduce to sites on chr20 
assoc_file <- "~/godmc/assoc_meta_all.csv"
cpg_file <- "~/F7/F7readyfor450kcor.Rdata"
snp_file <- "~/godmc/ARIES_F7_chr20freq.raw.gz"
freq_file <- "~/godmc/ARIESchr20_freq.frq.gz"

message("Reading CPG data")
load(cpg_file)
cpg <- F7Data.sort
cpg[1:5,1:5]
# reduce to chr20
library(meffil)
probeDetails <- meffil.featureset("450k")
rownames(probeDetails) <- probeDetails$name
head(probeDetails)
probeDetails <- probeDetails[probeDetails$chromosome == "chr20",]
dim(probeDetails)
chr20 <- as.character(probeDetails$name)
cpg <- cpg[,colnames(cpg) %in% chr20]
dim(cpg)
cpg[1:5,1:5]

message("Reading association list")
assoc <- fread(assoc_file)
dim(assoc)
# reduce to cis mQTLs only
assoc <- assoc[assoc$cistrans == TRUE,]
dim(assoc)

# SNP
message("Reading SNP data")
snp <- fread(paste0("zcat ", snp_file))
snp[1:5,1:10]
message("Organising SNP data")

message("Reading SNP frequencies")
# for ALSPAC 
freq <- fread(paste0("zcat ", freq_file))
freq[1:30,]


iid <- snp$IID
snp <- snp[, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"):=NULL]
snp[1:5,1:20]

snp <- as.matrix(snp)
rownames(snp) <- iid
snp[1:5,1:5]
# now take off the A so the ARIES ID matches DNAm
iid[1:50]
iid2 <- substr(iid, 1, nchar(iid)-1)
iid2[1:50]
rownames(snp) <- iid2
snp[1:5,1:5]
class(snp[,1])

# create snpinfo which will contain info on allele freqs etc
snpinfo <- as.character(colnames(snp))
snpinfo <- as.data.frame(snpinfo)
names(snpinfo) <- c("full")
head(snpinfo)
snpinfo$full <- as.character(snpinfo$full)
head(snpinfo)
# create EA column from last character of SNP name
snpinfo$SNP <- substr(snpinfo$full, 1, nchar(snpinfo$full)-2)
snpinfo$EA <- substr(snpinfo$full, nchar(snpinfo$full), nchar(snpinfo$full))
# check that the EA is the same in both freq and snpinfo
snpinfo[1:30,]
freq[1:30,]
snpinfo <- snpinfo[,c(2:3)]
head(snpinfo)
dim(snpinfo)
dim(freq)

colnames(snp) <- snpinfo$SNP
snp[1:5,1:10]
snpinfo <- merge(snpinfo, freq, by.x="SNP", all.x=TRUE)
head(snpinfo)
dim(snpinfo)
class(snpinfo$EA)
class(snpinfo$A1)

# where don't EA and A1 match?
snpinfo_nomatch <- snpinfo[! snpinfo$EA == snpinfo$A1,]
dim(snpinfo_nomatch)
snpinfo_nomatch[1:10,]
snpinfo_remove <- as.character(snpinfo_nomatch$SNP)
snpinfo <- snpinfo[! snpinfo$SNP %in% snpinfo_remove,]
dim(snpinfo)
stopifnot(all(snpinfo$EA == snpinfo$A1))
snpinfo <- dplyr::select(snpinfo, MARKERNAME=SNP, EA=EA, NEA=A2, EAF=MAF)
head(snpinfo)

# CPG

message("Organising CpG data")

dim(cpg)
cpg[1:5,1:5]
# (rows must be samples)
snp_ids <- as.character(rownames(snp))
length(snp_ids)
cpg_ids <- as.character(rownames(cpg))
length(cpg_ids)
matching_ids <- intersect(snp_ids, cpg_ids)
length(matching_ids)
# there are some DNAm samples that do not seem to have genetic data, so:
snp <- snp[rownames(snp) %in% matching_ids, ]
dim(snp)
cpg <- cpg[rownames(cpg) %in% matching_ids, ]
dim(cpg)
stopifnot(all(rownames(cpg) %in% rownames(snp)))
stopifnot(all(rownames(snp) %in% rownames(cpg)))

cpg <- cpg[match(rownames(snp), rownames(cpg)), ]
cpg[1:5,1:5]
rownames(snp)[1:5]
stopifnot(all(rownames(cpg) == rownames(snp)))


# ASSOC 

message("Organising association lists")
head(assoc)
dim(assoc)
# read in a list of SNPs and their rsids - assoc has chr:pos format
allSNPs <- fread("~/snps.csv", header=T)
head(allSNPs)
dim(allSNPs)
allSNPs <- allSNPs[ , -c(5:15), with = FALSE]
head(allSNPs)
assoc$snp_allSNPs <- allSNPs$rsid[match(assoc$snp, allSNPs$name)]
head(assoc)
class(assoc$cpg)
cpgnames <- as.character(colnames(cpg))
assoc <- assoc[assoc$cpg %in% cpgnames,]
dim(assoc)
class(assoc$snp)
snpnames <- as.character(colnames(snp))
assoc <- assoc[assoc$snp_allSNPs %in% snpnames,]
dim(assoc)
assoc$x <- match(assoc$snp_allSNPs, colnames(snp)) # matching snps to names of snp file
assoc$y <- match(assoc$cpg, colnames(cpg)) # matching cpgs to betas file

all(colnames(snp)[assoc$x] == assoc$snp_allSNPs)
all(colnames(cpg)[assoc$y] == assoc$cpg)
head(assoc)

message("\n\n")
message("Number of individuals: ", nrow(cpg))
message("Number of CpGs: ", ncol(cpg))
message("Number of SNPs: ", ncol(snp))
message("Number of associations: ", nrow(assoc))
message("\n\n")

# Get the best cis SNP for each CPG

message("Identifying best cis-SNPs for each CpG")

temp <- group_by(assoc, cpg) %>%
  slice(which.min(pval_are)) %>%
  dplyr::select(cpg, cissnp=snp_allSNPs)

head(temp)
dim(temp)
temp <- as.data.frame(temp)
head(temp)
dim(temp)

chr20_LD <- as.data.frame(chr20_LD)
temp$LDscore <- chr20_LD$L2[match(temp$cissnp,chr20_LD$SNP)]
class(temp$LDscore)
temp$LDscore <- as.numeric(temp$LDscore)
head(temp)
head(kTotal.df)
temp$kTotal <- kTotal.df$kTotal[match(temp$cpg,rownames(kTotal.df))]
head(temp)
save(temp, file="~/F7/temp_F7_LDandkTotal.Rdata")

reg<-lm(LDscore ~ kTotal, data = temp)
reg
coeff=coefficients(reg)
coeff
eq = paste0("y = ", round(coeff[2],4), "*x + ", round(coeff[1],4))
eq

jpeg(filename = "~/F7/chr20_kTotalvsLD_F7_smooth.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot() + 
  geom_point(data=temp,aes(x=kTotal, y=LDscore),color="#440154FF", alpha=0.4)+
  geom_smooth(data=temp,aes(x=kTotal, y=LDscore))+
  labs(title=paste("Chromosome 20: scatter of cpg connectivity vs LD of best\ncis SNP, in ARIES 7 year olds. Regression slope =",eq), x = "kTotal score", y = "Best cis SNP LD score")
dev.off()

# MASS will allow us to calculate density, so that we can plot density on
# the ggplot scatter plots
library(MASS)
library(viridis)
# Density function for ggplot2: 
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
temp$density <- get_density(temp$kTotal, temp$LDscore)


pdf(file = "~/F7/chr20_kTotalvsLD_F7_smooth_density.pdf", width = 7, height = 5)
ggplot() + 
  geom_point(data=temp,aes(x=kTotal, y=LDscore, color = density), alpha=0.4)+
  scale_color_viridis()+
  labs(color="Density")+
  geom_smooth(data=temp,aes(x=kTotal, y=LDscore))+
  labs(title=paste("Chromosome 20: scatter of cpg connectivity vs LD of best\ncis SNP, in ARIES 7 year olds. Regression slope =",eq), x = "kTotal score", y = "Best cis SNP LD score")
dev.off()

highscores <- temp[temp$LDscore > 400 & temp$kTotal > 500,]
dim(highscores)
head(highscores)
highscores$cpg.pos <- probeDetails$position[match(highscores$cpg, rownames(probeDetails))]
save(highscores, file="~/LD_kTotal_highscores.Rdata")

kTotalhighscores <- temp[temp$kTotal > 750,]
dim(kTotalhighscores)
head(kTotalhighscores)
kTotalhighscores$cpg.pos <- probeDetails$position[match(kTotalhighscores$cpg, rownames(probeDetails))]
save(kTotalhighscores, file="~/kTotal_highscores.Rdata")

LDhighscores <- temp[temp$LDscore > 550,]
dim(LDhighscores)
head(LDhighscores)
LDhighscores$cpg.pos <- probeDetails$position[match(LDhighscores$cpg, rownames(probeDetails))]
save(LDhighscores, file="~/LD_highscores.Rdata")
