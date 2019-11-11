library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(biovizBase)
library(meffil)
library(dplyr)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

# create list of all probes in the analysis:
load("/path/to/F7readyfor450kcor.Rdata")
allcgs <- as.character(colnames(F7Data.sort))
length(allcgs)
load("/path/to/F7_point9.Rdata")
load("/path/to/probeDetails_meffil.Rdata")
class(probeDetails$chromosome)
# character
table(probeDetails$chromosome)

dat$Var1.position <- probeDetails$position[match(dat$Var1, rownames(probeDetails))]
dat$Var2.position <- probeDetails$position[match(dat$Var2, rownames(probeDetails))]
dat$absdistance <- abs(dat$Var1.position - dat$Var2.position)
dat$Var1.chromosome <- probeDetails$chromosome[match(dat$Var1, rownames(probeDetails))]
dat$Var2.chromosome <- probeDetails$chromosome[match(dat$Var2, rownames(probeDetails))]
dat$cis_or_trans <- as.character(c("trans"))
dat$cis_or_trans[dat$absdistance <= 1000000 & dat$Var1.chromosome == dat$Var2.chromosome] <- "cis"
print(head(dat))

trans_chr <- dat[dat$cis_or_trans == "trans",]
dim(trans_chr)
trans_probes1 <- as.character(trans_chr$Var1)
trans_probes2 <- as.character(trans_chr$Var2)
trans_probes <- c(trans_probes1, trans_probes2)
trans_probes <- unique(trans_probes)
overpoint9 <- trans_probes
length(overpoint9)

##load cell type conversion and colors
cellType_conversions=fread("/path/to/CellTypes.tsv",drop="collection")
colors=fread("/path/to/color.tsv")

setwd("/path/to/LOLA")

##load regiondb - you can swap the enrichment set here by specifying a different collection
regionDB <- loadRegionDB("/path/to/LOLACore/hg19", collection="encode_tfbs")
head(regionDB)

retaincpg<-allcgs

#get 450k locations into a data table
Illumina450 <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
Illumina450_dt=as.data.table(Illumina450)
Illumina450_dt[,cpgID:=row.names(Illumina450),]
Illumina450_dt <- Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]
Illumina450_dt[,cpgstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
Illumina450_dt[,cpgend_pre:=ifelse(strand=="-",pos+500,pos+501),]
data <- Illumina450_dt[Illumina450_dt$cpgID%in%overpoint9,]

#collapse overlaps
gr_range = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre)))
gr_cpg = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

Illumina450_dt[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,cpgend:=end(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,NsubjectHits:=overlap_red$NsubjectHits]

head(data)
##prepare GoDMC data: merge CpGs and SNPs that are in proximity to eachother to avoid infalting the results, 1kb around cp
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
setnames(data,c("cpgID"),c("cpg"))
head(Illumina450_sub)
head(data)
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==pos,])

data[,chr:=gsub("23","X",chr),]
data[,chr:=gsub("24","Y",chr),]
####run with external background for CpGs

hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)

#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~3000)
Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

# add GC, CpG frequency of probes to the data table. Also add wheher or not the
# probe is in the list of top correlated probes.
Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isHighCor:=ifelse(cpgID%in%overpoint9,TRUE,FALSE),]
Illumina450_dt<-Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]

# Create quantiles of GC content and CpG frequency
Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$quantiles<-paste(Illumina450_dt$CpG_freqquantile,Illumina450_dt$GC_freqquantile)

#ALL
# create table of the proportions of probes overall in the two quantiles, and the 
# proportion of probes which are highly correlating.
t<-table(Illumina450_dt$quantiles,Illumina450_dt$isHighCor)
t
w<-which(t[,2]<5)
w
pr<-data.frame(t[-w,1]/t[-w,2])
pr
m<-min(pr)*t[,2]
m

controllist<-list()
for (i in 1:dim(t)[1]){
  cat(i,"\n")
  if(t[i,1]>5&t[i,2]>0){
    subgroup<-Illumina450_dt[which(Illumina450_dt$quantiles==row.names(t)[i]&Illumina450_dt$isHighCor==FALSE& Illumina450_dt$chr!="chrY"),]
    id<-subgroup[sample(nrow(subgroup), size=round(m[i],0), replace=FALSE),"cpgID"]
    controllist[[i]]<-id
  }}

# Create two data tables, one containing the matched control probes and one containing
# the highly correlating probes
bg.matched<-do.call("rbind",controllist)
Illumina450_bg_matched<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched$cpgID),]
Illumina450_F7spearman<-Illumina450_dt[which(Illumina450_dt$cpgID%in%overpoint9),]

# this is the matched control probes and the highly correlating probes together 
bg_matched<-rbind(Illumina450_bg_matched,Illumina450_F7spearman)

F7_cpg_gr=unique(with(data,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart_pre, end=cpgend_pre),strand=Rle("*"))))

head(F7_cpg_gr)

#plot CG and CpG frequency for GoDMC cpgs and background 
jpeg("~/LOLA/GCfreq_F7_trans_spearman.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(Illumina450_dt,aes(x=GC_freq,col=isHighCor))+geom_density()
dev.off()
#ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
jpeg("~/LOLA/CpGfreq_F7_trans_spearman.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(Illumina450_dt,aes(x=CpG_freq,col=isHighCor))+geom_density()
dev.off()
#ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()


jpeg("~/LOLA/GCfreq_F7_trans_spearmanMATCHED.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(bg_matched,aes(x=GC_freq,col=isHighCor))+geom_density()
dev.off()
#ggplot(bg_matched,aes(x=GC_freq))+geom_density()
jpeg("~/LOLA/CpGfreq_F7_trans_spearmanMATCHED.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(bg_matched,aes(x=CpG_freq,col=isHighCor))+geom_density()
#ggplot(bg_matched,aes(x=CpG_freq))+geom_density()
dev.off()


#use all for background
Illumina450_bg=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

cpg_bg_gr_matched=unique(c(Illumina450_bg_matched,F7_cpg_gr))

# run LOLA 
lola_res0_matched=runLOLA(F7_cpg_gr, cpg_bg_gr_matched, regionDB, cores=1)

lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
save(lola_res0_matched,file="~/LOLA/lola_encodetfbs_F7_trans_mqtlcpg_overpoint9.Rdata")

lola_res0_matched$logOddsRatio
# group all the matching TFBS together
lola_res0_matched$tf <- sapply(strsplit(lola_res0_matched$antibody,"_"), `[`, 1)

pdf(paste0("~/LOLA/F7_trans_overpoint9_encodetfbs.pdf"),height=12,width=24)
pl3=ggplot(lola_res0_matched,aes(x=tf,y=logOddsRatio))+
  geom_hline(yintercept=0, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=cellType,size=pValueLog))+
  facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  guides(fill = guide_legend(ncol=20))+
  theme(legend.text=element_text(size=12)) +
  labs(y="Odds ratio (log scale)",x="Encode segmentation") +
  scale_fill_brewer(type="qual")

print(pl3)
dev.off()

