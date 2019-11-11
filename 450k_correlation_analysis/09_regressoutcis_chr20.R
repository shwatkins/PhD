library(WGCNA)
library(stringr)
library(data.table)
library(ggplot2)
library(Hmisc)
library(reshape2)

load("/path/to/F7_cissnp_regressionresiduals.Rdata")
# construct chr20 correlation matrix
chr20cormat <- bicor(regression_residuals, maxPOutliers = 0.05)
dim(chr20cormat)
chr20cormat[lower.tri(chr20cormat, diag = TRUE)] <- NA
chr20cormat[1:5,1:5]
# melt to df with cpgs and values

chr20_melt <- melt(as.matrix(chr20cormat))
dim(chr20_melt)
head(chr20_melt)
chr20_melt <- chr20_melt[complete.cases(chr20_melt),]
dim(chr20_melt)
head(chr20_melt)

chr20_melt$Var1.position <- probeDetails$position[match(chr20_melt$Var1, probeDetails$name)]
chr20_melt$Var2.position <- probeDetails$position[match(chr20_melt$Var2, probeDetails$name)]
chr20_melt$absdistance <- abs(chr20_melt$Var1.position - chr20_melt$Var2.position)
chr20_melt$cis_or_trans <- as.character("trans")
chr20_melt$cis_or_trans[chr20_melt$absdistance < 1000000] <- "cis"

chr20melt_10kb <- chr20_melt[chr20_melt$absdistance <= 10000,]
dim(chr20melt_10kb)
head(chr20melt_10kb)
chr20_10k_sort <- chr20melt_10kb[order(chr20melt_10kb$absdistance),]
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
head(chr20_10k_pos)
# calculate the mean correlation in each bin
chr20_10k_pos_valuemean <- aggregate(chr20_10k_pos[, 3], list(chr20_10k_pos$group), mean) # [,3] is value
class(chr20_10k_pos_valuemean)
head(chr20_10k_pos_valuemean)
names(chr20_10k_pos_valuemean) <- c("Group.1", "value")
head(chr20_10k_pos_valuemean)
# calculate the median distance of each bin
chr20_10k_pos_mediandistance <- aggregate(chr20_10k_pos[, 6], list(chr20_10k_pos$group), median) # [,4] is absdistance
class(chr20_10k_pos_mediandistance)
head(chr20_10k_pos_mediandistance)
names(chr20_10k_pos_mediandistance) <- c("Group.1", "absdistance")
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
names(chr20_10k_neg_valuemean) <- c("Group.1", "value")
head(chr20_10k_neg_valuemean)
# calculate the median distance of each bin
chr20_10k_neg_mediandistance <- aggregate(chr20_10k_neg[, 6], list(chr20_10k_neg$group), median)
class(chr20_10k_neg_mediandistance)
head(chr20_10k_neg_mediandistance)
names(chr20_10k_neg_mediandistance) <- c("Group.1", "absdistance")
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

jpeg(filename = "/path/to/chr20_decayplot_F7_cisadjusted_10k_SD.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot()+
  geom_errorbar(data=chr20_10k_pos_df,aes(x=absdistance,ymin=value-x, ymax=value+x), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr20_10k_pos_df,aes(x=absdistance,y=value,colour="Positive"),size=0.3,alpha=1)+
  geom_errorbar(data=chr20_10k_neg_df,aes(x=absdistance,ymin=value-x, ymax=value+x), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr20_10k_neg_df,aes(x=absdistance,y=value,colour="Negative"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Positive", "Negative"),
                      values = c("Positive"="#440154FF", "Negative"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 20: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds, with\nstrongest cis mQTL adjusted for"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

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

jpeg(filename = "/path/to/chr20_histogram_F7_1k_cissnpadjusted.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(cor.group)+
  geom_bar(aes(x=Var1, y=Freq),stat="identity", position=position_dodge(),fill="#FDE725FF") +
  geom_text(aes(x=Var1, y=Freq,label=Freq), position=position_dodge(width=0), vjust=-0.1, hjust=0.5, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  labs(title=paste("Chromosome 20: values of cis correlations\nwithin 1kb in ARIES 7 year olds, adjusted\nfor cis SNPs"), x = "Correlation", y = "Frequency")
dev.off()

######

# looking at changes when we adjust for cis SNP / dont adjust for cis SNP

#####

# load and prepare the data that is unadjusted for cis SNPs: (code comes from the cis decay plots script)
load("/path/to/chr20_1.Rdata")
load("/path/to/chr20_2.Rdata")
load("/path/to/chr20_3.Rdata")

# unlist and combine these 3 files:
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

# reduce to 10kb:
chr20_10k_nosnp <- chr20[chr20$absdistance <= 10000]
dim(chr20_10k_nosnp)
head(chr20_10k_nosnp)
chr20_10k_sort_nosnp <- chr20_10k_nosnp[order(absdistance)]
head(chr20_10k_sort_nosnp)

# split to positive and negative correlations
chr20_10k_pos_nosnp <- chr20_10k_sort_nosnp[chr20_10k_sort_nosnp$value >= 0,]
dim(chr20_10k_pos_nosnp)
chr20_10k_neg_nosnp <- chr20_10k_sort_nosnp[chr20_10k_sort_nosnp$value < 0,]
dim(chr20_10k_neg_nosnp)

# bin the positive correlations on distance
chr20_10k_pos_nosnp$group <- as.numeric(cut2(chr20_10k_pos_nosnp$absdistance, m=100))
summary(chr20_10k_pos_nosnp$group)
chr20_10k_pos_nosnp$group[1:900]
# calculate the mean correlation in each bin
chr20_10k_pos_nosnp_valuemean <- aggregate(chr20_10k_pos_nosnp[, 3], list(chr20_10k_pos_nosnp$group), mean)
class(chr20_10k_pos_nosnp_valuemean)
head(chr20_10k_pos_nosnp_valuemean)
# calculate the median distance of each bin
chr20_10k_pos_nosnp_mediandistance <- aggregate(chr20_10k_pos_nosnp[, 4], list(chr20_10k_pos_nosnp$group), median)
class(chr20_10k_pos_nosnp_mediandistance)
head(chr20_10k_pos_nosnp_mediandistance)
# merge together the mean correlation and median distance
chr20_10k_pos_nosnp_df <- merge(chr20_10k_pos_nosnp_valuemean, chr20_10k_pos_nosnp_mediandistance, by=c("Group.1"))
head(chr20_10k_pos_nosnp_df)

# bin the negative correlations on distance
chr20_10k_neg_nosnp$group <- as.numeric(cut2(chr20_10k_neg_nosnp$absdistance, m=100))
summary(chr20_10k_neg_nosnp$group)
chr20_10k_neg_nosnp$group[1:900]
# calculate the mean correlation in each bin
chr20_10k_neg_nosnp_valuemean <- aggregate(chr20_10k_neg_nosnp[, 3], list(chr20_10k_neg_nosnp$group), mean)
class(chr20_10k_neg_nosnp_valuemean)
head(chr20_10k_neg_nosnp_valuemean)
# calculate the median distance of each bin
chr20_10k_neg_nosnp_mediandistance <- aggregate(chr20_10k_neg_nosnp[, 4], list(chr20_10k_neg_nosnp$group), median)
class(chr20_10k_neg_nosnp_mediandistance)
head(chr20_10k_neg_nosnp_mediandistance)
# merge together the mean correlation and median distance
chr20_10k_neg_nosnp_df <- merge(chr20_10k_neg_nosnp_valuemean, chr20_10k_neg_nosnp_mediandistance, by=c("Group.1"))
head(chr20_10k_neg_nosnp_df)

# calculate the SD in each positive bin
chr20_10k_pos_nosnp_sd <- aggregate(chr20_10k_pos_nosnp[, 3], list(chr20_10k_pos_nosnp$group), sd)
class(chr20_10k_pos_nosnp_sd)
head(chr20_10k_pos_nosnp_sd)
chr20_10k_pos_nosnp_df <- merge(chr20_10k_pos_nosnp_df, chr20_10k_pos_nosnp_sd, by=c("Group.1"))
head(chr20_10k_pos_nosnp_df)
chr20_10k_pos_nosnp_df$n_in_group <- as.numeric(table(chr20_10k_pos_nosnp$group))
head(chr20_10k_pos_nosnp_df)
table(chr20_10k_pos_nosnp$group)

# calculate the SD in each negative bin
chr20_10k_neg_nosnp_sd <- aggregate(chr20_10k_neg_nosnp[, 3], list(chr20_10k_neg_nosnp$group), sd)
class(chr20_10k_neg_nosnp_sd)
head(chr20_10k_neg_nosnp_sd)
chr20_10k_neg_nosnp_df <- merge(chr20_10k_neg_nosnp_df, chr20_10k_neg_nosnp_sd, by=c("Group.1"))
head(chr20_10k_neg_nosnp_df)


# Now we're ready to compare
# the 10kb dataframe should have the same correlations in it
# we need to sort the Var 1 and Var2 columns, and then match them up in a new
# dataframe to compare the correlation. Then to calculate the percentage change (ie
# correlation of 0.8 down to 0.4 is a 50% change) we calculate a column with 
# (adjusted_cor-unadjusted_cor)/unadjusted_cor * 100


# nope, different numbers of positive and negative now so there are a different number
# of bins. So we need to go back to the original df:

#chr20_10k_pos_nosnp
head(chr20_10k_pos_nosnp)
chr20_10k_pos_nosnp <- as.data.frame(chr20_10k_pos_nosnp)
head(chr20_10k_pos_nosnp)
#chr20_10k_sort
head(chr20_10k_sort)

# sort Var1 and Var2 in both df, so that var1 and var2 are the same way round for 
# both dfs:
chr20_10k_pos_nosnp$Var1 <- as.character(chr20_10k_pos_nosnp$Var1)
chr20_10k_pos_nosnp$Var2 <- as.character(chr20_10k_pos_nosnp$Var2)
chr20_10k_sort$Var1 <- as.character(chr20_10k_sort$Var1)
chr20_10k_sort$Var2 <- as.character(chr20_10k_sort$Var2)
chr20_10k_pos_nosnp_sort <- as.data.frame(t(apply(chr20_10k_pos_nosnp[1:2], 1, sort)))
chr20_10k_sort_sort <- as.data.frame(t(apply(chr20_10k_sort[1:2], 1, sort)))
head(chr20_10k_pos_nosnp_sort)
dim(chr20_10k_pos_nosnp_sort)
head(chr20_10k_sort_sort)
dim(chr20_10k_sort_sort)

# add the correlations back on:
# not adjusted for cis SNP
chr20_10k_pos_nosnp_sort$value <- chr20_10k_pos_nosnp$value
chr20_10k_pos_nosnp_sort$absdistance <- chr20_10k_pos_nosnp$absdistance
chr20_10k_pos_nosnp_sort$group <- chr20_10k_pos_nosnp$group
head(chr20_10k_pos_nosnp_sort)
dim(chr20_10k_pos_nosnp_sort)
# adjusted for cis SNP
chr20_10k_sort_sort$value <- chr20_10k_sort$value
chr20_10k_sort_sort$absdistance <- chr20_10k_sort$absdistance
head(chr20_10k_sort_sort)
dim(chr20_10k_sort_sort)

# now merge these DFs:
chr20_10k_merge <- merge(chr20_10k_pos_nosnp_sort, chr20_10k_sort_sort, by=c("V1", "V2"))
head(chr20_10k_merge)
dim(chr20_10k_merge)
chr20_10k_merge <- chr20_10k_merge[order(chr20_10k_merge$absdistance.x)]
#  V1         V2   value.x absdistance.x group     value.y     absdistance.y

# calc median distance and the 2 mean correlations for the df:
# calculate the mean unadjusted correlation in each bin
chr20_10k_pos_nosnp_mean <- aggregate(chr20_10k_merge[, 3], list(chr20_10k_merge$group), mean)
class(chr20_10k_pos_nosnp_mean)
head(chr20_10k_pos_nosnp_mean)
# calculate the mean cis SNP adjusted correlation in each bin
chr20_10k_pos_mean <- aggregate(chr20_10k_merge[, 6], list(chr20_10k_merge$group), mean)
class(chr20_10k_pos_mean)
head(chr20_10k_pos_mean)
# calculate the median distance of each bin
chr20_10k_pos_meddistance <- aggregate(chr20_10k_merge[, 4], list(chr20_10k_merge$group), median)
class(chr20_10k_pos_meddistance)
head(chr20_10k_pos_meddistance)
# merge together the mean correlation and median distance
chr20_10k_merge_df <- merge(chr20_10k_pos_nosnp_mean, chr20_10k_pos_meddistance, by=c("Group.1"))
head(chr20_10k_merge_df)
chr20_10k_merge_df <- merge(chr20_10k_merge_df, chr20_10k_pos_mean, by=c("Group.1"))
head(chr20_10k_merge_df)

jpeg(filename = "/path/to/chr20_decayplot_F7_adjustedvsnonadjusted_cissnp.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot()+
  #geom_errorbar(data=chr20_10k_merge_df,aes(x=absdistance,ymin=value-x, ymax=value+x), width=5, alpha=0.2, size=0.1,colour="#39568CFF") +
  geom_line(data=chr20_10k_merge_df,aes(x=x.y,y=x.x,colour="Unadjusted"),size=0.3,alpha=1)+
  #geom_errorbar(data=chr20_10k_neg_df,aes(x=absdistance,ymin=value-x, ymax=value+x), width=0.7, alpha=0.1, colour="#20A387FF") +
  geom_line(data=chr20_10k_merge_df,aes(x=x.y,y=x,colour="Cis SNP adjusted"),size=0.3,alpha=1)+
  scale_colour_manual("", 
                      breaks = c("Unadjusted", "Cis SNP adjusted"),
                      values = c("Unadjusted"="#440154FF", "Cis SNP adjusted"="#238A8DFF")) +
  labs(x="Distance",y="correlation")+
  theme_bw()+
  theme(aspect.ratio=1)+
  labs(title=paste("Chromosome 20: decay plot of pairwise correlations\nvs genomic distance in ARIES 7 year olds, comparison\nsof adjusting for cis mQTL vs not"), x = "Genomic distance (bp)", y = "Pairwise correlation")
dev.off()

# look at the change in correlation after adjustment:
chr20_10k_merge$mean_change <- chr20_10k_merge$value.y-chr20_10k_merge$value.x
# calculate the mean change in correlation in each bin
chr20_10k_meanchange <- aggregate(chr20_10k_merge[, 8], list(chr20_10k_merge$group), mean)
class(chr20_10k_meanchange)
head(chr20_10k_meanchange)
# calculate the median distance of each bin
chr20_10k_meanchangedistance <- aggregate(chr20_10k_merge[, 4], list(chr20_10k_merge$group), median)
class(chr20_10k_meanchangedistance)
head(chr20_10k_meanchangedistance)
# merge together the mean correlation and median distance
chr20_10k_change <- merge(chr20_10k_meanchange, chr20_10k_meanchangedistance, by=c("Group.1"))
head(chr20_10k_change)
chr20_10k_change$abschange <- abs(chr20_10k_change$x.x)

jpeg(filename = "/path/to/chr20_decayplot_F7_meanpercentchange.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot()+
  geom_line(data=chr20_10k_change,aes(x=x.y,y=x.x),colour="#440154FF",size=0.3,alpha=1)+
  labs(x="Distance",y="Mean change in correlation (r)")+
  theme_bw()+
  theme(aspect.ratio=1)+
  scale_y_reverse(limits=c(0,-0.2))+
  labs(title=paste("Chromosome 20: mean change in pairwise correlation R2\nvs genomic distance in ARIES 7 year olds, when\n cis mQTLs are adjusted for"), x = "Genomic distance (bp)", y = "Mean change in correlation R2")
dev.off()


# Finding the CpGs most affected by cis SNP

head(chr20_10k_merge)
dim(chr20_10k_merge)
summary(chr20_10k_merge$mean_change)

jpeg(filename = "/path/to/chr20_meanchange_hist.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot() + 
  geom_histogram(data=chr20_10k_merge,aes(x=mean_change),color="black", fill="#238A8DFF")+
  # geom_text(aes(label=Number), position=position_dodge(width=0), vjust=0.5, hjust=1.02, size=3)+
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  labs(title=paste("Chromosome 20: histogram of mean change in pairwise correlations,\nwhen best cis SNP is adjusted for, in ARIES 7 year olds"), x = "Mean change in correlation with cis SNP adjustment", y = "Frequency")
dev.off()

# let's subset to correlations that change by over 0.5:
biggestchange <- chr20_10k_merge[chr20_10k_merge$mean_change < -0.5,]
dim(biggestchange)
# there are 687 of them.
head(biggestchange)
# let's rename some of the variables for clarity:
#  V1         V2   value.x absdistance.x group     value.y     absdistance.y   mean_change
names(biggestchange) <- c("Var1","Var2","unadjusted_cor","absdistance","bin","adjusted_cor","absdistance2","mean_change")
summary(biggestchange$unadjusted_cor)
# Min.    1st Qu. Median    Mean  3rd Qu.  Max. 
# 0.4294  0.5558  0.6215  0.6336  0.7013  0.9501 
summary(biggestchange$adjusted_cor)
# so these are moderate to high correlations that are being reduced to:
#  Min.       1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.084266 -0.004364  0.025634  0.036202  0.068226  0.238784 

# load best cis SNP from F7_chr22_godmc_BC3 script
load("/path/to/bestcisSNP_godmc.Rdata")
head(temp)
biggestchange$Var1.snp <- temp$cissnp[match(biggestchange$Var1, temp$cpg)]
biggestchange$Var2.snp <- temp$cissnp[match(biggestchange$Var2, temp$cpg)]
head(biggestchange)

biggestchange_snps.1 <- as.character(biggestchange$Var1.snp)
length(biggestchange_snps.1)
biggestchange_snps.2 <- as.character(biggestchange$Var2.snp)
length(biggestchange_snps.2)
biggestchange_snps <- c(biggestchange_snps.1,biggestchange_snps.2)
length(biggestchange_snps)
biggestchange_snps <- unique(biggestchange_snps)
length(biggestchange_snps)
write.table(biggestchange_snps, file="/path/to/biggestchange_snps.txt", row.names = F, col.names = F, quote = F)