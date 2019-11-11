# get the dimensions of the files

library(data.table)

# create list of 450k probes inc gene and position
probeDetails <- meffil.featureset("450k")
rownames(probeDetails) <- probeDetails$name
head(probeDetails)
dim(probeDetails)
probeDetails <- na.omit(probeDetails)
dim(probeDetails)

# list files of correlation bands:
files <- list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

get_dim <- function(files){
  output <- load(files)
  print(files)
  print(head(dat))
  dat <- as.data.frame(dat)
  dat <- dat[complete.cases(dat),]
  dat$Var1.position <- probeDetails$position[match(dat$Var1, rownames(probeDetails))]
  dat$Var2.position <- probeDetails$position[match(dat$Var2, rownames(probeDetails))]
  dat$distance <- (dat$Var1.position - dat$Var2.position)
  dat$Var1.chromosome <- probeDetails$chromosome[match(dat$Var1, rownames(probeDetails))]
  dat$Var2.chromosome <- probeDetails$chromosome[match(dat$Var2, rownames(probeDetails))]
  dat$absdistance <- abs(dat$distance)
  dat$cis_or_trans <- as.character(c("trans"))
  dat$cis_or_trans[dat$absdistance <= 1000000 & dat$Var1.chromosome == dat$Var2.chromosome] <- as.character(c("cis"))
  print(head(dat))
  print(dim(dat))
  resultslist <- table(dat$cis_or_trans)
  return(resultslist)
}
corbands_cistrans <- lapply(files, FUN=get_dim)

# substring length needs to be altered - this will want to be the file names
# at the end of 'files'
names(corbands_cistrans) <- print(substring(files,41))

summary(corbands_cistrans)

save(corbands_cistrans, file="/path/to/correl_bands_cis_trans_dimensions.Rdata")

#############################

# plot 

library(ggplot2)
library(viridis)
library(scales)
options(scipen=999)
library(reshape2)

load("/path/to/correl_bands_cis_trans_dimensions.Rdata")

# put all the dimensions from a list into a single df:
x <- do.call("rbind", corbands_cistrans)
x <- as.data.frame(x)
x$band <- as.character(rownames(x))
#y <- x[, -grep("2$", colnames(x))]
#y <- as.data.frame(t(y))
x$cis <- as.numeric(x$cis)
x$trans <- as.numeric(x$trans)

# where there are 0 correlations in a band, it inserts the other number (cis/trans) instead of 0:
# manually check and replace them with 0 here.

# identify which bands this is a problem for
test <- x[x$cis==x$trans,]

# set the appropriate ones to 0
x$trans[x$band=="minus1.Rdata"] <- 0

# sanity check - total cors should match expected number.
totalciscors <- sum(x$cis)
totaltranscors <- sum(x$trans)

### cis

minuspoint3names <- as.character(grep("minuspoint3", x$band, value=TRUE))
minuspoint3 <- 0
for (i in minuspoint3names){
  y <- as.numeric(x$cis[x$band == i])
  minuspoint3 <- minuspoint3 + y
}

minuspoint2names <- as.character(grep("minuspoint2", x$band, value=TRUE))
minuspoint2 <- 0
for (i in minuspoint2names){
  y <- as.numeric(x$cis[x$band == i])
  minuspoint2 <- minuspoint2 + y
}

minuspoint1names <- as.character(grep("minuspoint1", x$band, value=TRUE))
minuspoint1 <- 0
for (i in minuspoint1names){
  y <- as.numeric(x$cis[x$band == i])
  minuspoint1 <- minuspoint1 + y
}

zeronames <- as.character(grep("zero", x$band, value=TRUE))
zero <- 0
for (i in zeronames){
  y <- as.numeric(x$cis[x$band == i])
  zero <- zero + y
}

point1names <- as.character(grep("^point1", x$band, value=TRUE))
point1 <- 0
for (i in point1names){
  y <- as.numeric(x$cis[x$band == i])
  point1 <- point1 + y
}

point2names <- as.character(grep("^point2", x$band, value=TRUE))
point2 <- 0
for (i in point2names){
  y <- as.numeric(x$cis[x$band == i])
  point2 <- point2 + y
}

point3names <- as.character(grep("^point3", x$band, value=TRUE))
point3 <- 0
for (i in point3names){
  y <- as.numeric(x$cis[x$band == i])
  point3 <- point3 + y
}

cis_correlations <- list()
cis_correlations$minus1 <- x$cis[x$band=="minus1.Rdata"]
cis_correlations$minuspoint9 <- x$cis[x$band=="minuspoint9.Rdata"]
cis_correlations$minuspoint8 <- x$cis[x$band=="minuspoint8.Rdata"]
cis_correlations$minuspoint7 <- x$cis[x$band=="minuspoint7.Rdata"]
cis_correlations$minuspoint6 <- x$cis[x$band=="minuspoint6.Rdata"]
cis_correlations$minuspoint5 <- x$cis[x$band=="minuspoint5.Rdata"]
cis_correlations$minuspoint4 <- x$cis[x$band=="minuspoint4.Rdata"]
cis_correlations$minuspoint3 <- minuspoint3
cis_correlations$minuspoint2 <- minuspoint2
cis_correlations$minuspoint1 <- minuspoint1
cis_correlations$zero <- zero
cis_correlations$point1 <- point1
cis_correlations$point2 <- point2
cis_correlations$point3 <- point3
cis_correlations$point4 <- x$cis[x$band=="point4.Rdata"]
cis_correlations$point5 <- x$cis[x$band=="point5.Rdata"]
cis_correlations$point6 <- x$cis[x$band=="point6.Rdata"]
cis_correlations$point7 <- x$cis[x$band=="point7.Rdata"]
cis_correlations$point8 <- x$cis[x$band=="point8.Rdata"]
cis_correlations$point9 <- x$cis[x$band=="point9.Rdata"]

correlations <- as.data.frame(cis_correlations)
correlations <- as.data.frame(t(correlations))
names(correlations) <- c("number")
correlations$percentage <- (correlations$number/totalciscors)*100


## trans

minuspoint3names <- as.character(grep("minuspoint3", x$band, value=TRUE))
minuspoint3 <- 0
for (i in minuspoint3names){
  y <- as.numeric(x$trans[x$band == i])
  minuspoint3 <- minuspoint3 + y
}

minuspoint2names <- as.character(grep("minuspoint2", x$band, value=TRUE))
minuspoint2 <- 0
for (i in minuspoint2names){
  y <- as.numeric(x$trans[x$band == i])
  minuspoint2 <- minuspoint2 + y
}

minuspoint1names <- as.character(grep("minuspoint1", x$band, value=TRUE))
minuspoint1 <- 0
for (i in minuspoint1names){
  y <- as.numeric(x$trans[x$band == i])
  minuspoint1 <- minuspoint1 + y
}

zeronames <- as.character(grep("zero", x$band, value=TRUE))
zero <- 0
for (i in zeronames){
  y <- as.numeric(x$trans[x$band == i])
  zero <- zero + y
}

point1names <- as.character(grep("^point1", x$band, value=TRUE))
point1 <- 0
for (i in point1names){
  y <- as.numeric(x$trans[x$band == i])
  point1 <- point1 + y
}

point2names <- as.character(grep("^point2", x$band, value=TRUE))
point2 <- 0
for (i in point2names){
  y <- as.numeric(x$trans[x$band == i])
  point2 <- point2 + y
}

point3names <- as.character(grep("^point3", x$band, value=TRUE))
point3 <- 0
for (i in point3names){
  y <- as.numeric(x$trans[x$band == i])
  point3 <- point3 + y
}

trans_correlations <- list()
trans_correlations$minus1 <- x$trans[x$band=="minus1.Rdata"]
trans_correlations$minuspoint9 <- x$trans[x$band=="minuspoint9.Rdata"]
trans_correlations$minuspoint8 <- x$trans[x$band=="minuspoint8.Rdata"]
trans_correlations$minuspoint7 <- x$trans[x$band=="minuspoint7.Rdata"]
trans_correlations$minuspoint6 <- x$trans[x$band=="minuspoint6.Rdata"]
trans_correlations$minuspoint5 <- x$trans[x$band=="minuspoint5.Rdata"]
trans_correlations$minuspoint4 <- x$trans[x$band=="minuspoint4.Rdata"]
trans_correlations$minuspoint3 <- minuspoint3
trans_correlations$minuspoint2 <- minuspoint2
trans_correlations$minuspoint1 <- minuspoint1
trans_correlations$zero <- zero
trans_correlations$point1 <- point1
trans_correlations$point2 <- point2
trans_correlations$point3 <- point3
trans_correlations$point4 <- x$trans[x$band=="point4.Rdata"]
trans_correlations$point5 <- x$trans[x$band=="point5.Rdata"]
trans_correlations$point6 <- x$trans[x$band=="point6.Rdata"]
trans_correlations$point7 <- x$trans[x$band=="point7.Rdata"]
trans_correlations$point8 <- x$trans[x$band=="point8.Rdata"]
trans_correlations$point9 <- x$trans[x$band=="point9.Rdata"]

transcorrelations <- as.data.frame(trans_correlations)
transcorrelations <- as.data.frame(t(transcorrelations))
names(transcorrelations) <- c("number")
transcorrelations$percentage <- (transcorrelations$number/totaltranscors)*100

correlations$cortype <- as.character("cis")
transcorrelations$cortype <- as.character("trans")

correlations$correlation.band <- as.character(rownames(correlations))
correlations$correlation.band <- as.character(c("-1 to -0.9", "-0.9 to -0.8", 
                                                "-0.8 to -0.7", "-0.7 to -0.6",
                                                "-0.6 to -0.5", "-0.5 to -0.4",
                                                "-0.4 to -0.3", "-0.3 to -0.2",
                                                "-0.2 to -0.1", "-0.1 to 0",
                                                "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                                "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                                "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                                "0.9 to 1"))
correlations$correlation.band <- factor(correlations$correlation.band, 
                                        levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                                 "-0.8 to -0.7", "-0.7 to -0.6",
                                                 "-0.6 to -0.5", "-0.5 to -0.4",
                                                 "-0.4 to -0.3", "-0.3 to -0.2",
                                                 "-0.2 to -0.1", "-0.1 to 0",
                                                 "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                                 "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                                 "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                                 "0.9 to 1"))

transcorrelations$correlation.band <- as.character(rownames(transcorrelations))
transcorrelations$correlation.band <- as.character(c("-1 to -0.9", "-0.9 to -0.8", 
                                                     "-0.8 to -0.7", "-0.7 to -0.6",
                                                     "-0.6 to -0.5", "-0.5 to -0.4",
                                                     "-0.4 to -0.3", "-0.3 to -0.2",
                                                     "-0.2 to -0.1", "-0.1 to 0",
                                                     "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                                     "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                                     "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                                     "0.9 to 1"))
transcorrelations$correlation.band <- factor(transcorrelations$correlation.band, 
                                             levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                                      "-0.8 to -0.7", "-0.7 to -0.6",
                                                      "-0.6 to -0.5", "-0.5 to -0.4",
                                                      "-0.4 to -0.3", "-0.3 to -0.2",
                                                      "-0.2 to -0.1", "-0.1 to 0",
                                                      "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                                      "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                                      "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                                      "0.9 to 1"))


allcorrelations <- rbind.data.frame(correlations, transcorrelations)


pdf(file ="/path/to/cistrans_proportions.pdf", width = 7, height = 5)
ggplot(allcorrelations, aes(x=correlation.band, y=percentage, fill=cortype)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_viridis(name = "Cis or\nTrans", discrete = T, begin=0.1, end = 0.5)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=12))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0))+
  #ylim(0,100)+
  labs(title="Bar plot showing the percentage of cis and trans\n correlations",x="Correlation band, from -1 to 1, in increments of 0.1", y="Percentage")
dev.off()

write.csv(correlations, file="/path/to/cis_percentages_table.csv")
write.csv(transcorrelations, file="/path/to/trans_percentages_table.csv")
