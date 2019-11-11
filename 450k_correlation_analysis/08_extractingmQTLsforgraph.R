################# 2. mQTLS ###########################

library(reshape2)
library(ggplot2)
library(viridis)
######################################################

# analysis looking at how many mQTLs (0,1 or 2) are in the correlations at each
# correlation band.

####### F7 ###########

load("/path/to/F7_mQTL_results_percorrelation.Rdata")
F7_mqtls_melt <- melt(F7_mqtls_correlations)

# split cis and trans correlations into 2 dataframes:
F7_mqtls_melt_cis <- F7_mqtls_melt[F7_mqtls_melt$Var2=="cis",]
F7_mqtls_melt_trans <- F7_mqtls_melt[F7_mqtls_melt$Var2=="trans",]


#### CIS 

# merge the outputs where there are multiple files per correlation band:
y <- names(F7_mqtls_melt)
# minus point 3
minuspoint3names <- as.character(grep("minuspoint3", F7_mqtls_melt_cis$L1, value=TRUE))
minuspoint3names <- unique(minuspoint3names)
minuspoint3_0 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  minuspoint3_0 <- minuspoint3_0 + x
}
minuspoint3_0 <- data.frame("0", "cis", "Freq", minuspoint3_0, "minuspoint3")
names(minuspoint3_0) <- y
minuspoint3_1 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  minuspoint3_1 <- minuspoint3_1 + x
}
minuspoint3_1 <- data.frame("1", "cis", "Freq", minuspoint3_1, "minuspoint3")
names(minuspoint3_1) <- y
minuspoint3_2 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  minuspoint3_2 <- minuspoint3_2 + x
}
minuspoint3_2 <- data.frame("2", "cis", "Freq", minuspoint3_2, "minuspoint3")
names(minuspoint3_2) <- y
minuspoint3 <- rbind.data.frame(minuspoint3_0, minuspoint3_1, minuspoint3_2)

# minus point 2
minuspoint2names <- as.character(grep("minuspoint2", F7_mqtls_melt_cis$L1, value=TRUE))
minuspoint2names <- unique(minuspoint2names)
minuspoint2_0 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  minuspoint2_0 <- minuspoint2_0 + x
}
minuspoint2_0 <- data.frame("0", "cis", "Freq", minuspoint2_0, "minuspoint2")
names(minuspoint2_0) <- y
minuspoint2_1 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  minuspoint2_1 <- minuspoint2_1 + x
}
minuspoint2_1 <- data.frame("1", "cis", "Freq", minuspoint2_1, "minuspoint2")
names(minuspoint2_1) <- y
minuspoint2_2 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  minuspoint2_2 <- minuspoint2_2 + x
}
minuspoint2_2 <- data.frame("2", "cis", "Freq", minuspoint2_2, "minuspoint2")
names(minuspoint2_2) <- y
minuspoint2 <- rbind.data.frame(minuspoint2_0, minuspoint2_1, minuspoint2_2)

# minus point 1
minuspoint1names <- as.character(grep("minuspoint1", F7_mqtls_melt_cis$L1, value=TRUE))
minuspoint1names <- unique(minuspoint1names)
minuspoint1_0 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  minuspoint1_0 <- minuspoint1_0 + x
}
minuspoint1_0 <- data.frame("0", "cis", "Freq", minuspoint1_0, "minuspoint1")
names(minuspoint1_0) <- y
minuspoint1_1 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  minuspoint1_1 <- minuspoint1_1 + x
}
minuspoint1_1 <- data.frame("1", "cis", "Freq", minuspoint1_1, "minuspoint1")
names(minuspoint1_1) <- y
minuspoint1_2 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  minuspoint1_2 <- minuspoint1_2 + x
}
minuspoint1_2 <- data.frame("2", "cis", "Freq", minuspoint1_2, "minuspoint1")
names(minuspoint1_2) <- y
minuspoint1 <- rbind.data.frame(minuspoint1_0, minuspoint1_1, minuspoint1_2)

# zero
zeronames <- as.character(grep("zero", F7_mqtls_melt_cis$L1, value=TRUE))
zeronames <- unique(zeronames)
zero_0 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  zero_0 <- zero_0 + x
}
zero_0 <- data.frame("0", "cis", "Freq", zero_0, "zero")
names(zero_0) <- y
zero_1 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  zero_1 <- zero_1 + x
}
zero_1 <- data.frame("1", "cis", "Freq", zero_1, "zero")
names(zero_1) <- y
zero_2 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  zero_2 <- zero_2 + x
}
zero_2 <- data.frame("2", "cis", "Freq", zero_2, "zero")
names(zero_2) <- y
zero <- rbind.data.frame(zero_0, zero_1, zero_2)

# point 1
point1names <- as.character(grep("^point1", F7_mqtls_melt_cis$L1, value=TRUE))
point1names <- unique(point1names)
point1_0 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  point1_0 <- point1_0 + x
}
point1_0 <- data.frame("0", "cis", "Freq", point1_0, "point1")
names(point1_0) <- y
point1_1 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  point1_1 <- point1_1 + x
}
point1_1 <- data.frame("1", "cis", "Freq", point1_1, "point1")
names(point1_1) <- y
point1_2 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  point1_2 <- point1_2 + x
}
point1_2 <- data.frame("2", "cis", "Freq", point1_2, "point1")
names(point1_2) <- y
point1 <- rbind.data.frame(point1_0, point1_1, point1_2)

# point 2
point2names <- as.character(grep("^point2", F7_mqtls_melt_cis$L1, value=TRUE))
point2names <- unique(point2names)
point2_0 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  point2_0 <- point2_0 + x
}
point2_0 <- data.frame("0", "cis", "Freq", point2_0, "point2")
names(point2_0) <- y
point2_1 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  point2_1 <- point2_1 + x
}
point2_1 <- data.frame("1", "cis", "Freq", point2_1, "point2")
names(point2_1) <- y
point2_2 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  point2_2 <- point2_2 + x
}
point2_2 <- data.frame("2", "cis", "Freq", point2_2, "point2")
names(point2_2) <- y
point2 <- rbind.data.frame(point2_0, point2_1, point2_2)

# point 3
point3names <- as.character(grep("^point3", F7_mqtls_melt_cis$L1, value=TRUE))
point3names <- unique(point3names)
point3_0 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="0"])
  point3_0 <- point3_0 + x
}
point3_0 <- data.frame("0", "cis", "Freq", point3_0, "point3")
names(point3_0) <- y
point3_1 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="1"])
  point3_1 <- point3_1 + x
}
point3_1 <- data.frame("1", "cis", "Freq", point3_1, "point3")
names(point3_1) <- y
point3_2 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_cis$value[F7_mqtls_melt_cis$L1==i & F7_mqtls_melt_cis$Var1=="2"])
  point3_2 <- point3_2 + x
}
point3_2 <- data.frame("2", "cis", "Freq", point3_2, "point3")
names(point3_2) <- y
point3 <- rbind.data.frame(point3_0, point3_1, point3_2)

# list all the complete ones together:
F7_mqtls_melt_cis_list <- list()
F7_mqtls_melt_cis_list$minus1 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minus1.Rdata",]
F7_mqtls_melt_cis_list$minuspoint9 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint9.Rdata",]
F7_mqtls_melt_cis_list$minuspoint8 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint8.Rdata",]
F7_mqtls_melt_cis_list$minuspoint7 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint7.Rdata",]
F7_mqtls_melt_cis_list$minuspoint6 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint6.Rdata",]
F7_mqtls_melt_cis_list$minuspoint5 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint5.Rdata",]
F7_mqtls_melt_cis_list$minuspoint4 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="minuspoint4.Rdata",]
F7_mqtls_melt_cis_list$minuspoint3 <- minuspoint3
F7_mqtls_melt_cis_list$minuspoint2 <- minuspoint2
F7_mqtls_melt_cis_list$minuspoint1 <- minuspoint1
F7_mqtls_melt_cis_list$zero <- zero
F7_mqtls_melt_cis_list$point1 <- point1
F7_mqtls_melt_cis_list$point2 <- point2
F7_mqtls_melt_cis_list$point3 <- point3
F7_mqtls_melt_cis_list$point4 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point4.Rdata",]
F7_mqtls_melt_cis_list$point5 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point5.Rdata",]
F7_mqtls_melt_cis_list$point6 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point6.Rdata",]
F7_mqtls_melt_cis_list$point7 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point7.Rdata",]
F7_mqtls_melt_cis_list$point8 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point8.Rdata",]
F7_mqtls_melt_cis_list$point9 <- F7_mqtls_melt_cis[F7_mqtls_melt_cis$L1=="point9.Rdata",]

F7_mqtls <- do.call(rbind.data.frame, F7_mqtls_melt_cis_list)

class(F7_mqtls$L1)
F7_mqtls$L1[F7_mqtls$L1 == "minus1.Rdata"] <- c("-1 to -0.9")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint1"] <- c("-0.1 to 0")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint2"] <- c("-0.2 to -0.1")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint3"] <- c("-0.3 to -0.2")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint4.Rdata"] <- c("-0.4 to -0.3")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint5.Rdata"] <- c("-0.5 to -0.4")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint6.Rdata"] <- c("-0.6 to -0.5")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint7.Rdata"] <- c("-0.7 to -0.6")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint8.Rdata"] <- c("-0.8 to -0.7")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint9.Rdata"] <- c("-0.9 to -0.8")
F7_mqtls$L1[F7_mqtls$L1 == "point1"] <- c("0.1 to 0.2")
F7_mqtls$L1[F7_mqtls$L1 == "point2"] <- c("0.2 to 0.3")
F7_mqtls$L1[F7_mqtls$L1 == "point3"] <- c("0.3 to 0.4")
F7_mqtls$L1[F7_mqtls$L1 == "point4.Rdata"] <- c("0.4 to 0.5")
F7_mqtls$L1[F7_mqtls$L1 == "point5.Rdata"] <- c("0.5 to 0.6")
F7_mqtls$L1[F7_mqtls$L1 == "point6.Rdata"] <- c("0.6 to 0.7")
F7_mqtls$L1[F7_mqtls$L1 == "point7.Rdata"] <- c("0.7 to 0.8")
F7_mqtls$L1[F7_mqtls$L1 == "point8.Rdata"] <- c("0.8 to 0.9")
F7_mqtls$L1[F7_mqtls$L1 == "point9.Rdata"] <- c("0.9 to 1")
F7_mqtls$L1[F7_mqtls$L1 == "zero"] <- c("0 to 0.1")

# there are no pairs with 0 or 1 mQTLs in F7 in the -1 to -0.9 band so we add a row
# to show that:
F7_mqtls[59,] <- c(0,"cis","Freq",0,"-1 to -0.9")
F7_mqtls[60,] <- c(1,"cis","Freq",0,"-1 to -0.9")

F7_mqtls$value <- as.numeric(F7_mqtls$value)
F7_mqtls$Var1 <- ordered(F7_mqtls$Var1,
                              levels=c("0","1","2"))

totalNcorrelations <- as.data.frame(rowsum(F7_mqtls$value, group = F7_mqtls$L1))
totalNcorrelations$L1 <- rownames(totalNcorrelations)
class(F7_mqtls$L1)
F7_mqtls$totalNcorrelations <- totalNcorrelations$V1[match(F7_mqtls$L1, totalNcorrelations$L1)]
F7_mqtls$percentperband <- (F7_mqtls$value/F7_mqtls$totalNcorrelations)*100
F7_mqtls$L1 <- as.factor(F7_mqtls$L1)
F7_mqtls$L1 <- ordered(F7_mqtls$L1, 
                            levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                     "-0.8 to -0.7", "-0.7 to -0.6",
                                     "-0.6 to -0.5", "-0.5 to -0.4",
                                     "-0.4 to -0.3", "-0.3 to -0.2",
                                     "-0.2 to -0.1", "-0.1 to 0",
                                     "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                     "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                     "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                     "0.9 to 1"))


jpeg(filename ="/path/to/F7_mQTLpercentage_percorrelation.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(F7_mqtls, aes(x=L1, y=percentperband, fill=Var1)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_viridis(name = "Number of\nmQTLs per\ncorrelation", discrete = T, begin=0.2, end = 0.85)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=12))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0))+
  ylim(0,100)+
  labs(title="Bar plot showing whether neither, one or both probes in\n a cis correlating pair have mQTLs, at 7 years old",x="Correlation band, from -1 to 1, in increments of 0.1", y="Percentage of correlations with an mQTL")
dev.off()


### TRANS

# merge the outputs where there are multiple files per correlation band:
y <- names(F7_mqtls_melt)
# minus point 3
minuspoint3names <- as.character(grep("minuspoint3", F7_mqtls_melt_trans$L1, value=TRUE))
minuspoint3names <- unique(minuspoint3names)
minuspoint3_0 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  minuspoint3_0 <- minuspoint3_0 + x
}
minuspoint3_0 <- data.frame("0", "trans", "Freq", minuspoint3_0, "minuspoint3")
names(minuspoint3_0) <- y
minuspoint3_1 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  minuspoint3_1 <- minuspoint3_1 + x
}
minuspoint3_1 <- data.frame("1", "trans", "Freq", minuspoint3_1, "minuspoint3")
names(minuspoint3_1) <- y
minuspoint3_2 <- 0
for (i in minuspoint3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  minuspoint3_2 <- minuspoint3_2 + x
}
minuspoint3_2 <- data.frame("2", "trans", "Freq", minuspoint3_2, "minuspoint3")
names(minuspoint3_2) <- y
minuspoint3 <- rbind.data.frame(minuspoint3_0, minuspoint3_1, minuspoint3_2)

# minus point 2
minuspoint2names <- as.character(grep("minuspoint2", F7_mqtls_melt_trans$L1, value=TRUE))
minuspoint2names <- unique(minuspoint2names)
minuspoint2_0 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  minuspoint2_0 <- minuspoint2_0 + x
}
minuspoint2_0 <- data.frame("0", "trans", "Freq", minuspoint2_0, "minuspoint2")
names(minuspoint2_0) <- y
minuspoint2_1 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  minuspoint2_1 <- minuspoint2_1 + x
}
minuspoint2_1 <- data.frame("1", "trans", "Freq", minuspoint2_1, "minuspoint2")
names(minuspoint2_1) <- y
minuspoint2_2 <- 0
for (i in minuspoint2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  minuspoint2_2 <- minuspoint2_2 + x
}
minuspoint2_2 <- data.frame("2", "trans", "Freq", minuspoint2_2, "minuspoint2")
names(minuspoint2_2) <- y
minuspoint2 <- rbind.data.frame(minuspoint2_0, minuspoint2_1, minuspoint2_2)

# minus point 1
minuspoint1names <- as.character(grep("minuspoint1", F7_mqtls_melt_trans$L1, value=TRUE))
minuspoint1names <- unique(minuspoint1names)
minuspoint1_0 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  minuspoint1_0 <- minuspoint1_0 + x
}
minuspoint1_0 <- data.frame("0", "trans", "Freq", minuspoint1_0, "minuspoint1")
names(minuspoint1_0) <- y
minuspoint1_1 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  minuspoint1_1 <- minuspoint1_1 + x
}
minuspoint1_1 <- data.frame("1", "trans", "Freq", minuspoint1_1, "minuspoint1")
names(minuspoint1_1) <- y
minuspoint1_2 <- 0
for (i in minuspoint1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  minuspoint1_2 <- minuspoint1_2 + x
}
minuspoint1_2 <- data.frame("2", "trans", "Freq", minuspoint1_2, "minuspoint1")
names(minuspoint1_2) <- y
minuspoint1 <- rbind.data.frame(minuspoint1_0, minuspoint1_1, minuspoint1_2)

# zero
zeronames <- as.character(grep("zero", F7_mqtls_melt_trans$L1, value=TRUE))
zeronames <- unique(zeronames)
zero_0 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  zero_0 <- zero_0 + x
}
zero_0 <- data.frame("0", "trans", "Freq", zero_0, "zero")
names(zero_0) <- y
zero_1 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  zero_1 <- zero_1 + x
}
zero_1 <- data.frame("1", "trans", "Freq", zero_1, "zero")
names(zero_1) <- y
zero_2 <- 0
for (i in zeronames){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  zero_2 <- zero_2 + x
}
zero_2 <- data.frame("2", "trans", "Freq", zero_2, "zero")
names(zero_2) <- y
zero <- rbind.data.frame(zero_0, zero_1, zero_2)

# point 1
point1names <- as.character(grep("^point1", F7_mqtls_melt_trans$L1, value=TRUE))
point1names <- unique(point1names)
point1_0 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  point1_0 <- point1_0 + x
}
point1_0 <- data.frame("0", "trans", "Freq", point1_0, "point1")
names(point1_0) <- y
point1_1 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  point1_1 <- point1_1 + x
}
point1_1 <- data.frame("1", "trans", "Freq", point1_1, "point1")
names(point1_1) <- y
point1_2 <- 0
for (i in point1names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  point1_2 <- point1_2 + x
}
point1_2 <- data.frame("2", "trans", "Freq", point1_2, "point1")
names(point1_2) <- y
point1 <- rbind.data.frame(point1_0, point1_1, point1_2)

# point 2
point2names <- as.character(grep("^point2", F7_mqtls_melt_trans$L1, value=TRUE))
point2names <- unique(point2names)
point2_0 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  point2_0 <- point2_0 + x
}
point2_0 <- data.frame("0", "trans", "Freq", point2_0, "point2")
names(point2_0) <- y
point2_1 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  point2_1 <- point2_1 + x
}
point2_1 <- data.frame("1", "trans", "Freq", point2_1, "point2")
names(point2_1) <- y
point2_2 <- 0
for (i in point2names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  point2_2 <- point2_2 + x
}
point2_2 <- data.frame("2", "trans", "Freq", point2_2, "point2")
names(point2_2) <- y
point2 <- rbind.data.frame(point2_0, point2_1, point2_2)

# point 3
point3names <- as.character(grep("^point3", F7_mqtls_melt_trans$L1, value=TRUE))
point3names <- unique(point3names)
point3_0 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="0"])
  point3_0 <- point3_0 + x
}
point3_0 <- data.frame("0", "trans", "Freq", point3_0, "point3")
names(point3_0) <- y
point3_1 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="1"])
  point3_1 <- point3_1 + x
}
point3_1 <- data.frame("1", "trans", "Freq", point3_1, "point3")
names(point3_1) <- y
point3_2 <- 0
for (i in point3names){
  x <- as.numeric(F7_mqtls_melt_trans$value[F7_mqtls_melt_trans$L1==i & F7_mqtls_melt_trans$Var1=="2"])
  point3_2 <- point3_2 + x
}
point3_2 <- data.frame("2", "trans", "Freq", point3_2, "point3")
names(point3_2) <- y
point3 <- rbind.data.frame(point3_0, point3_1, point3_2)

# list all the complete ones together:
F7_mqtls_melt_trans_list <- list()
F7_mqtls_melt_trans_list$minus1 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minus1.Rdata",]
F7_mqtls_melt_trans_list$minuspoint9 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint9.Rdata",]
F7_mqtls_melt_trans_list$minuspoint8 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint8.Rdata",]
F7_mqtls_melt_trans_list$minuspoint7 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint7.Rdata",]
F7_mqtls_melt_trans_list$minuspoint6 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint6.Rdata",]
F7_mqtls_melt_trans_list$minuspoint5 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint5.Rdata",]
F7_mqtls_melt_trans_list$minuspoint4 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="minuspoint4.Rdata",]
F7_mqtls_melt_trans_list$minuspoint3 <- minuspoint3
F7_mqtls_melt_trans_list$minuspoint2 <- minuspoint2
F7_mqtls_melt_trans_list$minuspoint1 <- minuspoint1
F7_mqtls_melt_trans_list$zero <- zero
F7_mqtls_melt_trans_list$point1 <- point1
F7_mqtls_melt_trans_list$point2 <- point2
F7_mqtls_melt_trans_list$point3 <- point3
F7_mqtls_melt_trans_list$point4 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point4.Rdata",]
F7_mqtls_melt_trans_list$point5 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point5.Rdata",]
F7_mqtls_melt_trans_list$point6 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point6.Rdata",]
F7_mqtls_melt_trans_list$point7 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point7.Rdata",]
F7_mqtls_melt_trans_list$point8 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point8.Rdata",]
F7_mqtls_melt_trans_list$point9 <- F7_mqtls_melt_trans[F7_mqtls_melt_trans$L1=="point9.Rdata",]

F7_mqtls <- do.call(rbind.data.frame, F7_mqtls_melt_trans_list)

# there are no trans pairs with 0 or 1 or 2 mQTLs in F7 in the -1 to -0.9 band so we add a row
# to show that:
F7_mqtls[58,] <- c(0,"trans","Freq",0,"-1 to -0.9")
F7_mqtls[59,] <- c(1,"trans","Freq",0,"-1 to -0.9")
F7_mqtls[60,] <- c(2,"trans","Freq",0,"-1 to -0.9")

class(F7_mqtls$L1)
F7_mqtls$L1[F7_mqtls$L1 == "minus1.Rdata"] <- c("-1 to -0.9")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint1"] <- c("-0.1 to 0")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint2"] <- c("-0.2 to -0.1")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint3"] <- c("-0.3 to -0.2")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint4.Rdata"] <- c("-0.4 to -0.3")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint5.Rdata"] <- c("-0.5 to -0.4")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint6.Rdata"] <- c("-0.6 to -0.5")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint7.Rdata"] <- c("-0.7 to -0.6")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint8.Rdata"] <- c("-0.8 to -0.7")
F7_mqtls$L1[F7_mqtls$L1 == "minuspoint9.Rdata"] <- c("-0.9 to -0.8")
F7_mqtls$L1[F7_mqtls$L1 == "point1"] <- c("0.1 to 0.2")
F7_mqtls$L1[F7_mqtls$L1 == "point2"] <- c("0.2 to 0.3")
F7_mqtls$L1[F7_mqtls$L1 == "point3"] <- c("0.3 to 0.4")
F7_mqtls$L1[F7_mqtls$L1 == "point4.Rdata"] <- c("0.4 to 0.5")
F7_mqtls$L1[F7_mqtls$L1 == "point5.Rdata"] <- c("0.5 to 0.6")
F7_mqtls$L1[F7_mqtls$L1 == "point6.Rdata"] <- c("0.6 to 0.7")
F7_mqtls$L1[F7_mqtls$L1 == "point7.Rdata"] <- c("0.7 to 0.8")
F7_mqtls$L1[F7_mqtls$L1 == "point8.Rdata"] <- c("0.8 to 0.9")
F7_mqtls$L1[F7_mqtls$L1 == "point9.Rdata"] <- c("0.9 to 1")
F7_mqtls$L1[F7_mqtls$L1 == "zero"] <- c("0 to 0.1")


F7_mqtls$value <- as.numeric(F7_mqtls$value)
F7_mqtls$Var1 <- ordered(F7_mqtls$Var1,
                         levels=c("0","1","2"))

totalNcorrelations <- as.data.frame(rowsum(F7_mqtls$value, group = F7_mqtls$L1))
totalNcorrelations$L1 <- rownames(totalNcorrelations)
class(F7_mqtls$L1)
F7_mqtls$totalNcorrelations <- totalNcorrelations$V1[match(F7_mqtls$L1, totalNcorrelations$L1)]
F7_mqtls$percentperband <- (F7_mqtls$value/F7_mqtls$totalNcorrelations)*100
F7_mqtls$L1 <- as.factor(F7_mqtls$L1)
F7_mqtls$L1 <- ordered(F7_mqtls$L1, 
                       levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                "-0.8 to -0.7", "-0.7 to -0.6",
                                "-0.6 to -0.5", "-0.5 to -0.4",
                                "-0.4 to -0.3", "-0.3 to -0.2",
                                "-0.2 to -0.1", "-0.1 to 0",
                                "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                "0.9 to 1"))


jpeg(filename ="/path/to/F7_mQTLpercentage_percorrelation_trans.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(F7_mqtls, aes(x=L1, y=percentperband, fill=Var1)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_viridis(name = "Number of\nmQTLs per\ncorrelation", discrete = T, begin=0.2, end = 0.85)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=12))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0))+
  ylim(0,100)+
  labs(title="Bar plot showing whether neither, one or both probes in\n a trans correlating pair have mQTLs, at 7 years old",x="Correlation band, from -1 to 1, in increments of 0.1", y="Percentage of correlations with an mQTL")
dev.off()
