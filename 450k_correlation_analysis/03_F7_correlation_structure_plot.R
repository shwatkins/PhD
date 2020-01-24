# get the dimensions of the files

library(data.table)

files <- list.files(path="/path/to/correl_bands", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)
get_dim <- function(files){
  output <- load(files)
  print(files)
  print(head(dat))
  dat <- as.data.frame(dat)
  dat <- dat[complete.cases(dat),]
  print(head(dat))
  print(dim(dat))
  return(as.data.frame(dim(dat)))
}
corbands_dim <- lapply(files, FUN=get_dim)

names(corbands_dim) <- print(substring(files,41))

summary(corbands_dim)

save(corbands_dim, file="/path/to/correl_bands_dimensions.Rdata")

load("/path/to/correl_bands_dimensions.Rdata")

# put all the dimensions from a list into a single df:
x <- do.call("rbind", corbands_dim)
x <- as.data.frame(t(x))
y <- x[, -grep("2$", colnames(x))]
y <- as.data.frame(t(y))
totalcors <- sum(y$`dim(dat)`)
y <- as.data.frame(t(y))

# combine the dimensions where there was more than one file for a correlation 
# band (in this case between -0.3 and 0.3):

minuspoint3names <- as.character(grep("minuspoint3", colnames(y), value=TRUE))
minuspoint3 <- 0
for (i in minuspoint3names){
  x <- as.numeric(y[,i])
  minuspoint3 <- minuspoint3 + x
}

minuspoint2names <- as.character(grep("minuspoint2", colnames(y), value=TRUE))
minuspoint2 <- 0
for (i in minuspoint2names){
  x <- as.numeric(y[,i])
  minuspoint2 <- minuspoint2 + x
}

minuspoint1names <- as.character(grep("minuspoint1", colnames(y), value=TRUE))
minuspoint1 <- 0
for (i in minuspoint1names){
  x <- as.numeric(y[,i])
  minuspoint1 <- minuspoint1 + x
}

zeronames <- as.character(grep("zero", colnames(y), value=TRUE))
zero <- 0
for (i in zeronames){
  x <- as.numeric(y[,i])
  zero <- zero + x
}

point1names <- as.character(grep("^point1", colnames(y), value=TRUE))
point1 <- 0
for (i in point1names){
  x <- as.numeric(y[,i])
  point1 <- point1 + x
}

point2names <- as.character(grep("^point2", colnames(y), value=TRUE))
point2 <- 0
for (i in point2names){
  x <- as.numeric(y[,i])
  point2 <- point2 + x
}

point3names <- as.character(grep("^point3", colnames(y), value=TRUE))
point3 <- 0
for (i in point3names){
  x <- as.numeric(y[,i])
  point3 <- point3 + x
}


library(ggplot2)
library(viridis)
library(scales)
options(scipen=999)

# aggregate the dimensions together:
corbands_dim_F7 <- list()
corbands_dim_F7$minus1 <- as.numeric(corbands_dim$minus1.Rdata[1,])
corbands_dim_F7$minuspoint9 <- as.numeric(corbands_dim$minuspoint9.Rdata[1,])
corbands_dim_F7$minuspoint8 <- as.numeric(corbands_dim$minuspoint8.Rdata[1,])
corbands_dim_F7$minuspoint7 <- as.numeric(corbands_dim$minuspoint7.Rdata[1,])
corbands_dim_F7$minuspoint6 <- as.numeric(corbands_dim$minuspoint6.Rdata[1,])
corbands_dim_F7$minuspoint5 <- as.numeric(corbands_dim$minuspoint5.Rdata[1,])
corbands_dim_F7$minuspoint4 <- as.numeric(corbands_dim$minuspoint4.Rdata[1,])
corbands_dim_F7$minuspoint3 <- minuspoint3
corbands_dim_F7$minuspoint2 <- minuspoint2
corbands_dim_F7$minuspoint1 <- minuspoint1
corbands_dim_F7$zero <- zero
corbands_dim_F7$point1 <- point1
corbands_dim_F7$point2 <- point2
corbands_dim_F7$point3 <- point3
corbands_dim_F7$point4 <- as.numeric(corbands_dim$point4.Rdata[1,])
corbands_dim_F7$point5 <- as.numeric(corbands_dim$point5.Rdata[1,])
corbands_dim_F7$point6 <- as.numeric(corbands_dim$point6.Rdata[1,])
corbands_dim_F7$point7 <- as.numeric(corbands_dim$point7.Rdata[1,])
corbands_dim_F7$point8 <- as.numeric(corbands_dim$point8.Rdata[1,])
corbands_dim_F7$point9 <- as.numeric(corbands_dim$point9.Rdata[1,])

numbers.df <- as.data.frame(corbands_dim_F7)
numbers.df <- as.data.frame(t(numbers.df))
numbers.df$correlationband_dim <- as.character(rownames(numbers.df))
names(numbers.df) <- c("Number", "correlationband_dim")

# relabel the dimensions so they come out in order on the plot:
numbers.df$correlation.band <- as.character(c("-1 to -0.9", "-0.9 to -0.8", 
                                              "-0.8 to -0.7", "-0.7 to -0.6",
                                              "-0.6 to -0.5", "-0.5 to -0.4",
                                              "-0.4 to -0.3", "-0.3 to -0.2",
                                              "-0.2 to -0.1", "-0.1 to 0",
                                              "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                              "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                              "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                              "0.9 to 1"))
numbers.df$correlation.band <- factor(numbers.df$correlation.band, 
                                      levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                               "-0.8 to -0.7", "-0.7 to -0.6",
                                               "-0.6 to -0.5", "-0.5 to -0.4",
                                               "-0.4 to -0.3", "-0.3 to -0.2",
                                               "-0.2 to -0.1", "-0.1 to 0",
                                               "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                               "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                               "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                               "0.9 to 1"))

# create plot:
jpeg(filename ="path/to/F7_numberofcorrelationsperband.jpg", width = 7, height = 5, units = 'in', res=300)
ggplot(numbers.df, aes(x=correlation.band, y=Number)) + 
  geom_bar(stat="identity", position=position_dodge(), fill="#238A8DFF") +
  scale_fill_viridis(discrete = T, begin=0.2, end=0.6)+
  geom_text(aes(label=Number), position=position_dodge(width=0), vjust=0.5, angle = -90, hjust=1.02, size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0))+
  coord_cartesian(ylim = c(0,30000000000))+
  scale_y_continuous("Number of correlations", sec.axis = sec_axis(~./sum(numbers.df$Number)*100, name = "Percentage"), labels = comma)+
  labs(title="Distribution of pairwise correlation\nvalues in ARIES at 7 years old",x="Correlation from -1 to 1, in increments of 0.1", y="Number of correlations")
dev.off()

# and the total number of correlations (for sanity check) is...

totalcorrelations <- sum(numbers.df$Number)




