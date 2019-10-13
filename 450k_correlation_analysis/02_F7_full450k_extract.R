print(Sys.time())

library(reshape2)
library(WGCNA)
enableWGCNAThreads()

# list the paths to the files containing the correlation matrices:
files = list.files(path="/path/to/correlationmatrix/bicor", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

# Now extract the correlations between probe pairs, in bands of 0.1 of 
# correlation. 

dat <- data.frame(matrix(ncol=3))
names <- c("Var1", "Var2", "value")
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -1 & value < -0.9)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
sum(anyNA(dat))
save(dat, file="/path/to/correl_bands/minus1.Rdata")

gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.9 & value < -0.8)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint9.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.8 & value < -0.7)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint8.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.7 & value < -0.6)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint7.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.6 & value < -0.5)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint6.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.5 & value < -0.4)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint5.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.4 & value < -0.3)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint4.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.3 & value < -0.2)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint3.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.2 & value < -0.1)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint2.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= -0.1 & value < 0)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/minuspoint1.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0 & value < 0.1)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/zero.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.1 & value < 0.2)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point1.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.2 & value < 0.3)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point2.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.3 & value < 0.4)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point3.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.4 & value < 0.5)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point4.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.5 & value < 0.6)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point5.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.6 & value < 0.7)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point6.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.7 & value < 0.8)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point7.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.8 & value < 0.9)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
save(dat, file="/path/to/correl_bands/point8.Rdata")
gc()

dat <- data.frame(matrix(ncol=3))
colnames(dat) <- names
for (j in files) {
  load(j)
  print(dim(cormat))
  df <- subset(melt(as.matrix(cormat)), value >= 0.9 & value < 1)
  print(dim(df))
  dat <- rbind(dat, df)
  gc()
}
dat <- dat[complete.cases(dat), ]
dim(dat)
dat <- dat[!(dat$Var1 == dat$Var2),]
dim(dat)
save(dat, file="/path/to/correl_bands/point9.Rdata")
gc()

rm(dat)
