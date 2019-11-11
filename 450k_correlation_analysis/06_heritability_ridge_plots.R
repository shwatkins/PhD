# Heritability of probes

library(data.table)
library(ggplot2)
library(viridis)
library(ggridges)

# the input from the correlation analysis is a list of unique probes that are 
# in each correlation band/the high correlations and the full 394k.

# load file which has a list of unique probes for each correlation band
load("/path/to/corbands_cgnames_F7.Rdata")

# choose the group of probes to estimate heritability from
correlatedProbes <- corbands_cgnames_F7$point9

# Load in Hannon et al heritability estimates
# this file has all CpGs (including HLA and sex chromosomes), so we need to 
# slim it down:
Hannon.heritability <- fread("/path/to/HeritabilityEstimatesFromACMModel_All.csv")
names(Hannon.heritability) <- c("CpG", "A", "C", "E")

# genetic 

genetic.list <- list()

for (i in names(corbands_cgnames_F7)) {
  print(i)
  x <- Hannon.heritability[Hannon.heritability$CpG %in% corbands_cgnames_F7[[i]]]
  print("dim x:")
  print(dim(x))
  x <- x[,c(1,2)]
  print(names(x))
  x$band <- paste(i)
  print(head(x))
  names(x) <- c("CpG", "Heritability", "cor_band")
  print(head(x))
  genetic.list[[i]] <- x
}

dim(genetic.list)
names(genetic.list)

# plot genetic
genetic <- rbindlist(genetic.list)
dim(genetic)
head(genetic)

genetic$cor_band[genetic$cor_band == "minus1"] <- c("-1 to -0.9")
genetic$cor_band[genetic$cor_band == "minuspoint1"] <- c("-0.1 to 0")
genetic$cor_band[genetic$cor_band == "minuspoint2"] <- c("-0.2 to -0.1")
genetic$cor_band[genetic$cor_band == "minuspoint3"] <- c("-0.3 to -0.2")
genetic$cor_band[genetic$cor_band == "minuspoint4"] <- c("-0.4 to -0.3")
genetic$cor_band[genetic$cor_band == "minuspoint5"] <- c("-0.5 to -0.4")
genetic$cor_band[genetic$cor_band == "minuspoint6"] <- c("-0.6 to -0.5")
genetic$cor_band[genetic$cor_band == "minuspoint7"] <- c("-0.7 to -0.6")
genetic$cor_band[genetic$cor_band == "minuspoint8"] <- c("-0.8 to -0.7")
genetic$cor_band[genetic$cor_band == "minuspoint9"] <- c("-0.9 to -0.8")
genetic$cor_band[genetic$cor_band == "point1"] <- c("0.1 to 0.2")
genetic$cor_band[genetic$cor_band == "point2"] <- c("0.2 to 0.3")
genetic$cor_band[genetic$cor_band == "point3"] <- c("0.3 to 0.4")
genetic$cor_band[genetic$cor_band == "point4"] <- c("0.4 to 0.5")
genetic$cor_band[genetic$cor_band == "point5"] <- c("0.5 to 0.6")
genetic$cor_band[genetic$cor_band == "point6"] <- c("0.6 to 0.7")
genetic$cor_band[genetic$cor_band == "point7"] <- c("0.7 to 0.8")
genetic$cor_band[genetic$cor_band == "point8"] <- c("0.8 to 0.9")
genetic$cor_band[genetic$cor_band == "point9"] <- c("0.9 to 1")
genetic$cor_band[genetic$cor_band == "zero"] <- c("0 to 0.1")

genetic$cor_band <- as.factor(genetic$cor_band)
genetic$cor_band <- ordered(genetic$cor_band, 
                            levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                     "-0.8 to -0.7", "-0.7 to -0.6",
                                     "-0.6 to -0.5", "-0.5 to -0.4",
                                     "-0.4 to -0.3", "-0.3 to -0.2",
                                     "-0.2 to -0.1", "-0.1 to 0",
                                     "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                     "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                     "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                     "0.9 to 1"))

pdf(file = "~/Heritability_ridges_genetic_F7.pdf", width = 7, height = 5)
ggplot(genetic, aes(x = `Heritability`, y = `cor_band`, fill = ..x..)) +
  #  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, gradient_lwd = 1.) +
  geom_density_ridges_gradient(scale = 0.95, rel_min_height = 0.005) +
  #  scale_x_continuous(expand = c(0.01, 0)) +
  # scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis(name = "Heritability", option = "C") +
  labs(title = 'Contribution of heritability to DNAm variation', x="Proportion of DNAm variation accounted for by genetic factors", y="Correlation band") 
#theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank())
dev.off()


# unique environment 

unique.list <- list()

for (i in names(corbands_cgnames_F7)) {
  print(i)
  x <- Hannon.heritability[Hannon.heritability$CpG %in% corbands_cgnames_F7[[i]]]
  print("dim x:")
  print(dim(x))
  x <- x[,c(1,4)]
  print(names(x))
  x$band <- paste(i)
  print(head(x))
  names(x) <- c("CpG", "UniqueEnvironment", "cor_band")
  print(head(x))
  unique.list[[i]] <- x
}

dim(unique.list)
names(unique.list)

# plot unique.e
unique.e <- rbindlist(unique.list)
dim(unique.e)
head(unique.e)

unique.e$cor_band[unique.e$cor_band == "minus1"] <- c("-1 to -0.9")
unique.e$cor_band[unique.e$cor_band == "minuspoint1"] <- c("-0.1 to 0")
unique.e$cor_band[unique.e$cor_band == "minuspoint2"] <- c("-0.2 to -0.1")
unique.e$cor_band[unique.e$cor_band == "minuspoint3"] <- c("-0.3 to -0.2")
unique.e$cor_band[unique.e$cor_band == "minuspoint4"] <- c("-0.4 to -0.3")
unique.e$cor_band[unique.e$cor_band == "minuspoint5"] <- c("-0.5 to -0.4")
unique.e$cor_band[unique.e$cor_band == "minuspoint6"] <- c("-0.6 to -0.5")
unique.e$cor_band[unique.e$cor_band == "minuspoint7"] <- c("-0.7 to -0.6")
unique.e$cor_band[unique.e$cor_band == "minuspoint8"] <- c("-0.8 to -0.7")
unique.e$cor_band[unique.e$cor_band == "minuspoint9"] <- c("-0.9 to -0.8")
unique.e$cor_band[unique.e$cor_band == "point1"] <- c("0.1 to 0.2")
unique.e$cor_band[unique.e$cor_band == "point2"] <- c("0.2 to 0.3")
unique.e$cor_band[unique.e$cor_band == "point3"] <- c("0.3 to 0.4")
unique.e$cor_band[unique.e$cor_band == "point4"] <- c("0.4 to 0.5")
unique.e$cor_band[unique.e$cor_band == "point5"] <- c("0.5 to 0.6")
unique.e$cor_band[unique.e$cor_band == "point6"] <- c("0.6 to 0.7")
unique.e$cor_band[unique.e$cor_band == "point7"] <- c("0.7 to 0.8")
unique.e$cor_band[unique.e$cor_band == "point8"] <- c("0.8 to 0.9")
unique.e$cor_band[unique.e$cor_band == "point9"] <- c("0.9 to 1")
unique.e$cor_band[unique.e$cor_band == "zero"] <- c("0 to 0.1")

unique.e$cor_band <- as.factor(unique.e$cor_band)
unique.e$cor_band <- ordered(unique.e$cor_band, 
                             levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                      "-0.8 to -0.7", "-0.7 to -0.6",
                                      "-0.6 to -0.5", "-0.5 to -0.4",
                                      "-0.4 to -0.3", "-0.3 to -0.2",
                                      "-0.2 to -0.1", "-0.1 to 0",
                                      "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                      "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                      "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                      "0.9 to 1"))

pdf(file = "/path/to/Heritability_ridges_uniqueE_F7.pdf", width = 7, height = 5)
ggplot(unique.e, aes(x = `UniqueEnvironment`, y = `cor_band`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.005, gradient_lwd = 1.) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis(name = "Unique\nEnvironment", option = "C") +
  labs(title = 'Contribution of unique environment to DNAm variation', x="Proportion of DNAm variation accounted for by unique\nenvironmental factors (including measurement error)", y="Correlation band") 
#theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank())
dev.off()



# common environment 

common.list <- list()

for (i in names(corbands_cgnames_F7)) {
  print(i)
  x <- Hannon.heritability[Hannon.heritability$CpG %in% corbands_cgnames_F7[[i]]]
  print("dim x:")
  print(dim(x))
  x <- x[,c(1,3)]
  print(names(x))
  x$band <- paste(i)
  print(head(x))
  names(x) <- c("CpG", "CommonEnvironment", "cor_band")
  print(head(x))
  common.list[[i]] <- x
}

dim(common.list)
names(common.list)

# plot common.e
common.e <- rbindlist(common.list)
dim(common.e)
head(common.e)

common.e$cor_band[common.e$cor_band == "minus1"] <- c("-1 to -0.9")
common.e$cor_band[common.e$cor_band == "minuspoint1"] <- c("-0.1 to 0")
common.e$cor_band[common.e$cor_band == "minuspoint2"] <- c("-0.2 to -0.1")
common.e$cor_band[common.e$cor_band == "minuspoint3"] <- c("-0.3 to -0.2")
common.e$cor_band[common.e$cor_band == "minuspoint4"] <- c("-0.4 to -0.3")
common.e$cor_band[common.e$cor_band == "minuspoint5"] <- c("-0.5 to -0.4")
common.e$cor_band[common.e$cor_band == "minuspoint6"] <- c("-0.6 to -0.5")
common.e$cor_band[common.e$cor_band == "minuspoint7"] <- c("-0.7 to -0.6")
common.e$cor_band[common.e$cor_band == "minuspoint8"] <- c("-0.8 to -0.7")
common.e$cor_band[common.e$cor_band == "minuspoint9"] <- c("-0.9 to -0.8")
common.e$cor_band[common.e$cor_band == "point1"] <- c("0.1 to 0.2")
common.e$cor_band[common.e$cor_band == "point2"] <- c("0.2 to 0.3")
common.e$cor_band[common.e$cor_band == "point3"] <- c("0.3 to 0.4")
common.e$cor_band[common.e$cor_band == "point4"] <- c("0.4 to 0.5")
common.e$cor_band[common.e$cor_band == "point5"] <- c("0.5 to 0.6")
common.e$cor_band[common.e$cor_band == "point6"] <- c("0.6 to 0.7")
common.e$cor_band[common.e$cor_band == "point7"] <- c("0.7 to 0.8")
common.e$cor_band[common.e$cor_band == "point8"] <- c("0.8 to 0.9")
common.e$cor_band[common.e$cor_band == "point9"] <- c("0.9 to 1")
common.e$cor_band[common.e$cor_band == "zero"] <- c("0 to 0.1")

common.e$cor_band <- as.factor(common.e$cor_band)
common.e$cor_band <- ordered(common.e$cor_band, 
                             levels=c("-1 to -0.9", "-0.9 to -0.8", 
                                      "-0.8 to -0.7", "-0.7 to -0.6",
                                      "-0.6 to -0.5", "-0.5 to -0.4",
                                      "-0.4 to -0.3", "-0.3 to -0.2",
                                      "-0.2 to -0.1", "-0.1 to 0",
                                      "0 to 0.1", "0.1 to 0.2", "0.2 to 0.3",
                                      "0.3 to 0.4","0.4 to 0.5","0.5 to 0.6",
                                      "0.6 to 0.7","0.7 to 0.8","0.8 to 0.9",
                                      "0.9 to 1"))

pdf(file = "/path/to/Heritability_ridges_commonE_F7.pdf", width = 7, height = 5)
ggplot(common.e, aes(x = `CommonEnvironment`, y = `cor_band`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.005, gradient_lwd = 1.) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis(name = "Common\nEnvironment", option = "C") +
  labs(title = 'Contribution of common environment to DNAm variation', x="Proportion of DNAm variation accounted for by\ncommon environmental factors", y="Correlation band") 
#theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank())
dev.off()