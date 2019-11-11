### WGCNA code for AIRES
## F7 timepoint
# Standalone network construction
# adapted from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

## DNAm betas have been normalised, filtered, and adjusted for age, sex and cell counts

# Loading data

# install.packages(WGCNA)
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

## DATA LOADING
# Load the AIRES F7 methylation dataset
lnames = load("/path/to/F7_450k/F7readyfor450kcor.Rdata")
lnames
F7Data.sort[1:10,1:10]

## 2. DATA CLEANING

# If there are columns etc that are not needed (eg probe quality) 
# then remove here. Only methyl values in the matrix.

F7data1 = as.data.frame(F7Data.sort)
colnames(F7data1)[1:5]
# CpGs
rownames(F7data1)[1:5]
# samples
F7data1[1:10,1:10]

# We first check for genes and samples with too many missing values.
# If the last statement returns TRUE, all genes have passed the cuts. 
gsg = goodSamplesGenes(F7data1, verbose = 3);
gsg$allOK

# If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(F7data1)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(F7data1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  F7data1 = F7data1[gsg$goodSamples, gsg$goodGenes]
}
# How many CpGs are left?
dim(F7data1)

#======================================================
#  Code chunk 7
# We now read in the trait data and match the samples
#=====================================================================================


traitData = read.csv("/path/to/F7NetworkTraits.txt", header = T);
dim(traitData)
names(traitData)

# Form a data frame (F7Samples) analogous to expression data that will hold the clinical traits.
# first create a vector of the sample IDs in the metylation data:
F7Samples = rownames(F7data1);
# make a vector traitRows that identifies sample IDs that are present in both the 
# methylation and trait datasets:
traitRows = match(F7Samples, traitData$ID);
length(traitRows)
#Remove the 4 NAs which pop up where the siblings were:
traitRows <- na.omit(traitRows)
length(traitRows)

# makes a vector datTraits which is the trait data read in from the csv (traitData),
# the samples are restricted to those who also have methylation data,
# and the sample IDs are removed as a column and made into the rownames:

datTraits = traitData[traitRows,]
head(datTraits)
dim(datTraits)
datTraits = traitData[traitRows, -1];
head(datTraits)
dim(datTraits)
rownames(datTraits) = traitData[traitRows, 1];
head(datTraits)
dim(datTraits)

collectGarbage();

#=====================================================================================

# finally, make F7data1 the same length as datTraits (as F7 Data currently has more
# samples)
#=====================================================================================

dim(F7data1)
traitSamples = rownames(datTraits);
length(traitSamples)
F7data1 <- F7data1[traitSamples,]
dim(F7data1)

head(datTraits)

save(F7data1, datTraits, file = "/path/to/F7_WGCNA/F7_full450k_WGCNA.RData")

#=====================================================================================

# Sample clustering

#=====================================================================================

# Next we cluster the samples (in contrast to clustering genes that will 
# come later) to see if there are any obvious outliers.

sampleTree = hclust(dist(F7data1), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "/path/to/F7_WGCNA/sampleClusteringF7_WGCNA_full450k.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,4,0))
plot(sampleTree, main = "Sample clustering to detect outliers, F7,\n outliers removed, and adjusted for covariates in linear model", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

#================================================
#  Code chunk 6
# Are there outlier samples?
# Choose a height cut that will remove the offending sample, say [15] (the red 
# line in the plot), and use a branch cut at that height.
#=====================================================================================


# Plot a line to show the cut 
sizeGrWindow(12,9)
pdf(file = "/path/to/F7_WGCNA/sampleClusteringF7WithProposedHeightCut-full450k_WGCNA.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers, F7, with height cut", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2);
abline(h = 40, col = "red")
dev.off()

# Determine cluster under the line 
clust = cutreeStatic(sampleTree, cutHeight = 40, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
F7DataExperimental = F7data1[keepSamples, ]
nGenes = ncol(F7DataExperimental)
nSamples = nrow(F7DataExperimental)
F7DataExperimental[1:10,1:10]

sampleTree = hclust(dist(F7DataExperimental), method = "average");

sizeGrWindow(12,9)
pdf(file = "/path/to/F7_WGCNA/sampleClusteringF7WithHeightCut-full450k_WGCNA.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers, F7, with height cut", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2);
abline(h = 42.5, col = "red")
dev.off()

F7data1 <- F7DataExperimental
names <- rownames(F7data1)
length(names)
datTraits <- datTraits[match(names, rownames(datTraits)),]

save(F7data1, datTraits, file = "/path/to/F7_WGCNA/F7_heightcut_WGCNA_full450k.RData")


# F7WGCNAdata now contains the expression data ready for network analysis.


#**************************************
#**************************************
## PHASE 2 - CONSTRUCTING THE NETWORK 
#**************************************
#**************************************

Sys.time()

#workingDir = ".";
#setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "/path/to/F7_WGCNA/F7_heightcut_WGCNA_full450k.RData");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
# Choosing the soft threshold power
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:20))
# Call the network topology analysis function
sft = pickSoftThreshold(F7data1, powerVector = powers, networkType = "signed", verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "/path/to/F7_WGCNA/F7_BWScaleFreeTopologyPlot_full450k.pdf", width = 12, height = 9);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence, F7 full 450k"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()
# Mean connectivity as a function of the soft-thresholding power
sizeGrWindow(9, 5)
pdf(file = "/path/to/F7_WGCNA/F7_BWConnectivityPlot_full450k.pdf", width = 12, height = 9);
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity, F7 full 450k"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
