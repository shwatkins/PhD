#=====================================================================================
#
#  Section 3 - Code chunk 1. F7-HighV-lm-noout-power 2-Deep Split 2
#
#=====================================================================================
# adapted from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

Sys.time()

#***********************
#PHASE 3 - RELATING MODULES TO EXTERNAL INFORMATION AND IDENTIFYING 
#                       IMPORTANT GENES
#***********************
#***********************
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "/path/to/F7_heightcut_WGCNA_full450k.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "/path/to/F7-BW-P7-ModS30-DS2-PAM-signed.RData");
lnames
ls()

### 1 - relating the modules to clinical traits ####
#=====================================================================================
#
#  Code chunk 2
# correlating modules with traits
#=====================================================================================

datTraits[1:10,]
F7data1[1:10,1:10]
blockwiseMEs[1:10,1:10]

samples <- as.character(rownames(datTraits))
length(samples)
head(datTraits)

# Define numbers of genes and samples
nGenes = ncol(F7data1);
nSamples = nrow(F7data1);
# Recalculate MEs with color labels (this is needed for blockwise analyses)
MEs0 = moduleEigengenes(F7data1, bwModuleColors)$eigengenes
MEs = orderMEs(MEs0)

# Regression of traits and MEs
EG <- names(MEs)
EG
TraitNames <- names(datTraits)

MEbmi_transfReg <- lm(as.matrix(MEs) ~ datTraits$bmi_transf + datTraits$householdclass)
MEbmi_transfReg.sum <- summary(MEbmi_transfReg)
MEbmi_transfReg.sum
MEbmi_transfReg.t <- as.data.frame(lapply(MEbmi_transfReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEhouseholdclassReg <- lm(as.matrix(MEs) ~ datTraits$householdclass)
MEhouseholdclassReg.sum <- summary(MEhouseholdclassReg)
MEhouseholdclassReg.t <- as.data.frame(lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEmatage_monthsReg <- lm(as.matrix(MEs) ~ datTraits$matage)
MEmatage_monthsReg.sum <- summary(MEmatage_monthsReg)
MEmatage_monthsReg.t <- as.data.frame(lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEASTHMAReg <- lm(as.matrix(MEs) ~ datTraits$ASTHMA + datTraits$matage + datTraits$householdclass + datTraits$antenatalAsthma + datTraits$sustained)
MEASTHMAReg.sum <- summary(MEASTHMAReg)
MEASTHMAReg.sum
MEASTHMAReg.t <- as.data.frame(lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEsustainedReg <- lm(as.matrix(MEs) ~ datTraits$sustained)
MEsustainedReg.sum <- summary(MEsustainedReg)
MEsustainedReg.t <- as.data.frame(lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(3)]))

singleTvalues <- rbind.data.frame(MEmatage_monthsReg.t, MEbmi_transfReg.t, MEhouseholdclassReg.t, MEASTHMAReg.t, MEsustainedReg.t)
singleTvalues
names(singleTvalues) <- EG
rownames(singleTvalues) <- c("Maternal age", "bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")
TraitNames
# and for p values

MEbmi_transfReg.p <- (lapply(MEbmi_transfReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEhouseholdclassReg.p <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatage_monthsReg.p <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEASTHMAReg.p <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEsustainedReg.p <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(4)]))

singlePvalues <- rbind.data.frame(MEmatage_monthsReg.p, MEbmi_transfReg.p, MEhouseholdclassReg.p, MEASTHMAReg.p, MEsustainedReg.p)
singlePvalues
names(singlePvalues) <- EG
rownames(singlePvalues) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEbmi_transfReg.r2 <- (lapply(MEbmi_transfReg.sum, function(x) x$r.squared))
MEhouseholdclassReg.r2 <- (lapply(MEhouseholdclassReg.sum, function(x) x$r.squared))
MEmatage_monthsReg.r2 <- (lapply(MEmatage_monthsReg.sum, function(x) x$r.squared))
MEASTHMAReg.r2 <- (lapply(MEASTHMAReg.sum, function(x) x$r.squared))
MEsustainedReg.r2 <- (lapply(MEsustainedReg.sum, function(x) x$r.squared))
R2 <- rbind.data.frame(MEmatage_monthsReg.r2, MEbmi_transfReg.r2, MEhouseholdclassReg.r2, MEASTHMAReg.r2, MEsustainedReg.r2)
R2
names(R2) <- EG
rownames(R2) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEbmi_transfReg.beta <- (lapply(MEbmi_transfReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEhouseholdclassReg.beta <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatage_monthsReg.beta <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEASTHMAReg.beta <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEsustainedReg.beta <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(1)]))
singlebeta <- rbind.data.frame(MEmatage_monthsReg.beta, MEbmi_transfReg.beta, MEhouseholdclassReg.beta, MEASTHMAReg.beta, MEsustainedReg.beta)
singlebeta
names(singlebeta) <- EG
rownames(singlebeta) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEbmi_transfReg.se <- (lapply(MEbmi_transfReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEhouseholdclassReg.se <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatage_monthsReg.se <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEASTHMAReg.se <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEsustainedReg.se <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(2)]))
singleSE <- rbind.data.frame(MEmatage_monthsReg.se, MEbmi_transfReg.se, MEhouseholdclassReg.se, MEASTHMAReg.se, MEsustainedReg.se)
singleSE
names(singleSE) <- EG
rownames(singleSE) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

#=====================================================================================
#
#  Code chunk 3
# graphical representation of module-trait relationships:
#=====================================================================================


sizeGrWindow(12,9)
pdf(file="/path/to/betas_F7-P7-ModS30-DS2-signed-Module-traitRelationshipHeatmap-singleregressions_test.pdf", width = 12, height = 9);
# Will display regressions and their p-values
SingleTScoresM <- as.matrix(singleTvalues)
SinglePValuesM <- as.matrix(singlePvalues)
singleR2.m <- as.matrix(R2)
singlebeta.m <- as.matrix(singlebeta)
singleSE.m <- as.matrix(singleSE)
# textMatrix holds both the t scores and the p values as a single value so each 
# square has both values per regression
textMatrix2 =  paste(signif(singlebeta.m, 1), "\n\n",
                     signif(singleSE.m, 1), "\n\n",
                     signif(SingleTScoresM, 2), "\n\n(",
                     signif(SinglePValuesM, 2),")\n\n",
                     signif(singleR2.m, 1), sep = "");
# Checking data hasn't shuffled somewhere it shouldn't...
textMatrix2[1:10]
singleTvalues[1:5,1:5]
singlePvalues[1:5,1:5]

dim(textMatrix2) = dim(singleTvalues)
textMatrix2

# Extend margins to fit all labels
par(mar = c(9, 9, 6, 9));
labeledHeatmap(singleTvalues, xLabels = names(singleTvalues), yLabels = rownames(singleTvalues),
               xSymbols = names(singleTvalues),
               colorLabels = TRUE,
               colors = blueWhiteRed(100), zlim=c(-7,7),
               setStdMargins = FALSE,
               textMatrix = textMatrix2,
               cex.text = 0.5,
               main = "F7 methylation heatmap@P7, with each square displaying t-score (p-value)\n for regression of each trait separately on each module eigengene");
dev.off()


#################### adjust for Eos:

# load predicted cell counts:
load("/path/to/F7_houseman_withEos.Rdata")

houseman$Bcell <- as.numeric(as.character(houseman$Bcell))
houseman$CD4T <- as.numeric(as.character(houseman$CD4T))
houseman$CD8T <- as.numeric(as.character(houseman$CD8T))
houseman$Eos <- as.numeric(as.character(houseman$Eos))
houseman$Mono <- as.numeric(as.character(houseman$Mono))
houseman$Neu <- as.numeric(as.character(houseman$Neu))
houseman$NK <- as.numeric(as.character(houseman$NK))

samples <- as.character(rownames(datTraits))
length(samples)
houseman <- houseman[samples,]
datTraits$Eos <- houseman$Eos[match(rownames(datTraits),rownames(houseman))]
datTraits$Neu <- houseman$Neu[match(rownames(datTraits),rownames(houseman))]

class(datTraits$Eos)
head(datTraits)
datTraits$Eos <- as.numeric(datTraits$Eos)
datTraits$Neu <- as.numeric(datTraits$Neu)
head(datTraits)

MEASTHMAReg <- lm(as.matrix(MEs) ~ datTraits$ASTHMA + datTraits$matage + datTraits$householdclass + datTraits$antenatalAsthma + datTraits$sustained + datTraits$Eos + datTraits$Neu)
MEASTHMAReg.sum <- summary(MEASTHMAReg)
MEASTHMAReg.sum
MEASTHMAReg.t <- as.data.frame(lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(3)]))
MEsustainedReg.t <- as.data.frame(lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(3)]))
singleTvalues <- rbind.data.frame(MEmatage_monthsReg.t, MEbmi_transfReg.t, MEhouseholdclassReg.t, MEASTHMAReg.t, MEsustainedReg.t)
singleTvalues
names(singleTvalues) <- EG
rownames(singleTvalues) <- c("Maternal age", "bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEASTHMAReg.p <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(4)]))
singlePvalues <- rbind.data.frame(MEmatage_monthsReg.p, MEbmi_transfReg.p, MEhouseholdclassReg.p, MEASTHMAReg.p, MEsustainedReg.p)
singlePvalues
names(singlePvalues) <- EG
rownames(singlePvalues) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEASTHMAReg.r2 <- (lapply(MEASTHMAReg.sum, function(x) x$r.squared))
R2 <- rbind.data.frame(MEmatage_monthsReg.r2, MEbmi_transfReg.r2, MEhouseholdclassReg.r2, MEASTHMAReg.r2, MEsustainedReg.r2)
R2
names(R2) <- EG
rownames(R2) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEASTHMAReg.beta <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(1)]))
singlebeta <- rbind.data.frame(MEmatage_monthsReg.beta, MEbmi_transfReg.beta, MEhouseholdclassReg.beta, MEASTHMAReg.beta, MEsustainedReg.beta)
singlebeta
names(singlebeta) <- EG
rownames(singlebeta) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")

MEASTHMAReg.se <- (lapply(MEASTHMAReg.sum, function(x) x$coefficients[c(2), c(2)]))
singleSE <- rbind.data.frame(MEmatage_monthsReg.se, MEbmi_transfReg.se, MEhouseholdclassReg.se, MEASTHMAReg.se, MEsustainedReg.se)
singleSE
names(singleSE) <- EG
rownames(singleSE) <- c("Maternal age","bmi_transf", "Household\nclass", "Asthma", "Sustained\nmaternal\nsmoking")


sizeGrWindow(12,9)
pdf(file="/path/to/F7_traits_adjforeos.pdf", width = 12, height = 9);
# Will display regressions and their p-values
SingleTScoresM <- as.matrix(singleTvalues)
SinglePValuesM <- as.matrix(singlePvalues)
singleR2.m <- as.matrix(R2)
singlebeta.m <- as.matrix(singlebeta)
singleSE.m <- as.matrix(singleSE)
# textMatrix holds both the t scores and the p values as a single value so each 
# square has both values per regression
textMatrix2 =  paste(signif(singlebeta.m, 1), "\n\n",
                     signif(singleSE.m, 1), "\n\n",
                     signif(SingleTScoresM, 2), "\n\n(",
                     signif(SinglePValuesM, 2),")\n\n",
                     signif(singleR2.m, 1), sep = "");
# Checking data hasn't shuffled somewhere it shouldn't...
textMatrix2[1:10]
singleTvalues[1:5,1:5]
singlePvalues[1:5,1:5]
dim(textMatrix2) = dim(singleTvalues)
textMatrix2

# Extend margins to fit all labels
par(mar = c(9, 9, 6, 9));
labeledHeatmap(singleTvalues, xLabels = names(singleTvalues), yLabels = rownames(singleTvalues),
               xSymbols = names(singleTvalues),
               colorLabels = TRUE,
               colors = blueWhiteRed(100), zlim=c(-7,7),
               setStdMargins = FALSE,
               textMatrix = textMatrix2,
               cex.text = 0.5,
               main = "F7 methylation heatmap@P7, with each square displaying t-score (p-value)\n for regression of each trait separately on each module eigengene");
dev.off()

