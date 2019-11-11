
#  4 Relating consensus module to external microarray sample traits
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
# adapted from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "/path/to/ARIES_BiB_cord_bothethnicities_Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "/path/to/ARIES_BiB_cord_bothethnicities_Consensus-Network_BW.RData");
lnames
exprSize = checkSets(multiExpr);
nSets = exprSize$nSets;
nSets
consMEs = bnet$multiMEs

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
setLabels = c("ARIES", "BiB_WB", "BiB_P")
names(Traits[[1]])
names(Traits[[1]]$data)

# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]], use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}
names(moduleTraitCor)
moduleTraitCor[[1]]
names(moduleTraitPvalue)
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
names(consMEs[[1]]$data) <- MEColorNames
names(consMEs[[1]]$data)
set = 1

datTraits_ARIES <- read.csv("/path/to/CordNetworkTraits.txt", header = T)
datTraits_ARIES <- datTraits_ARIES[datTraits_ARIES$aln %in% rownames(multiExpr[[set]]$data),]
dim(datTraits_ARIES)
head(datTraits_ARIES)
rownames(datTraits_ARIES) <- datTraits_ARIES$aln

# Regression of traits and MEs
EG <- names(consMEs[[1]]$data)
TraitNames <- names(datTraits_ARIES)

# Create two data frames - one with the summary t values and one with the 
# p values. 
tScores <- as.data.frame(lapply(lmSummary, function(x) x$coefficients[, c(3)]))
names(tScores)
names(tScores) <- EG
rownames(tScores)
tScores <- tScores[-1,]
pValues <- as.data.frame(lapply(lmSummary, function(x) x$coefficients[, c(4)]))
names(pValues)
names(pValues) <- EG
pValues <- pValues[-1,]

tScores

MEgestageReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$gestage_transf + datTraits_ARIES$sustained)
MEgestageReg.sum <- summary(MEgestageReg)
MEgestageReg.sum
MEgestageReg.t <- as.data.frame(lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEbwtReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$bwt + datTraits_ARIES$gestage_transf + datTraits_ARIES$matbmi_transf + datTraits_ARIES$matage_months + datTraits_ARIES$sustained + datTraits_ARIES$householdclass)
MEbwtReg.sum <- summary(MEbwtReg)
MEbwtReg.sum
MEbwtReg.t <- as.data.frame(lapply(MEbwtReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEmatbmiReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$matbmi_transf + datTraits_ARIES$householdclass)
MEmatbmiReg.sum <- summary(MEmatbmiReg)
MEmatbmiReg.t <- as.data.frame(lapply(MEmatbmiReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEhouseholdclassReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$householdclass)
MEhouseholdclassReg.sum <- summary(MEhouseholdclassReg)
MEhouseholdclassReg.t <- as.data.frame(lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEcordAsthmaReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$asthma7 + datTraits_ARIES$antenatalAsthma + datTraits_ARIES$householdclass + datTraits_ARIES$sustained + datTraits_ARIES$matage_months)
MEcordAsthmaReg.sum <- summary(MEcordAsthmaReg)
MEcordAsthmaReg.t <- as.data.frame(lapply(MEcordAsthmaReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEmatage_monthsReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$matage_months)
MEmatage_monthsReg.sum <- summary(MEmatage_monthsReg)
MEmatage_monthsReg.t <- as.data.frame(lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(3)]))
MEmatage_monthsReg.sum

MEsustainedReg <- lm(as.matrix(consMEs[[1]]$data) ~ datTraits_ARIES$sustained + datTraits_ARIES$householdclass)
MEsustainedReg.sum <- summary(MEsustainedReg)
MEsustainedReg.t <- as.data.frame(lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(3)]))
MEsustainedReg.sum

singleTvalues <- rbind.data.frame(MEgestageReg.t, MEbwtReg.t, MEmatbmiReg.t, MEhouseholdclassReg.t, MEcordAsthmaReg.t, MEmatage_monthsReg.t, MEsustainedReg.t)
singleTvalues
names(singleTvalues) <- EG
rownames(singleTvalues) <- c("Gestational age", "Birthweight", "Maternal BMI", "Household class", "Asthma@7", "Maternal age", "Sustained\nmaternal smoking")
TraitNames

MEgestageReg.p <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEbwtReg.p <- (lapply(MEbwtReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatbmiReg.p <- (lapply(MEmatbmiReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEhouseholdclassReg.p <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEcordAsthmaReg.p <- (lapply(MEcordAsthmaReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatage_monthsReg.p <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEsustainedReg.p <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(4)]))

singlePvalues <- rbind.data.frame(MEgestageReg.p, MEbwtReg.p, MEmatbmiReg.p, MEhouseholdclassReg.p, MEcordAsthmaReg.p, MEmatage_monthsReg.p, MEsustainedReg.p)
singlePvalues
names(singlePvalues) <- EG
rownames(singlePvalues) <- c("Gestational age", "Birthweight", "Maternal BMI", "Household class", "Asthma@7", "Maternal age", "Sustained\nmaternal smoking")

MEgestageReg.r2 <- (lapply(MEgestageReg.sum, function(x) x$r.squared))
MEbwtReg.r2 <- (lapply(MEbwtReg.sum, function(x) x$r.squared))
MEmatbmiReg.r2 <- (lapply(MEmatbmiReg.sum, function(x) x$r.squared))
MEhouseholdclassReg.r2 <- (lapply(MEhouseholdclassReg.sum, function(x) x$r.squared))
MEcordAsthmaReg.r2 <- (lapply(MEcordAsthmaReg.sum, function(x) x$r.squared))
MEmatage_monthsReg.r2 <- (lapply(MEmatage_monthsReg.sum, function(x) x$r.squared))
MEsustainedReg.r2 <- (lapply(MEsustainedReg.sum, function(x) x$r.squared))
R2 <- rbind.data.frame(MEgestageReg.r2, MEbwtReg.r2, MEmatbmiReg.r2, MEhouseholdclassReg.r2, MEcordAsthmaReg.r2, MEmatage_monthsReg.p, MEsustainedReg.r2)
R2
names(R2) <- EG
rownames(R2) <- c("Gestational age", "Birthweight", "Maternal BMI", "Household class", "Asthma@7", "Maternal age", "Sustained\nmaternal smoking")

MEgestageReg.beta <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEbwtReg.beta <- (lapply(MEbwtReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatbmiReg.beta <- (lapply(MEmatbmiReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEhouseholdclassReg.beta <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEcordAsthmaReg.beta <- (lapply(MEcordAsthmaReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatage_monthsReg.beta <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEsustainedReg.beta <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(1)]))
singlebeta <- rbind.data.frame(MEgestageReg.beta, MEbwtReg.beta, MEmatbmiReg.beta, MEhouseholdclassReg.beta, MEcordAsthmaReg.beta, MEmatage_monthsReg.beta, MEsustainedReg.beta)
singlebeta
names(singlebeta) <- EG
rownames(singlebeta) <- c("Gestational age", "Birthweight", "Maternal BMI", "Household class", "Asthma@7", "Maternal age", "Sustained\nmaternal smoking")

MEgestageReg.se <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEbwtReg.se <- (lapply(MEbwtReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatbmiReg.se <- (lapply(MEmatbmiReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEhouseholdclassReg.se <- (lapply(MEhouseholdclassReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEcordAsthma.se <- (lapply(MEcordAsthmaReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatage_monthsReg.se <- (lapply(MEmatage_monthsReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEsustainedReg.se <- (lapply(MEsustainedReg.sum, function(x) x$coefficients[c(2), c(2)]))
singleSE <- rbind.data.frame(MEgestageReg.se, MEbwtReg.se, MEmatbmiReg.se, MEhouseholdclassReg.se, MEcordAsthma.se, MEmatage_monthsReg.se, MEsustainedReg.se)
singleSE
names(singleSE) <- EG
rownames(singleSE) <- c("Gestational age", "Birthweight", "Maternal BMI", "Household class", "Asthma@7", "Maternal age", "Sustained\nmaternal smoking")

#=====================================================================================
#
#  Code chunk 3
# graphical representation of module-trait relationships:
# Regression with single variables 
#=====================================================================================

sizeGrWindow(12,9)
pdf(file="/path/to/ARIES_traits_consensus_bothethnicities_regression.pdf", width = 12, height = 9);
# Will display regressions and their p-values
SingleTScoresM <- as.matrix(singleTvalues)

abs(max(SingleTScoresM, na.rm = T))
max(abs(SingleTScoresM),na.rm = T)
SinglePValuesM <- as.matrix(singlePvalues)
singleR2.m <- as.matrix(R2)
singlebeta.m <- as.matrix(singlebeta)
singleSE.m <- as.matrix(singleSE)
# textMatrix holds both the t scores and the p values as a single value so each 
# square has both values per regression
textMatrix2 =  paste(signif(singlebeta.m, 1), "\n",
                     signif(singleSE.m, 1), "\n",
                     signif(SingleTScoresM, 2), "\n(",
                     signif(SinglePValuesM, 2),")\n",
                     signif(singleR2.m, 1), sep = "");
# Checking data hasn't shuffled somewhere it shouldn't...
textMatrix2[1:10]
singleTvalues[1:5,1:5]
singlePvalues[1:5,1:5]
# Transformations so that the heatmap is the right way round...
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
               main = "cord methylation heatmap, with each square displaying t-score (p-value) \nfor regression of each trait separately on each module eigengene");
dev.off()

###########################################################

# Plot the module-trait relationship table for set number 2
set = 2

load("/path/to/BiB_cord_WB_traits.Rdata")
ls()

samples <- as.character(rownames(multiExpr[[2]]$data))
length(samples)
BiB_cord_traits_WB <- BiB_cord_traits_WB[samples,]

head(BiB_cord_traits_WB)

MEColors = labels2colors(as.numeric(substring(names(consMEs[[2]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
names(consMEs[[2]]$data) <- MEColorNames
names(consMEs[[2]]$data)

EG <- names(consMEs[[2]]$data)
EG
TraitNames <- names(BiB_cord_traits_WB)

MEsmokeReg <- lm(as.matrix(consMEs[[2]]$data) ~ BiB_cord_traits_WB$smoke + BiB_cord_traits_WB$IMD)
MEsmokeReg.sum <- summary(MEsmokeReg)
MEsmokeReg.t <- as.data.frame(lapply(MEsmokeReg.sum, function(x) x$coefficients[c(2), c(3)]))
MEsmokeReg.t
class(MEsmokeReg.t)
# list
summary(MEsmokeReg)

MEmatBMIReg <- lm(as.matrix(consMEs[[2]]$data) ~ BiB_cord_traits_WB$matBMI)
MEmatBMIReg.sum <- summary(MEmatBMIReg)
MEmatBMIReg.t <- as.data.frame(lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEIMDReg <- lm(as.matrix(consMEs[[2]]$data) ~ BiB_cord_traits_WB$IMD)
MEIMDReg.sum <- summary(MEIMDReg)
MEIMDReg.t <- as.data.frame(lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEgestageReg <- lm(as.matrix(consMEs[[2]]$data) ~ BiB_cord_traits_WB$gestage + BiB_cord_traits_WB$smoke)
MEgestageReg.sum <- summary(MEgestageReg)
MEgestageReg.t <- as.data.frame(lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEmatageReg <- lm(as.matrix(consMEs[[2]]$data) ~ BiB_cord_traits_WB$hdl)
MEmatageReg.sum <- summary(MEmatageReg)
MEmatageReg.t <- as.data.frame(lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(3)]))

singleTvalues <- rbind.data.frame(MEsmokeReg.t, MEmatBMIReg.t, MEIMDReg.t, MEgestageReg.t, MEmatageReg.t)
singleTvalues
names(singleTvalues) <- EG
rownames(singleTvalues) <- c("Maternal\nsmoking", "Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")
TraitNames
# and for p values

MEsmokeReg.p <- (lapply(MEsmokeReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatBMIReg.p <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEIMDReg.p <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEgestageReg.p <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatageReg.p <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(4)]))

singlePvalues <- rbind.data.frame(MEsmokeReg.p, MEmatBMIReg.p, MEIMDReg.p, MEgestageReg.p, MEmatageReg.p)
singlePvalues
names(singlePvalues) <- EG
rownames(singlePvalues) <- c("Maternal\nsmoking", "Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEsmokeReg.r2 <- (lapply(MEsmokeReg.sum, function(x) x$r.squared))
MEmatBMIReg.r2 <- (lapply(MEmatBMIReg.sum, function(x) x$r.squared))
MEIMDReg.r2 <- (lapply(MEIMDReg.sum, function(x) x$r.squared))
MEgestageReg.r2 <- (lapply(MEgestageReg.sum, function(x) x$r.squared))
MEmatageReg.r2 <- (lapply(MEmatageReg.sum, function(x) x$r.squared))
R2 <- rbind.data.frame(MEsmokeReg.r2, MEmatBMIReg.r2, MEIMDReg.r2, MEgestageReg.r2, MEmatageReg.r2)
R2
names(R2) <- EG
rownames(R2) <- c("Maternal\nsmoking", "Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEsmokeReg.beta <- (lapply(MEsmokeReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatBMIReg.beta <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEIMDReg.beta <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEgestageReg.beta <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatageReg.beta <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(1)]))
singlebeta <- rbind.data.frame(MEsmokeReg.beta, MEmatBMIReg.beta, MEIMDReg.beta, MEgestageReg.beta, MEmatageReg.beta)
singlebeta
names(singlebeta) <- EG
rownames(singlebeta) <- c("Maternal\nsmoking", "Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEsmokeReg.se <- (lapply(MEsmokeReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatBMIReg.se <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEIMDReg.se <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEgestageReg.se <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatageReg.se <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(2)]))
singleSE <- rbind.data.frame(MEsmokeReg.se, MEmatBMIReg.se, MEIMDReg.se, MEgestageReg.se, MEmatageReg.se)
singleSE
names(singleSE) <- EG
rownames(singleSE) <- c("Maternal\nsmoking", "Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

#=====================================================================================
#
#  Code chunk 3
# graphical representation of module-trait relationships:
# Regression with single variables 
#=====================================================================================

sizeGrWindow(12,9)
pdf(file="/path/to/BiB_WB_consensus_traits_regressions.pdf", width = 12, height = 9);
# Will display regressions and their p-values
SingleTScoresM <- as.matrix(singleTvalues)

abs(max(SingleTScoresM, na.rm = T))
max(abs(SingleTScoresM),na.rm = T)
SinglePValuesM <- as.matrix(singlePvalues)
singleR2.m <- as.matrix(R2)
singlebeta.m <- as.matrix(singlebeta)
singleSE.m <- as.matrix(singleSE)
# textMatrix holds both the t scores and the p values as a single value so each 
# square has both values per regression
textMatrix2 =  paste(signif(singlebeta.m, 1), "\n",
                     signif(singleSE.m, 1), "\n",
                     signif(SingleTScoresM, 2), "\n(",
                     signif(SinglePValuesM, 2),")\n",
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
               main = "BiB cord WB methylation heatmap@P6, with each square displaying t-score (p-value)\n for regression of each trait separately on each module eigengene");
dev.off()


##############################################################

# Plot the module-trait relationship table for set number 3
set = 3

load("/path/to/BiB_cord_P_traits.Rdata")
ls()

samples <- as.character(rownames(multiExpr[[3]]$data))
length(samples)
BiB_cord_traits_P <- BiB_cord_traits_P[samples,]

head(BiB_cord_traits_P)

MEColors = labels2colors(as.numeric(substring(names(consMEs[[3]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
names(consMEs[[3]]$data) <- MEColorNames
names(consMEs[[3]]$data)

EG <- names(consMEs[[3]]$data)
EG
TraitNames <- names(BiB_cord_traits_P)

MEmatBMIReg <- lm(as.matrix(consMEs[[3]]$data) ~ BiB_cord_traits_P$matBMI)
MEmatBMIReg.sum <- summary(MEmatBMIReg)
MEmatBMIReg.t <- as.data.frame(lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEIMDReg <- lm(as.matrix(consMEs[[3]]$data) ~ BiB_cord_traits_P$IMD)
MEIMDReg.sum <- summary(MEIMDReg)
MEIMDReg.t <- as.data.frame(lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEgestageReg <- lm(as.matrix(consMEs[[3]]$data) ~ BiB_cord_traits_P$gestage)
MEgestageReg.sum <- summary(MEgestageReg)
MEgestageReg.t <- as.data.frame(lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(3)]))

MEmatageReg <- lm(as.matrix(consMEs[[3]]$data) ~ BiB_cord_traits_P$hdl)
MEmatageReg.sum <- summary(MEmatageReg)
MEmatageReg.t <- as.data.frame(lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(3)]))

singleTvalues <- rbind.data.frame(MEmatBMIReg.t, MEIMDReg.t, MEgestageReg.t, MEmatageReg.t)
singleTvalues
names(singleTvalues) <- EG
rownames(singleTvalues) <- c("Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")
TraitNames
# and for p values

MEmatBMIReg.p <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEIMDReg.p <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEgestageReg.p <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(4)]))
MEmatageReg.p <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(4)]))

singlePvalues <- rbind.data.frame(MEmatBMIReg.p, MEIMDReg.p, MEgestageReg.p, MEmatageReg.p)
singlePvalues
names(singlePvalues) <- EG
rownames(singlePvalues) <- c("Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEmatBMIReg.r2 <- (lapply(MEmatBMIReg.sum, function(x) x$r.squared))
MEIMDReg.r2 <- (lapply(MEIMDReg.sum, function(x) x$r.squared))
MEgestageReg.r2 <- (lapply(MEgestageReg.sum, function(x) x$r.squared))
MEmatageReg.r2 <- (lapply(MEmatageReg.sum, function(x) x$r.squared))
R2 <- rbind.data.frame(MEmatBMIReg.r2, MEIMDReg.r2, MEgestageReg.r2, MEmatageReg.r2)
R2
names(R2) <- EG
rownames(R2) <- c("Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEmatBMIReg.beta <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEIMDReg.beta <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEgestageReg.beta <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(1)]))
MEmatageReg.beta <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(1)]))
singlebeta <- rbind.data.frame(MEmatBMIReg.beta, MEIMDReg.beta, MEgestageReg.beta, MEmatageReg.beta)
singlebeta
names(singlebeta) <- EG
rownames(singlebeta) <- c("Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

MEmatBMIReg.se <- (lapply(MEmatBMIReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEIMDReg.se <- (lapply(MEIMDReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEgestageReg.se <- (lapply(MEgestageReg.sum, function(x) x$coefficients[c(2), c(2)]))
MEmatageReg.se <- (lapply(MEmatageReg.sum, function(x) x$coefficients[c(2), c(2)]))
singleSE <- rbind.data.frame(MEmatBMIReg.se, MEIMDReg.se, MEgestageReg.se, MEmatageReg.se)
singleSE
names(singleSE) <- EG
rownames(singleSE) <- c("Maternal\nBMI", "IMD", "Gestational\nage", "Maternal age")

#=====================================================================================
#
#  Code chunk 3
# graphical representation of module-trait relationships:
# Regression with single variables 
#=====================================================================================

sizeGrWindow(12,9)
pdf(file="/path/to/BiB_P_consensus_traits_regressions.pdf", width = 12, height = 9);
# Will display regressions and their p-values
SingleTScoresM <- as.matrix(singleTvalues)

abs(max(SingleTScoresM, na.rm = T))
max(abs(SingleTScoresM),na.rm = T)
SinglePValuesM <- as.matrix(singlePvalues)
singleR2.m <- as.matrix(R2)
singlebeta.m <- as.matrix(singlebeta)
singleSE.m <- as.matrix(singleSE)
# textMatrix holds both the t scores and the p values as a single value so each 
# square has both values per regression
textMatrix2 =  paste(signif(singlebeta.m, 1), "\n",
                     signif(singleSE.m, 1), "\n",
                     signif(SingleTScoresM, 2), "\n(",
                     signif(SinglePValuesM, 2),")\n",
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
               main = "BiB cord P methylation heatmap@P6, with each square displaying t-score (p-value)\n for regression of each trait separately on each module eigengene");
dev.off()

