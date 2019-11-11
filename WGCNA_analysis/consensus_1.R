
# 2c Blockwise network construction
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
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "/path/to/ARIES_BiB_cord_bothethnicities_Consensus-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

#=====================================================================================
#
#  Code chunk 3
# Construct the blockwise network
#=====================================================================================

bnet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 45000, power = 7, minModuleSize = 30,
  deepSplit = 2, corType = "bicor", maxPOutliers = 0.05,
  networkType = "signed", TOMType = "signed",
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, saveTOMFileBase = "/path/to/ARIES_BiB_cord_bothethnicities_TOMfilebase",
  verbose = 5)


#=====================================================================================
#
#  Code chunk 4
# Comparing to automatic network
#=====================================================================================

# convert block labels to colors
bwLabels = bnet$colors
table(bwLabels)

bwColors = labels2colors(bwLabels)
table(bwColors)

consMEs = bnet$multiMEs;

#=====================================================================================
#
#  Code chunk 5
# Plotting dendrogram for each of the blocks
#=====================================================================================

# Here we show a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
pdf(file = "/path/to/ARIES_BiB_cord_bothethnicities_BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(bnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks)
  plotDendroAndColors(bnet$dendrograms[[block]], bwColors[bnet$blockGenes[[block]]],
                      "Module colors", 
                      main = paste("Gene dendrogram and module colors in block", block), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off()

save(bnet, bwLabels, bwColors, file = "/path/to/ARIES_BiB_cord_bothethnicities_Consensus-Network_BW.RData")

