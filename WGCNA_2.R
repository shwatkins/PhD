#**************************************
#**************************************
## PHASE 2 - CONSTRUCTING THE NETWORK 
#**************************************
#**************************************
# adapted from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

Sys.time()

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
#  Code chunk 3
# Block wise network construction
# to change deepSplit or PAM, this is the place to do it
# lots of other options in the help...
#=====================================================================================


bwnet = blockwiseModules(F7data1, maxBlockSize = 45000, corType = "bicor",
                         power = 7, networkType = "signed", TOMType = "signed", 
                         minModuleSize = 30, deepSplit = 2, maxPOutliers = 0.05,
                         pamStage = TRUE, pamRespectsDendro = FALSE,
                         reassignThreshold = 0, mergeCutHeight = 0.25, 
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "/path/to/F7_WGCNA/F7-BW-P7-ModS30-DS2-PAM-signed",
                         verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# convert block labels to colors
bwLabels = bwnet$colors
bwModuleColors = labels2colors(bwLabels)

# Print and save lists of the modules, 
# 1 labelled by number
table(bwLabels)
Modulesizes <- table(bwLabels)
write.csv(Modulesizes, file="/path/to/F7_WGCNA/ModulesF7-BW-P7-ModS30-DS2-PAM-signed.csv")
# 2 labelled by color
table(bwModuleColors)
ModulesizesAndColors <- table(bwModuleColors)
write.csv(ModulesizesAndColors, file="/path/to/F7_WGCNA/ModuleColorsAndSizesF7-BW-P7-ModS30-DS2-PAM-signed.csv")



#=====================================================================================
#
#  Code chunk 5
# Plot the 5 block dendrograms separately
#=====================================================================================

sizeGrWindow(6,6)
pdf(file = "/path/to/F7_WGCNA/Dendrograms-F7-BW-P7-ModS30-DS2-PAM-signed.pdf", width = 12, height = 9);
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))

nBWblocks = length(bwnet$dendrograms)

for (block in 1:nBWblocks)
  plotDendroAndColors(bwnet$dendrograms[[block]], bwModuleColors[bwnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors, F7-BW-P7-ModS30-DS2-PAM-signed, in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
dev.off();

#==================================
# SAVE
#==================================

blockwiseMEs = bwnet$MEs
geneTree = bwnet$dendrograms[[1]]


# Save module colors and labels for use in subsequent parts
save(blockwiseMEs, geneTree, bwLabels, bwModuleColors, bwnet, file = "/path/to/F7_WGCNA/F7-BW-P7-ModS30-DS2-PAM-signed.RData")