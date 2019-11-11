
# with adjustment for cell counts

################################################################
# adapted from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

library(WGCNA)
enableWGCNAThreads()

# Taking the modules found at birth and seeing how they are preserved
# at 7 years
# load data (rows are samples and columns are probes)
load(file = "~/F7-lm-NoOut-heightcut-full450k.RData");
lnames = load(file = "~/cord_heightcut_WGCNA_full450k.RData");
cordModuleData = read.csv(file = "~/birth_GeneInfo.csv") # module assignment info in here

colorscord <- as.character(cordModuleData$moduleColor)

setLabels = c("ARIEScord", "ARIESF7");
# set the datasets as a list
multiExpr = list(ARIEScord = list(data = corddata1), ARIESF7 = list(data = F7data1));
# set the Horvath WB450k module colours as a list:
multiColor = list(ARIEScord = colorscord);
#multiColor

mp = modulePreservation(multiExpr, multiColor,
                        networkType = "signed", corFnc="bicor", 
                        referenceNetworks = 1,
                        nPermutations = 200,
                        randomSeed = 1,
                        maxModuleSize = 4000,
                        quickCor = 0,
                        verbose = 3)

# Save the results
save(mp, file = "~/cordModulePreservationInF7.Rdata")

# We now analyze the data. Isolate the observed statistics and their Z scores:
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# We look at the main output: the preservation medianRank and Zsummary statistics.
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# The numbers are nice, but (to paraphrase a saying) a picture is worth a thousand 
# numbers. We plot the preservation medianRank and Zsummary for the female modules 
# as a function of module size. The plotting code may seem a bit involved, 
# but it is worth going through

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="~/cordModulesInF7-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 4000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
dev.off();

# We now plot the density and connectivity statistics all in one plot. 
# We include the module quality measures for comparison:
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
pdf(file="~/cordModulesInF7-modulePreservation-moduleQualityAndPreservation.pdf", wi=10, h=5)
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 4000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()

data.frame(color = modColors[plotMods], label = labs)



