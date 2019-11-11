F7MMGeneInfo <- read.table("/path/to/F7-GeneInfo.csv", header=T, sep = ",")
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)

modules <- as.character(F7MMGeneInfo$moduleColor)
modules <- unique(modules)
modules

for (i in modules){
  print(i)
  icol <- paste("MM.",i,sep = "")
  print(icol)
  over0.7 <- F7MMGeneInfo[F7MMGeneInfo[[icol]] > 0.7,]
  print(dim(over0.7))
  if (nrow(over0.7) > 1) {
    probeNames <- as.character(over0.7$TargetID)
    cpgs <- as.character(F7MMGeneInfo$TargetID)
    # GO testing with prior probabilities taken into account
    # Plot of bias due to differing numbers of CpG sites per gene
    gst <- gometh(sig.cpg = probeNames, all.cpg = cpgs, collection = "GO", plot.bias = TRUE, prior.prob = TRUE)
    # Total number of GO categories significant at 5% FDR
    print(table(gst$FDR<0.05))
    # Table of top GO results
    topGOterms <- topGO(gst)
    #topGOterms
    kegg <- gometh(sig.cpg = probeNames, all.cpg = cpgs, collection = "KEGG", prior.prob=TRUE)
    # Table of top KEGG results
    tableKEGG <- topKEGG(kegg)
    write.csv(topGOterms, file=paste("/path/to/F7_topGO_",i,"module.csv",sep = ""))
    write.csv(tableKEGG, file=paste("/path/to/F7_topKEGG_",i,"module.csv",sep = ""))
  } else {
    print(paste("no probes over kME 0.7 in ",i," module",sep = ""))
  }
  
}
