library(circlize)
library(viridis)

load("/path/to/F7_point9.Rdata")
load("/path/to/probeDetails_meffil.Rdata")

dat$Var1.position <- probeDetails$position[match(dat$Var1, rownames(probeDetails))]
dat$Var2.position <- probeDetails$position[match(dat$Var2, rownames(probeDetails))]
dat$absdistance <- abs(dat$Var1.position - dat$Var2.position)
dat$Var1.chromosome <- probeDetails$chromosome[match(dat$Var1, rownames(probeDetails))]
dat$Var2.chromosome <- probeDetails$chromosome[match(dat$Var2, rownames(probeDetails))]
dat$cis_or_trans <- as.character(c("trans"))
dat$cis_or_trans[dat$absdistance <= 1000000 & dat$Var1.chromosome == dat$Var2.chromosome] <- "cis"
dat$cis_or_trans[dat$Var1.chromosome == dat$Var2.chromosome] <- "cis"
print(head(dat))

trans_chr <- dat[dat$cis_or_trans == "trans",]

bedforcircos1 <- trans_chr[,c(7,4)]
names(bedforcircos1) <- c("chr","start")
bedforcircos1$start <- bedforcircos1$start - 200000
bedforcircos1$end <- bedforcircos1$start + 200000
bedforcircos2 <- trans_chr[,c(8,5)]
names(bedforcircos2) <- c("chr","start")
bedforcircos2$start <- bedforcircos2$start - 200000
bedforcircos2$end <- bedforcircos2$start + 200000

circos.clear()

jpeg(filename = "/path/to/F7_trans_circos.jpg", width = 7, height = 5, units = 'in', res=600)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(plotType = c("axis", "labels"), chromosome.index = paste0("chr", c(1:22)))
circos.track(ylim = c(0, 1), 
             bg.col = viridis(22), 
             bg.border = NA, track.height = 0.05)
circos.genomicLink(bedforcircos1, bedforcircos2, col = "black", 
                   border = NA)
dev.off()