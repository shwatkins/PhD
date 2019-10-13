Sys.time()
load(file="/path/to/adjusted_data.Rdata")
F7Data.sort <- t(F7Data.sort)
class(F7Data.sort)

library(WGCNA)
Sys.time()

## number of CpG sites to calculate correlations with
## at one time
block.size <- 25000

## correlation function
CORFUN <- function(x,y) bicor(x,y, maxPOutliers = 0.05)

Sys.time()

## divide sites into <=block.size groups
# Total no of blocks (ceiling ensures the remainder is also a block)
n.blocks <- ceiling(nrow(F7Data.sort)/block.size)
n.blocks
# assign the rows to a block
blocks <- rep(1:n.blocks, each=block.size, len=nrow(F7Data.sort))
blocks <- split(1:nrow(F7Data.sort), blocks)
# which row does the block start at?
block.start <- sapply(blocks, head, n=1)
# Size of each block:
block.size <- sapply(blocks, length)

Sys.time()

## calculate correlations between each pair of blocks
## and save correlations to the file
system.time(sapply(1:length(blocks), function(b1) {
  cat(date(), b1, " of ", length(blocks), "\n")
  sapply(b1:length(blocks), function(b2) {
    cormat <- CORFUN(t(F7Data.sort[blocks[[b1]],]),
                t(F7Data.sort[blocks[[b2]],]))
    if(b1==b2) {
      print(paste("b1=",b1))
      print(paste("b2=",b2))
      print(cormat[1:5,1:5])
      cormat[lower.tri(cormat)] <- NA
      print(cormat[1:5,1:5])
    } else {
      print(paste("b1=",b1))
      print(paste("b2=",b2))
    }
    save(cormat, file=(paste("/path/to/bicor",b1,"-",b2,".Rdata", sep="")))
    NULL
    gc()
  })
  NULL
})) ## 25s for a 5320x1000 matrix
Sys.time()

