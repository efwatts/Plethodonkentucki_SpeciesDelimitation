## Identify maximum heterozygosity to flag paralogs

# Loci with sites in which most/all individuals are heterozygous means paralogs maybe present
# This script identifies the per-site number of heterozygous individuals and reports the maximum per locus in a table ("maximumHet.txt")
# This does not guarantee paralogy, but it flags loci to examine further; Plus, SNPs w/ >0.5 heterozgyosity should probably be removed anyways


# check out the RAxML trees of flagged loci. In at least a couple cases so far for P. cinereus (L204, L211), the tree was essentially one major split, with each individual having one phased haplotype on each side of the split. Pretty good evidence for a gene duplicate. Alignments also show overabundance of heterozygous sites  


# have a folder of unphased loci where ambiguity codes mean heterozygosity
setwd("G:/My Drive/OUwork/Research/Pkentucki_Filemodifications/paralogs")

library(ape)

#amb.codes <- c("m", "r", "w", "s", "y", "k")
## count number of heterozygotes per site, and sum them
## returns vector of number of heterzogous individuals per site for a given locus
het.count <- function(site){
  m <- length(which(site == "m"))
  r <- length(which(site == "r"))
  w <- length(which(site == "w"))
  s <- length(which(site == "s"))
  y <- length(which(site == "y"))
  k <- length(which(site == "k"))
  SUM <- sum(m,r,w,s,y,k)
  return(SUM)
}


## calculate heterozygosity per site for each locus, and report max in a table
## i.e. if 39 individuals in data set, then "39" for a locus means every individual is heterozygous for  at least one position, and locus should be investigated further for paralogy
x <- dir(pattern = ".nexus")
locusNum <- regmatches(x, regexpr("L\\d+", x))
max.heterozygosity <- matrix(nrow = length(locusNum), ncol = 1)
row.names(max.heterozygosity) <- locusNum
colnames(max.heterozygosity) <- "max.num.heterozygotes"
for (i in 1:length(x)) {
  nex <- read.nexus.data(x[i])
  df <- as.data.frame(nex)
  heterozygosity <- apply(df, 1, het.count)
  max.heterozygosity[i,] <- max(heterozygosity)
}
## write table
write.table(max.heterozygosity, file = "maximumHet.txt", sep = "\t", quote = FALSE, col.names = TRUE)

## show max number of heterozygotes in R
## more heterozygotes => higher possibility of parology
maxHet.sorted <- max.heterozygosity[order(max.heterozygosity[,1], decreasing = TRUE),]
as.data.frame(maxHet.sorted)

