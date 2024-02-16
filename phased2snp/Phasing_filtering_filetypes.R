################## 
## SNP Pipeline ## 
##################


# Phased haplotypes are converted to unphased data where ambiguity codes indicate heterozygosity. Ambiguity codes in the phased data (including 'N'), as well as gaps are recoded as '?'. This can be modified if desired.

## for input, have each locus in a separate fasta file. Each individual should be on consecutive lines. Assumes haplotypes are named as "..seq1" and "..seq2"

## Assumes that locus files have an 'L' and the number somewhere in the file name (e.g. "..L15...fasta")

##devtools::install_github("bbanbury/phrynomics") ##in case need to install again

library(ape)
library(phrynomics)

## set location of script with modified functions: "phrynomics_SRKFunctions.R" 
#source("G:/My Drive/OUwork/Research/Pkentucki_FileTypes/phrynomics_myFunctions.R")
source("/Volumes/GoogleDrive-114548620038544830238/My Drive/OUwork/Research/Pkentucki_FileTypes/phrynomics_myFunctions.R")

## set your working directory, which has each locus in its own fasta file
#setwd("G:/My Drive/OUwork/Research/Pkentucki_FileTypes/FASTA_filtered")
setwd("/Volumes/GoogleDrive-114548620038544830238/My Drive/OUwork/Research/Pkentucki_FileTypes/FASTA_filtered")

## set missing data threshold
missingAllowed = 21        # adjust to specify maximum number of individuals for each site that can have missing data. For loop, this will test all values from 0 through missingAllowed




###########################################################
#### Pipeline begins below: run entire rest of script #####
###########################################################

## ---- Unphasing Sequence Data ----

## to get taxa names:
x <- list.files()
loc1 <- read.FASTA(x[1]) # read any locus
seq1 <- grep(".+seq1", names(loc1), value=TRUE) # get names of just seq1
taxa <- gsub("_seq1", " ", seq1) # drop the _seq1 from the name, replace with space; space following name important for making phylip file later     

## **this will take a while to run**
## note that, as is, this turns anything that is not a base to a '?', including gaps, '?', and all ambiguity codes
x <- list.files()
datalist <- list()
for(z in 1:length(x)){
  file1 <- read.FASTA(x[z])
  matrix <- as.matrix(file1)
  fullMatrix <- as.character(matrix)
  df1 <- fullMatrix
  newDF <- matrix(nrow = nrow(df1)/2, ncol = ncol(df1), 0)
  for (i in 1:ncol(df1)){
    for (j in 1:length(newDF[,i])){
      newDF[j,i] <- paste0(df1[(j*2)-1,i], df1[j*2,i])  # concatenates each individual's two bases for each site
    }                                                   # e.g. two lines each "A" and "A" becomes one line "AA"
  }
  Na <- newDF
  unphased <- data.frame(apply(Na, c(1,2), function(x){ # make base pairs uppercase (optional aesthetic choice)
    toupper(x)
  }))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub(".M|.R|.W|.S|.Y|.K|.V|.H|.D|.B|M.|R.|W.|S.|Y.|K.|V.|H.|D.|B.|.N|N.|.\\?|\\?.|.-|-.", "?", x)}))
  #unphased <- data.frame(apply(unphased, c(1,2), function(x){
  #  gsub(".M|.R|.W|.S|.Y|.K|.V|.H|.D|.B|M.|R.|W.|S.|Y.|K.|V.|H.|D.|B.|.N|N.", "N", x)}))  ## adjust how missing data coded, if desired
  #unphased <- data.frame(apply(unphased, c(1,2), function(x){
  #  gsub(".\\?|\\?.|.-|-.", "?", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("AA", "A", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("TT", "T", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("CC", "C", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("GG", "G", x)}))
  ## Ambiguity for heterozygotes
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("AC|CA", "M", x)}))      
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("AG|GA", "R", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("AT|TA", "W", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("CG|GC", "S", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("CT|TC", "Y", x)}))
  unphased <- data.frame(apply(unphased, c(1,2), function(x){
    gsub("GT|TG", "K", x)}))
  rownames(unphased) <- taxa
  datalist[[z]] <- unphased
}

## add row names
for (i in 1:length(datalist)){
  rownames(datalist[[i]]) <- taxa
}

## define function to create a table with phylip format header
write.table_with_header <- function(x, file, header){
  cat(header, '\n',  file = file)
  write.table(x, file, append = T, sep = "", row.names = T, quote = F, col.names = F)
}


## create directory and write phylip files to it
dir.create("../unphasedphylip_filtered")
for (i in 1:length(datalist)){
  write.table_with_header(datalist[[i]], paste0("../unphasedphylip_filtered/", x[i], ".phylip"), header = paste(nrow(datalist[[i]]), ncol(datalist[[i]]), sep= "\t"))
}

## move to the new directory with phylip files
setwd("../unphasedphylip_filtered/")


## create new folder, convert phylip to nexus (or unhash the two lines below for fasta files)
dir.create("../unphased_nexus_filtered/")
dir.create("../unphased_fasta_filtered")  
x <- list.files()
for (i in 1:length(x)) {
  y <- ape::read.dna(x[i], format = "sequential")
  ape::write.nexus.data(y, interleaved = F, paste0("../unphased_nexus_filtered/", x[i], ".nexus"))
  ape::write.FASTA(y, file = paste0("../unphased_fasta_filtered/", x[i], ".fasta"))
}

## end in new folder with nexus files 
setwd("../unphased_nexus_filtered/")

## Proceed to phrynomics SNP pipeline



## ---- Phyrnomics pipeline ----
## data should be unphased, with ambiguity = heterozygosity. Because of this, ambiguity codes as well as gaps were all coded as question marks for the unphased data, while phased heterozygosity was turned into a single ambiguity code. Have each locus in a separate nexus file. First provide a list of taxa and loci below

## note that when phrynomics takes a single SNP per locus, it does so for SNPs with best coverage, and if >1 are equal, it takes a random SNP. This is OK, but note that it is technically not a random SNP per locus 
#setwd("G:/My Drive/OUwork/Research/Pkentucki_FileTypes/kentuckiunphased_nexus_noparalogs_above1dist")

## create and set output folder
dir.create("../SNP-data-phrynomics_test/")
outputFolder <- "../SNP-data-phrynomics_test/"
## or, set output folder of your choosing
#outputFolder <- "G:/My Drive/OUwork/Research/Pkentucki_FileTypes/unfiltered-SNP-data"

## load one locus to get taxa names:
x <- list.files()
loc1 <- read.nexus.data(x[1]) 
taxa <- names(loc1)

## to get locus names; assuming files contain name of locus as L1, L2...:
locusNum <- regmatches(x, regexpr("L\\d+", x))

## parameters
ntax = length(loc1)        # number of individuals
nloci = length(locusNum)   # number of loci 
loci <- locusNum           # Locus names; seems to not matter because of how Phrynomics processes data



## ---- perform SNP selection ----
library(adegenet)

## phrynomics accepts an odd format of SNP data; have your unphased loci in individual nexus files; for my data, all missing characters I coded as '?' (including gaps and ambiguity codes that didn't mean heterozygosity); if desired, see 'data(fakeData)' for example of phrynomics raw data format
x <- list.files()
newMat <- matrix(nrow = ntax, ncol = nloci, NA) # empty matrix
newMat <- data.frame(newMat, row.names = taxa)
for (i in 1:ncol(newMat)){
  nex <- ReadSNP(x[i], fileFormat = "nex")
  data <- unlist(nex$data$locus1)       
  newMat[,i] <- data
}
rownames(newMat) <- taxa
colnames(newMat) <- loci

snpData <- ReadSNP(newMat) # read SNP data

## remove non-variant sites, then non-biallelic, then non-parsimony-informative
notVariant <- RemoveInvariantSites(snpData, chatty = TRUE) ### this step takes a while
notBinary <- RemoveNonBinary(notVariant)
informativeOnly <- RemoveUninformative(notBinary)

## First reduce dataset by maximum most relaxed missing data threshold 
reduced <- ReduceMinInd(informativeOnly, threshold = ntax - missingAllowed, calc = 'sites')
single_Snp <- MODTakeSingleSNPfromEachLocus(reduced)$snpdata

## ---- create list containing datasets with varying missing data allowed ----
all.list <- list()
one.list <- list()
all.list[[1]] <- reduced
one.list[[1]] <- single_Snp
missingLoop <- c(0:(missingAllowed-1))
allowed = rev(missingLoop)
for (i in 1:length(missingLoop)) {
  red <- ReduceMinInd(all.list[[i]], threshold = ntax - allowed[i], calc = 'sites')
  #sin <- MODTakeSingleSNPfromEachLocus(red)$snpdata
  one.red <- ReduceMinInd(one.list[[i]], threshold = ntax - allowed[i], calc = 'sites')
  all.list[[i+1]] <- red
  one.list[[i+1]] <- one.red
}


## ---- Make output files ----

## ---- ALL SNPs ----
set <- c(missingAllowed, allowed) #  descending sequence from maximum missing to zero
for (i in 1:length(all.list)){
  run <- set[i]
  path.output <- paste0(outputFolder, run, "_missing_allowed/")
  dir.create(path.output)
  allSnp <- TranslateBases(all.list[[i]], ordered = T)          # turn "ACTG" format to "0,1,2"
  allsnpTable <- strsplit(allSnp$data$locus1, split = "") # gives each base its own column
  alldf <- data.frame(matrix(unlist(allsnpTable), nrow=ntax, byrow=T),stringsAsFactors=FALSE)
  rownames(alldf) <- taxa
  write.table(alldf, file= paste0(path.output, run, "_missing_allSnpTable_012.txt") ,quote = F, sep = "\t", col.names = F)
  ## Structure file:
  structure <- MakeStructureFormat(all.list[[i]])
  pops <- rep(1:ntax, each = 2)
  structure_pops <- cbind(pops, structure)
  colnames(structure_pops)[1] <- "\t"
  structure_pops[structure_pops == "?"] <- "-9"
  write.table(structure_pops, file= paste0(path.output, run, "_missing_structure_allSnps.txt"), quote=F, sep = "\t", col.names = T)
  ## Phylip file:
  WriteSNP(all.list[[i]], file = paste0(path.output, run, "_missing_allSnps.phy"), format = "phylip")
  ## Nexus file:
  WriteSNP(all.list[[i]], file = paste0(path.output, run, "_missing_allSnps.nex"), missing = "?", format = "nexus")
  
  ## ---- One SNP per locus ----
  oneSnp <- TranslateBases(one.list[[i]], ordered = T)
  onesnpTable <- strsplit(oneSnp$data$locus1, split = "")
  onedf <- data.frame(matrix(unlist(onesnpTable), nrow=ntax, byrow=T),stringsAsFactors=FALSE)
  rownames(onedf) <- taxa
  write.table(onedf, file= paste0(path.output, run, "_missing_oneSnpTable_012.txt") ,quote = F, sep = "\t", col.names = F)
  ## Structure file:
  structure <- MakeStructureFormat(one.list[[i]])
  pops <- rep(1:ntax, each = 2)
  structure_pops <- cbind(pops, structure)
  colnames(structure_pops)[1] <- "\t"
  structure_pops[structure_pops == "?"] <- "-9"
  write.table(structure_pops,  file= paste0(path.output, run, "_missing_structure_oneSnp.txt"), quote=F, sep = "\t", col.names = T)
  ## Phylip file:
  WriteSNP(one.list[[i]], file = paste0(path.output, run, "_missing_oneSnp.phy"), format = "phylip")
  ## Nexus file:
  WriteSNP(one.list[[i]], file = paste0(path.output, run, "_missing_oneSnp.nex"), missing = "?", format = "nexus")
  ## SNAPP file:
  WriteSNP(oneSnp, file = paste0(path.output, run, "_missing_SNAPP_oneSnp.nex"), missing = "?", format = "nexus")
}

#########################################################
## Calculate some summary statistics about the pipeline
#########################################################

## change to output folder
setwd(outputFolder)

## ---- Missing data per individual given missing data allowed----
missingData <- matrix(nrow = length(taxa), ncol = missingAllowed+1)
for (i in 1:missingAllowed+1) {
  missing <- i - 1
  snps <- read.table(paste0("./", missing, "_missing_allowed/", missing, "_missing_onesnpTable_012.txt"), header = F, row.names = 1)
  for (j in 1:length(taxa)) {
    missingData[j,i] <- length(which(snps[j,] == "?"))
  }
}
colnames(missingData) <- seq(0,missingAllowed)
rownames(missingData) <- taxa
#rownames(missingData) <- sample_ID
write.table(missingData, "missingDataIndividual.txt", sep = '\t', quote = F)


## ---- Number of loci remaining, given missing data allowed ----
numLoci_remaining <- matrix(nrow=length(one.list))
numLoci_remaining
rownames(numLoci_remaining) <- rev(seq(0,missingAllowed))
colnames(numLoci_remaining) <- "remainingLoci"
for (i in 1:missingAllowed+1) {
  numLoci_remaining[i,] <- (one.list[[i]]$nsites)
}
write.table(numLoci_remaining, "numberLoci_remaining.txt", sep = '\t', quote = F)

