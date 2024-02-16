## Phrynomics modified and added funcitons
# my script will only work with these funcitons loaded

## ---- Modify/Add phyrnomics Functions ----
# some functions must be modified so that phyrnomics treats ambiguity codes as heterozygosity, not ambiguity

## ---- change IsVariable() to MODIsVar() ----
MODIsVar <- function(SNP){
  var <- FALSE
  bases <- c("A", "C", "G", "T")
  ambs <-  c("M", "R", "W", "S", "Y", "K")   # added list of heterozygote ambiguiuty codes
  basesInSNP <- bases[which(bases %in% SNP)]
  ambsInSNP <- ambs[which(ambs %in% SNP)]    
  if(all(SNP == " "))  # if row empty, true because it is a break between loci
    return(TRUE)
  if(length(basesInSNP) == 0 && length(ambsInSNP) == 0)  # all "N"
    return(FALSE)
  if(length(basesInSNP) == 1 && length(ambsInSNP) == 0)  # fixed
    return(FALSE)
  if(length(basesInSNP) == 1 && length(ambsInSNP) > 0)    # one base and ambiguity
    return(TRUE)
  if(length(basesInSNP) == 0 && length(ambsInSNP) > 0)    # one base and ambiguity
    return(TRUE)
  if(length(basesInSNP) > 1)  # more than one base 
    return(TRUE)
}

## ---- Change RemoveInvariantSites() to MODRemoveInvariantSitees() ----
MODRemoveInvariantSites <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1,]))
  initialLociLengths <- nchar(SNPdataset[1,])
  splitdata <- SplitSNP(SNPdataset)
  KeepVector <- apply(splitdata, 2, MODIsVar)  # changed to adopt my modified IsVar() function
  breaks <- which(splitdata[1,] == " ")
  newSNPdataset <- cSNP(splitdata, KeepVector=KeepVector, maintainLoci=TRUE)
  newsnps <- sum(nchar(newSNPdataset[1,]))
  if(chatty)
    message(paste("removed", snps-newsnps, "of", snps, "sites"))
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}

## ---- Modify Take single SNP ---- 
## there is an error in original function: did NOT choose SNP with the least missing data.
MODTakeSingleSNPfromEachLocus <- function(SNPdataset) { 
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  nuSites <- GetNumberOfSitesForLocus(SNPdataset)
  singleSNPfromLocus <- matrix(nrow=dim(SNPdataset)[1], ncol=length(nuSites))
  rownames(singleSNPfromLocus) <- rownames(SNPdataset)
  whichRandomSites <- NULL
  randomOrMax <- NULL
  for(locus in sequence(length(nuSites))){
    singleLocus <- as.matrix(SNPdataset[,locus])
    keepSNPs <- which(CalculateMissingData(singleLocus, "sites") == min(CalculateMissingData(singleLocus, "sites")))
    if(length(keepSNPs) > 1) {
      randomSite <- unname(sample(keepSNPs, 1)) # Modified to fix error
      whichRandomSites <- c(whichRandomSites, randomSite)
      randomOrMax <- c(randomOrMax, "R")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, randomSite]    
    }
    else {
      whichRandomSites <- c(whichRandomSites, keepSNPs)
      randomOrMax <- c(randomOrMax, "M")
      singleSNPfromLocus[,locus] <- SplitSNP(singleLocus)[, keepSNPs]
    }
  }
  singleSNPfromLocus <- cSNP(singleSNPfromLocus)
  if(snpclass == "snp")
    return(list(snpdata=ReadSNP(singleSNPfromLocus), position=whichRandomSites, randomOrMax=randomOrMax))
  return(list(snpdata=singleSNPfromLocus, position=whichRandomSites, randomOrMax=randomOrMax))
}

## ---- Remove parsimony-uninformative sites ----
RemoveUninformative <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  dropSingletons <- function(tab) {
    bad <- c("?") # for now, limited to one type of missing data code ("?")
    if (any(rownames(tab) == bad)) {
      tab <- tab[-which(rownames(tab)=="?")]
    }
    if (length(tab) == 2 && any(tab == "1")) {
      tab <- tab[-which(tab == "1")]
    }
    return(tab)
  }
  ## After dropping singletons (and missing data), if length > 1, site is informative; if site had only singleton, Returns NULL
  Informative <- function(tab) {
    if (length(tab) > 1) {
      return(tab)
    } else if  (any(tab==ntax)) {
      return(tab)
    }
  }
  ## Credit: https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists/26540063
  ## A helper function that tests whether an object is either NULL _or_ 
  ## a list of NULLs
  is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
  ## Recursively step down into list, removing all such objects 
  rmNullObs <- function(x) {
    x <- Filter(Negate(is.NullOb), x)
    lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
  }
  #  snps <- sum(nchar(SNPdataset[1,]))
  splitdata <- SplitSNP(SNPdataset)
  counts <- lapply(splitdata, table)
  dropMiss <- lapply(counts, dropSingletons)
  dropSin <- lapply(dropMiss, Informative)
  InformativeSites <- rmNullObs(dropSin)
  subset <- names(InformativeSites)
  bools <- colnames(splitdata) %in% subset
  newSNPdataset <- cSNP(splitdata, KeepVector = bools)
  #  newsnps <- sum(nchar(newSNPdataset[1,]))
  #  if(chatty)
  #    message(paste("removed", snps-newsnps, "of", snps, "sites"))
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}
