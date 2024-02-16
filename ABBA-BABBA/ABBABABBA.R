setwd("G:/My Drive/OUwork/Research/ABBABABBA_anchored phylogenomics")
setwd("/Volumes/GoogleDrive-114548620038544830238/My Drive/OUwork/Research/ABBABABBA_anchored phylogenomics")

library('evobiR')

##the alignment should have three species and one outgroup, so i combined Kentucky/Cumberland
#species 1 = kanawha
#species 2 = cumberland/kentucky
#species 3 = new
#outgroup = petraeus
#(((K,C),N),P)
CalcD(alignment = "ABBABABA_input.fasta", sig.test = "J", block.size = 100, replicate = 100)
#chose "sig.test = "J" because jackknifing is more appropriate for alignment data ("B" for bootstrap would be good for ddRAD) (other option is "N" for none)
#block.size is the number of sites to be dropped for jackknifing
#replicate is number of replicates used in estimating variance

setwd("/Volumes/GoogleDrive-114548620038544830238/My Drive/OUwork/Research/ABBABABBA_anchored phylogenomics")
##reading about ABBA-BABA, I should actually test a few hypotheses (shows introgression between the )
CalcD(alignment = "((C,KY),K).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((C,KY),N).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((K,N),C).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((K,N),KY).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((KY,C),K).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((KY,C),N).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((N,K),C).fasta", sig.test = "J", block.size = 10, replicate = 100)

CalcD(alignment = "((N,K),KY).fasta", sig.test = "J", block.size = 10, replicate = 100)
 

#the above was one SNP per locus, so I am going to do all SNPs to see if that gives a different result...here, the scores are weirdly high
setwd("/Volumes/GoogleDrive-114548620038544830238/My Drive/OUwork/Research/ABBABABBA_anchored phylogenomics/all_SNPs")

library('evobiR')

CalcD(alignment = "((C,KY),K).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((C,KY),N).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((K,N),C).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((K,N),KY).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((KY,C),K).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((KY,C),N).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((N,K),C).fasta", sig.test = "J", block.size = 25, replicate = 100)

CalcD(alignment = "((N,K),KY).fasta", sig.test = "J", block.size = 25, replicate = 100)





