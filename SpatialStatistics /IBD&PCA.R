
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")
library("ade4")
library("apex")
library("mmod")
library("poppr")

# Creating DNAbin objects
bin <- fasta2DNAbin("0_missing_allSnps_filtered_nooutgroup.fasta")
head(bin)

k_data <- DNAbin2genind(bin)
k_data

k_data$tab[10:15, 12:15]

k_data@pop #empty, so adding pop numbers
# pop number key: "kanawha" = 1, "cumberland" = 2, "new" = 3, "kentucky" = 4
k_data@pop <- as.factor(c(1,1,2,4,1,4,2,3,3,4,2,2,2,4,1,2,4,3,3,3,4))

k_data@pop <- as.factor(c("kanawha","kanawha","cumberland","kentucky","kanawha","kentucky","cumberland","new","new","kentucky","cumberland","cumberland","cumberland","kentucky","kanawha","cumberland","kentucky","new","new","new","kentucky"))
k_data@pop #now pops are listed



##Using summaries
to <- summary(k_data)
to
names(to)

par(mfrow=c(1,1))
plot(to$n.by.pop, to$pop.n.all, xlab="Populations sample size",
     ylab="Number of alleles", main="Alleles numbers and sample sizes",
     type="n") 
text(to$n.by.pop, to$pop.n.all, lab=names(to$n.by.pop))

barplot(to$loc.n.all, ylab="Number of alleles", main="Number of alleles per locus")
barplot(to$Hexp-to$Hobs, ylab ="Heterozygosity: expected-observed",
        main = "Number of alleles per locus") #this doesn't work--no Hexp?
barplot(to$n.by.pop, main = "Sample sizes per population",
        ylab = "Number of genotypes", las=3)


##those aren't working, so trying some in hierfstat
bs <- basic.stats(k_data,diploid=TRUE)
head(bs)

#classical genetic distances estimation
genet.dist(k_data, diploid=TRUE, method="Dch") #Dch here 
genet.dist(k_data, diploid=TRUE, method="Fst") #Latter's Fst here ********THIS HAS POPULATION FST**************
genet.dist(k_data, diploid=TRUE, method="Da")


####Performing a PCA on genind objects
sum(is.na(k_data$tab)) #shows amount of missing data

X <- scaleGen(k_data, NA.method="mean") #scales to account for missing data
class(X)
dim(X)
X[1:5,1:5]

pca1 <- dudi.pca(X,cent=FALSE, scale=TRUE, scannf=FALSE, nf=3) #here, disabled scaling from before
barplot(pca1$eig[1:50], main="PCA eigenvalues", col=heat.colors(50))
pca1

s.label(pca1$li)
title("PCA of kentucki dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

#improve figure; set ellipses for populations
s.class(pca1$li, pop(k_data), grid=FALSE)
title("PCA of kentucki dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

#add third axis
s.class(pca1$li, pop(k_data), xax=1,yax=3,sub="PCA 1-3", csub=2, grid=FALSE)
title("PCA of kentucki dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20], nf=3, xax=1, yax=3)

pca2 <- dudi.pca(X,cent=FALSE, scale=FALSE, scannf=FALSE, nf=3) #here, disabled scaling from before
barplot(pca2$eig[1:50], main="PCA eigenvalues", col=heat.colors(50))
pca2

s.label(pca2$li)
title("PCA of kentucki dataset\naxes 1-2")
add.scatter.eig(pca2$eig[1:20], 3,1,2)

#improve figure; set ellipses for populations
s.class(pca2$li, pop(k_data), grid=FALSE)
title("PCA of kentucki dataset\naxes 1-2")
add.scatter.eig(pca2$eig[1:20], 3,1,2)

#add third axis
s.class(pca1$li, pop(k_data), xax=1,yax=3,sub="PCA 1-3", csub=2, grid=FALSE, addaxes = TRUE)
title("PCA of kentucki dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20], nf=3, xax=1, yax=3)

percent_variance <- (pca2$eig / sum(pca2$eig)) * 100
percent_variance[1]

percent_variance1 <- (pca1$eig / sum(pca1$eig)) * 100
percent_variance1

##############################################################################
######################Spatial analysis########################################
#############################################################################

##Isolation by distance
bin <- fasta2DNAbin("0_missing_allSnps_filtered_nooutgroup.fasta")
head(bin)

k_ind <- DNAbin2genind(bin)
k_ind

k_ind$tab[10:15, 12:15]

k_ind@pop #empty, so adding pop numbers
# treating each ind as a pop
k_ind@pop <- as.factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))
k_ind@pop #now pops are listed

library(readr)
Coord<-read_csv("Locations.csv")
head(Coord)
colnames(Coord)<-c("X", "Y")
k_data@other$xy<-Coord
k_data@other$xy

k_pop <- genind2genpop(k_data)
obj <- seppop(k_data)
new <- obj$new
cumberland <- obj$cumberland
kentucky <- obj$kentucky
kanawha <- obj$kanawha

Dgen <- dist.genpop(k_pop, method=1) #Neis
#Dgen <- dist(k_data, method = "euclidean", diag=FALSE, upper=FALSE, p=2) #this is euclidean...want Fst
Dgeo <- dist(other(k_pop)$xy) #Euclidean is default here
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
plot(ibd) #IBD significant


#cline or distant patches?
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo))
#scatterplot to view data\

#2-dimensional kernerl density estimation 
library(MASS)
dens <- kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20, cex=1)
image(dens, col=transp(myPal(300), .7), add=TRUE)
abline(lm(Dgen~Dgeo))

#plot shows single consistent cloud witout discontinuities, which would have indicated patached
#indicative of IBD


#as is above k_pop is separated into 4 pops, but k_ind shows all 21 as individuals

#############################################################################################################
##################START HERE FOR NEW FIGURES (once data is loaded in)########################################
#############################################################################################################

k_ind@other$xy<-Coord
k_ind@other$xy

k_pop <- genind2genpop(k_ind)

Dgen <- dist.genpop(k_pop, method=1) 
#Dgen <- dist(k_data, method = "euclidean", diag=FALSE, upper=FALSE, p=2) #this is euclidean...want Fst
Dgeo <- dist(other(k_pop)$xy) #Euclidean is default here
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
plot(ibd) #IBD significant


#cline or distant patches?
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo))
#scatterplot to view data\

#2-dimensional kernerl density estimation 
library(MASS)
dens <- kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20, cex=1)
image(dens, col=transp(myPal(300), .7), add=TRUE)
abline(lm(Dgen~Dgeo))


######IBD for new#####

new$pop <- as.factor(c(1,2,3,4,5))

new_pop <- genind2genpop(new)
new_pop

new_Dgen <- dist.genpop(new_pop, method=1) 
new_Dgeo <- dist(other(new_pop)$xy) #Euclidean is default here
new_ibd <- mantel.randtest(new_Dgen,new_Dgeo)
new_ibd
plot(new_ibd) #IBD significant


#cline or distant patches?
plot(new_Dgeo, new_Dgen)
abline(lm(new_Dgen~new_Dgeo))
#scatterplot to view data

#2-dimensional kernerl density estimation 
library(MASS)
new_dens <- kde2d(new_Dgeo, new_Dgen, n=300)
new_myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(new_Dgeo, new_Dgen, pch=20, cex=1)
image(new_dens, col=transp(new_myPal(300), .7), add=TRUE)
abline(lm(new_Dgen~new_Dgeo))


######IBD for cumberland####
cumberland$pop <- as.factor(c(1,2,3,4,5,6))

cumberland_pop <- genind2genpop(cumberland)
cumberland

cumberland_Dgen <- dist.genpop(cumberland_pop, method=1) 
#Dgen <- dist(k_data, method = "euclidean", diag=FALSE, upper=FALSE, p=2) #this is euclidean...want Fst
cumberland_Dgeo <- dist(other(cumberland_pop)$xy) #Euclidean is default here
cumberland_ibd <- mantel.randtest(cumberland_Dgen,cumberland_Dgeo)
cumberland_ibd
plot(cumberland_ibd) #IBD significant


#cline or distant patches?
plot(cumberland_Dgeo, cumberland_Dgen)
abline(lm(cumberland_Dgen~cumberland_Dgeo))
#scatterplot to view data

#2-dimensional kernerl density estimation 
library(MASS)
cumberland_dens <- kde2d(cumberland_Dgeo, cumberland_Dgen, n=300)
cumberland_myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(cumberland_Dgeo, cumberland_Dgen, pch=20, cex=1)
image(cumberland_dens, col=transp(cumberland_myPal(300), .7), add=TRUE)
abline(lm(cumberland_Dgen~cumberland_Dgeo))


#####IBD for kentucky####
kentucky$pop <- as.factor(c(1,2,3,4,5,6))

kentucky_pop <- genind2genpop(kentucky)
kentucky

kentucky_Dgen <- dist.genpop(kentucky_pop, method=1) 
#Dgen <- dist(k_data, method = "euclidean", diag=FALSE, upper=FALSE, p=2) #this is euclidean...want Fst
kentucky_Dgeo <- dist(other(kentucky_pop)$xy) #Euclidean is default here
kentucky_ibd <- mantel.randtest(kentucky_Dgen,kentucky_Dgeo)
kentucky_ibd
plot(kentucky_ibd) #IBD significant


#cline or distant patches?
plot(kentucky_Dgeo, kentucky_Dgen)
abline(lm(kentucky_Dgen~kentucky_Dgeo))
#scatterplot to view data

#2-dimensional kernerl density estimation 
library(MASS)
kentucky_dens <- kde2d(kentucky_Dgeo, kentucky_Dgen, n=300)
kentucky_myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(kentucky_Dgeo, kentucky_Dgen, pch=20, cex=1)
image(kentucky_dens, col=transp(kentucky_myPal(300), .7), add=TRUE)
abline(lm(kentucky_Dgen~kentucky_Dgeo))



####IBD for kanawha####
kanawha$pop <- as.factor(c(1,2,3,4))

kanawha_pop <- genind2genpop(kanawha)
kanawha

kanawha_Dgen <- dist.genpop(kanawha_pop, method=1) 
#Dgen <- dist(k_data, method = "euclidean", diag=FALSE, upper=FALSE, p=2) #this is euclidean...want Fst
kanawha_Dgeo <- dist(other(kanawha_pop)$xy) #Euclidean is default here
kanawha_ibd <- mantel.randtest(kanawha_Dgen,kanawha_Dgeo)
kanawha_ibd
plot(kanawha_ibd) #IBD significant


#cline or distant patches?
plot(kanawha_Dgeo, kanawha_Dgen)
abline(lm(kanawha_Dgen~kanawha_Dgeo))
#scatterplot to view data

#2-dimensional kernerl density estimation 
library(MASS)
kanawha_dens <- kde2d(kanawha_Dgeo, kanawha_Dgen, n=300)
kanawha_myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(kanawha_Dgeo, kanawha_Dgen, pch=20, cex=1)
image(kanawha_dens, col=transp(kanawha_myPal(300), .7), add=TRUE)
abline(lm(kanawha_Dgen~kanawha_Dgeo))

####making four panels######
par(mfrow=c(2,2))

plot(new_Dgeo, new_Dgen, pch=20, cex=1)
image(new_dens, col=transp(new_myPal(300), .7), add=TRUE)
abline(lm(new_Dgen~new_Dgeo))

plot(cumberland_Dgeo, cumberland_Dgen, pch=20, cex=1)
image(cumberland_dens, col=transp(cumberland_myPal(300), .7), add=TRUE)
abline(lm(cumberland_Dgen~cumberland_Dgeo))

plot(kentucky_Dgeo, kentucky_Dgen, pch=20, cex=1)
image(kentucky_dens, col=transp(kentucky_myPal(300), .7), add=TRUE)
abline(lm(kentucky_Dgen~kentucky_Dgeo))

plot(kanawha_Dgeo, kanawha_Dgen, pch=20, cex=1)
image(kanawha_dens, col=transp(kanawha_myPal(300), .7), add=TRUE)
abline(lm(kanawha_Dgen~kanawha_Dgeo))

par(mfrow=c(1,1))



