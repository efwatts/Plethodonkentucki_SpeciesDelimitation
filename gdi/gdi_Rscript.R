#set color palette
green <- "#488A5A"
blue <- "#7DBEDA"
plum <- "#8F5582"
tan <- "#E6DDB0"
darkblue <- "#053e4c"
brown <- "#471212"
black <- "#000000"



#set working directory
setwd("~/Desktop/BPP/gdi/A00_runs")

#import MCMC file
kfiles <- list.files(pattern="kentucki_A00.mcmc[1-5].txt")
ktip <- do.call(`rbind`,lapply(kfiles, read.table, header=T))


#generate figure GDI
densityKENTUCKY <-density(1-exp(-2*ktip$tau_8KENTUCKYCUMBERLAND/ktip$theta_1KENTUCKY))
densityCUMBERLAND <-density(1-exp(-2*ktip$tau_8KENTUCKYCUMBERLAND/ktip$theta_2CUMBERLAND))
densityNEW <- density(1-exp(-2*ktip$tau_9NEWKANAWHA/ktip$theta_3NEW))
densityKANAWHA <- density(1-exp(-2*ktip$tau_9NEWKANAWHA/ktip$theta_4KANAWHA))
densityPETRAEUS <- density(1-exp(-2*ktip$tau_6KENTUCKYCUMBERLANDNEWKANAWHAPETRAEUS/ktip$theta_5PETRAEUS))
densityKENTUCKYCUMBERLAND <- density(1-exp(-2*ktip$tau_7KENTUCKYCUMBERLANDNEWKANAWHA/ktip$theta_8KENTUCKYCUMBERLAND))
densityNEWKANAWHA <- density(1-exp(-2*ktip$tau_7KENTUCKYCUMBERLANDNEWKANAWHA/ktip$theta_9NEWKANAWHA))

plot((densityKENTUCKY), xlim=c(0,1), ylim=c(0,100), col=green, lwd=2, xlab="gdi", main="Plethodon kentucki")
lines(densityCUMBERLAND, lw=2, col=blue)
lines(densityNEW, lw=2, col=plum)
lines(densityKANAWHA, lw=2, col=tan)
lines(densityPETRAEUS, lw=2, col=black)
lines(densityKENTUCKYCUMBERLAND, lw=2, col=darkblue)
lines(densityNEWKANAWHA, lw=2, col = brown)
abline(v=0.2,col="grey50", lty=2)
abline(v=0.7, col="grey50", lty=2)
legend("topright", 10,legend=c("KENTUCKY", "CUMBERLAND", "NEW", "KANAWHA", "PETRAEUS","KENTUCKYCUMBERLAND", "NEWKANAWHA"),
       col=c(green,blue,plum,tan,black,darkblue, brown), pch=15, pt.cex=2, bty="n")


#import Sarracenia MCMC file
sctfiles <- list.files(pattern="mcmc[1-4]tip.txt")
Scintip <- do.call(`rbind`,lapply(sctfiles, read.table, header=T))

scmfiles <- list.files(pattern="mcmc[1-4]mid.txt")
Scinmid <- do.call(`rbind`,lapply(scmfiles, read.table, header=T))

#generate Scincella tip figures
densityWestern <-density(1-exp(-2*Scintip$tau_7WesternCentral/Scintip$theta_1Western))
densityCentral <-density(1-exp(-2*Scintip$tau_7WesternCentral/Scintip$theta_2Central))
densityWC <-density(1-exp(-2*Scinmid$tau_5WesternCentralEastern/Scinmid$theta_1WesternCentral))
densityE <-density(1-exp(-2*Scinmid$tau_5WesternCentralEastern/Scinmid$theta_2Eastern))
plot((densityWestern), xlim=c(0,1), ylim=c(0,40), col=dgr, lwd=2, main="Scincella", xlab="gdi")
lines(densityCentral, lw=2, col=mgr)
lines(densityWC, lw=2, col=lbr)
lines(densityE, lw=2, col=dbr)
abline(v=0.2,col="grey50", lty=2)
abline(v=0.7, col="grey50", lty=2)
legend("topright", 10,legend=c("Western", "Central", "Western+Central", "Eastern"),
       col=c(dgr,mgr,lbr,dbr), pch=15, pt.cex=2, bty="n")