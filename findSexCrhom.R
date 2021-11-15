# identify the sex chromosome from depth data
library(data.table)
library(scales)
dp<-as.matrix(fread("depthKnulli.txt",header=FALSE))

## snp and ind. ids
load("../redwood_gwa/gp.rdat")
snps<-as.matrix(read.table("../redwood_gwa/Snps.txt",sep=":",header=FALSE))
Mal<-which(ph$Sex.at.Day.15=="MALE")
Fem<-which(ph$Sex.at.Day.15=="FEMALE")

Mdep<-apply(dp[,Mal],1,mean)
Fdep<-apply(dp[,Fem],1,mean)

ex<-1+(snps[,1]==775) ## == LG 
cs<-alpha(c("darkgray","firebrick"),c(.5,.7))
pdf("SexChrom.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
plot(Mdep,Fdep,pch=19,col=cs[sex],xlab="Coverage in males",ylab="Coverage in females",cex.lab=1.5)
abline(a=0,b=1)
legend(1,16,c("Chrom. 13","Other chrom."),rev(cs),bty='n')
dev.off()
