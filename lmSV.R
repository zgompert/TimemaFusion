## R script to test SV/fusion genotype effect on performance
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)

## pca RW without BCTURN
pc<-prcomp(t(g_RW_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted


## lm fits, RW only 15 d weight significant
summary(lm(gemma_phRW_sub[,1] ~ gen))

#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.019278   0.008424   2.289   0.0273 *
#gen         -0.008555   0.003821  -2.239   0.0306 *
#---
#Residual standard error: 0.01902 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.109,	Adjusted R-squared:  0.08722 
#F-statistic: 5.013 on 1 and 41 DF,  p-value: 0.03064

## pca C without BCTURN
pc<-prcomp(t(g_C_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted ... ran till they were

## lm fits, C 15 and 21 d weight and survival all significant
summary(lm(gemma_phC_sub[,1] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.04160    0.01269   3.279  0.00194 **
#gen         -0.01845    0.00572  -3.226  0.00226 **
#---
#Residual standard error: 0.03116 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1782,	Adjusted R-squared:  0.1611 
#F-statistic: 10.41 on 1 and 48 DF,  p-value: 0.002261

summary(lm(gemma_phC_sub[,2] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.047148   0.012501   3.772 0.000471 ***
#gen         -0.020302   0.005498  -3.693 0.000598 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02801 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.2326,	Adjusted R-squared:  0.2155 
#F-statistic: 13.64 on 1 and 45 DF,  p-value: 0.0005979

summary(lm(gemma_phC_sub[,3] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.59451    0.12126   4.903 1.04e-05 ***
#gen          0.14099    0.05519   2.555   0.0137 *  
#---
#Residual standard error: 0.3064 on 50 degrees of freedom
#Multiple R-squared:  0.1154,	Adjusted R-squared:  0.09775 
#F-statistic: 6.526 on 1 and 50 DF,  p-value: 0.01373
tapply(X=gemma_phC_sub[,3],INDEX=gen,mean)
#        1         2         3 
#0.7142857 0.9047619 1.0000000 
## out of 10, 19, and 17 inds.

## combined
g_comb_sub<-cbind(g_RW_sub[a,],g_C_sub[a,])
pc<-prcomp(t(g_comb_sub),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted
save(list=ls(),file="svLm.rdat")

N_rw<-dim(gemma_phRW_sub)[1]
N_c<-dim(gemma_phC_sub)[1]
host_trt<-rep(c(1,2),c(N_rw,N_c))
summary(lm(gemma_phRW_sub[,1] ~ gen[host_trt==1]))
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)     0.019278   0.008424   2.289   0.0273 *
#gen[host == 1] -0.008555   0.003821  -2.239   0.0306 *
#Residual standard error: 0.01902 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.109,	Adjusted R-squared:  0.08722 
#F-statistic: 5.013 on 1 and 41 DF,  p-value: 0.03064
summary(lm(gemma_phC_sub[,1] ~ gen[host_trt==2]))
#               Estimate Std. Error t value Pr(>|t|)   
#(Intercept)     0.04160    0.01269   3.279  0.00194 **
#gen[host == 2] -0.01845    0.00572  -3.226  0.00226 **
#---
#Residual standard error: 0.03116 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1782,	Adjusted R-squared:  0.1611 
#F-statistic: 10.41 on 1 and 48 DF,  p-value: 0.002261

summary(lm(gemma_phC_sub[,3] ~ gen[host_trt==2]))
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.59451    0.12126   4.903 1.04e-05 ***
#gen[host == 2]  0.14099    0.05519   2.555   0.0137 *  
#Residual standard error: 0.3064 on 50 degrees of freedom
#Multiple R-squared:  0.1154,	Adjusted R-squared:  0.09775 
#F-statistic: 6.526 on 1 and 50 DF,  p-value: 0.01373

##################################################
csC<-c("cadetblue2","cadetblue3","cadetblue4")
csRW<-c("firebrick2","firebrick3","firebrick4")
cl<-1.5;lx<-1.5;cm<-1.5;ca<-1.1

pdf("knulliPerformanceSV.pdf",width=8,height=12)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,1],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca)
mns<-tapply(X=gemma_phC_sub[,1],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.002",side=3,adj=.1,line=-2)
title("(a) 15d weight, C",cex.main=cm)
box()

plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,1],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca)
mns<-tapply(X=gemma_phRW_sub[,1],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.031",side=3,adj=.1,line=-2)
title("(b) 15d weight, RW",cex.main=cm)
box()

plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,2],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca)
mns<-tapply(X=gemma_phC_sub[,2],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P < 0.001",side=3,adj=.1,line=-2)
title("(c) 21d weight, C",cex.main=cm)
box()

plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,2],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca)
mns<-tapply(X=gemma_phRW_sub[,2],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.138",side=3,adj=.1,line=-2)
title("(d) 21d weight, RW",cex.main=cm)
box()

surv<-tapply(INDEX=gen[host_trt==2]-1,gemma_phC_sub[,3],mean)
barplot(surv,col=csC,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca)
mtext("P = 0.014",side=3,adj=.5,line=-2)
title("(e) Survival, C",cex.main=cm)

surv<-tapply(INDEX=gen[host_trt==1]-1,gemma_phRW_sub[,3],mean)
barplot(surv,col=csRW,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca)
mtext("P = 0.759",side=3,adj=.5,line=-2)
title("(f) Survival, RW",cex.main=cm)
dev.off()


###################################################
dat<-dat[1:138,]
host_source_RW<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli")])
host_source_C<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli")])
host_source<-c(host_source_RW,host_source_C)

dat_sub_C<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli"),]
dat_sub_RW<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli"),]


###################################################

tapply(X=gen==1,INDEX=host_source,mean)
#        C        RW 
#0.1052632 0.6800000 
tapply(X=gen==2,INDEX=host_source,mean)
#        C        RW 
#0.4605263 0.3200000 
tapply(X=gen==3,INDEX=host_source,mean)
#        C        RW 
#0.4342105 0.0000000 

