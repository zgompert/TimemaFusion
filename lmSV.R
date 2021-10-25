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

