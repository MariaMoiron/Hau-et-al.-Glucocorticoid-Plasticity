# Code for "Great tits differ in glucocorticoid plasticity in response to spring temperature"
# BioRxiv (doi: 10.1101/2022.04.21.489013)
# Hau M, Deimel C, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper.

######################################################
# DATA ANALYSIS OF WITHIN-INDIVIDUAL PLASTICITY
######################################################

# Loading packages
library(tidyr)
library(dplyr)
library(lme4)
library(arm)
library(ggplot2)
library (broom)
library(MCMCglmm)
library(nlme)

# Loading phenotypic data
data- read.csv("model1.csv", header=TRUE, na.strings=c(""," ","NA"))

# Formating response variabe
data$C0=as.numeric(data$CORT0/1000) #choose either CORT0 or CORT30

# Formating random effects
data$ID=as.factor(data$Ring_number) #Individual effect

# Formating fixed effects
data$year=as.factor(data$Year)
data$C0T=as.numeric(scale(data$time_CORT0sec))
data$SEX=as.factor(data$sex-1) #two levels: male or female

# To estimate mean and delta Temperature values
#Extract mean individual values of temperature
meansTemp <- as.data.frame(tapply(data$Temp, data$ID, mean, rm.na=TRUE))
ID <- as.data.frame(rownames(meansTemp))
df1 <- cbind(ID, meansTemp)
names(df1) <- c("ID","IDmeanTemp")
df1$IDmeanTemp=as.numeric(df1$IDmeanTemp)

# To merge both dataframes (individual mean dataframe and original dataframe) and extract deviations from individual mean values of Temperature
dataMerged <- as.data.frame(merge(data, df1,by="ID"))
dataMerged$Tempsd <- as.numeric(dataMerged$Temp-dataMerged$IDmeanTemp)
data=dataMerged

# To set prior
prior<- list(R = list(V = 1, nu = 1.002),
              G = list(G1 = list(V = 1, nu = 1.002, alpha.mu = 0, alpha.V = 25^2)))
              
# To set number of iterations
a=10000           

# To run model
mod<-MCMCglmm(C0~SEX+year+C0T+Tempsd+IDmeanTemp,
                 random=~ID,
                 family="gaussian",data=dat,
                 nitt=1300*a,thin=1*a,burnin=300*a,
                 prior=prior,verbose=TRUE)
                 
modbi$DIC #extract DIC value
summary(mod)

# Assesing convergence and auto-correlation of posterior distributions 
plot(mod$VCV)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

#Fixed effects
posterior.mode(mod$Sol)
HPDinterval(mod$Sol)
#plot(mod$Sol)
 
difference=mod$Sol[,9]-mod$Sol[,8] #extract the difference between mean and delta Temperature
round(posterior.mode(difference),3)
round(HPDinterval(difference), 3)
plot(difference)

#Random effects
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)
plot(mod$VCV)

