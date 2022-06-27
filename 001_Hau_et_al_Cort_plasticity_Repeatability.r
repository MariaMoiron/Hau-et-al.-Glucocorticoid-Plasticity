# Code for "Great tits differ in glucocorticoid plasticity in response to spring temperature"
# BioRxiv (doi: 10.1101/2022.04.21.489013)
# Hau M, Deimel C, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper.

######################################################
# DATA ANALYSIS OF REPEATABILITY
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
data$Temp=as.numeric(scale(data$Temperature_at_capture)) #continuous variable
data$SEX=as.factor(data$sex-1) #two levels: male or female

# To set prior for Bayesian model
prior<- list(R = list(V = 1, nu = 1.002),
              G = list(G1 = list(V = 1, nu = 1.002,alpha.mu = 0,alpha.V = 25^2)))

# To set number of iterations
a=10000

# To run model
mod <-MCMCglmm(C0~SEX+year+C0T+Temp,
                 random=~ID,
                 family="gaussian",data=data,
                 nitt=1300*a,thin=1*a,burnin=300*a,
                 prior=prior,verbose=TRUE)

modbi$DIC #extract DIC value
summary(mod)

# Assesing convergence and auto-correlation of posterior distributions 
plot(mod$VCV)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

# To calculate repeatability
rpt=mod$VCV[,"ID"]/rowSums(mod$VCV)  
mean(rpt) ###mean
HPDinterval(rpt) ###Calculate 95%CI
plot(rpt)
