# Code for "Great tits differ in glucocorticoid plasticity in response to spring temperature"
# BioRxiv (doi: 10.1101/2022.04.21.489013)
# Hau M, Deimel C, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper.

######################################################
# DATA ANALYSIS OF REACTION NORMS
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
data$Temp=as.numeric(scale(data$Temperature_at_capture)) #continuous variable

#Model with residuals divided into 5 blocks (1 per year)
nblocks5 <- 5	# number of 'residual blocks'
data$YEAR=as.numeric(data$Year)
data$envclass5 <- as.numeric(arules::discretize(data$YEAR,breaks= nblocks5, method='interval'))
data$envclass5=as.factor(data$envclass5)

#set prior
prior<- list(R = list(V = diag(5), nu = 2.002),
               G = list(G1 = list(V = diag(2), nu =2.002,alpha.mu = rep(0, 2), alpha.V= diag(25^2, 2, 2))))


# To set number of iterations
a=1000

# To run model
mod<-MCMCglmm(C0~SEX+year+C0T+Temp,
                 random=~us(1+Temp):ID,
                 rcov= ~idh(envclass5):units,    
                 family="gaussian",data=dat,
                 nitt=13000*a,thin=10*a,burnin=3000*a,
                 prior=prior,verbose=TRUE, pr=TRUE)

#Fixed effects
posterior.mode(mod$Sol)
HPDinterval(mod$Sol)
#plot(mod$Sol)

#Random effects
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)
plot(mod$VCV)

# To calculate the intercept-slope correlation
int.slope.cor <- mod$VCV[, 2]/sqrt(mod$VCV[, 1] * mod$VCV[, 4])
posterior.mode(int.slope.cor)
HPDinterval(int.slope.cor)
plot(int.slope.cor)
