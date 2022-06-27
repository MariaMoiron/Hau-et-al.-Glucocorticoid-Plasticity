# Code for "Great tits differ in glucocorticoid plasticity in response to spring temperature"
# BioRxiv (doi: 10.1101/2022.04.21.489013)
# Hau M, Deimel C, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper.

######################################################
# DATA ANALYSIS OF BEWTEEN CORT0 AND CORT30 COVARIANCE
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
data$C0=as.numeric(data$CORT0/1000)
datA$C30=as.numeric(dat$CORT30/1000)

# Formating random effects
data$ID=as.factor(data$Ring_number) #Individual effect

# Formating fixed effects
data$year=as.factor(data$Year)
data$C0T=as.numeric(scale(data$time_CORT0sec))
dat$C30T=as.numeric(scale(dat$timeCORT30))
data$SEX=as.factor(data$sex-1) #two levels: male or female
data$Temp=as.numeric(scale(data$Temperature_at_capture)) #continuous variable

# Setting prior
prior= list(R = list(V = diag(2), nu = 3.002),
                     G = list(G1 = list(V = diag(2), nu =3.002,
                                        alpha.mu = rep(0,2),
                                        alpha.V = diag(1000,2,2))))
                                        

# To set number of iterations
a=1000

# To set model
model<- MCMCglmm(cbind(C0,C30) ~ trait-1+ #fixed effects
                  trait:year+  #fixed effects
                  trait:SEX+
                  trait:Temp+
                  at.level(trait,1):C0T + 
                  at.level(trait,2):C30T, 
                  random=~ us(trait):ID,                        #random effect ID
                  rcov=~ us(trait):units,                   #residual variance, "us" fit unstructured residual covariance
                  family=c("gaussian","gaussian"), #distribution of data
                  prior=prior, #to specify priors
                  nitt=13000*a,thin=10*a,burnin=3000*a, #to specify the number of iterations
                  verbose=TRUE, #to plot the iterations on the screen when running
                  data=data, #load data
                  pr=TRUE #to save the random effect estimates for each ID and fit later to a plot
)

##Summary of model
summary(model)

#Extract DIC value of the model
model$DIC

#Fixed effects
posterior.mode(model$Sol)
HPDinterval(model$Sol)
#plot(model$Sol)

#Random effects
round(posterior.mode(model$VCV),3)
round(HPDinterval(model$VCV), 3)
plot(model$VCV)

###Correlation among individuals
posterior.mode(posterior.cor(model$VCV[,1:4]))
HPDinterval(posterior.cor(model$VCV[,1:4]))

###Correlation within individuals
posterior.mode(posterior.cor(model$VCV[,5:8]))
HPDinterval(posterior.cor(model$VCV[,5:8]))
