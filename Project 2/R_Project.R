rm(list=ls())
getwd()
setwd("C:/Users/juach/Documents\\UHasselt\\Generalized Linear Models\\Project")

#############################################################
###################Chapter 1#################################
#############################################################

#####################
#   Baseline Logit  #
#####################

mhdori=read.table("EG.dat",sep=" ",header=T)
attach(mhdori)
normal=as.numeric(response==1) 
malformed=as.numeric(response==2) 
dead=as.numeric(response==3)
mhd1<-data.frame(normal,malformed,dead, dose) 
mhd2<-data.frame(responsen=factor(response), dose) 
detach(mhdori)
head(mhd1)
View(mhd1)
str(mhd2)

# baseline multinomial model 
# first option: function vglm (0-1 per category response)
library(VGAM) 

mnfit1=vglm(cbind(normal,malformed,dead)~dose,multinomial(refLevel=1),mhd1) 
summary(mnfit1) #same in sas 

# baseline multinomial model 
# second option: function multinom (summarized/1-2-3-4)
library(nnet) 
mnfit3b=multinom(responsen~dose,mhd2)
summary(mnfit3b)

mhd2b=data.frame(responsen=relevel(x=mhd2$responsen, ref="1"),mhd2) 
mnfit3c=multinom(responsen~dose,mhd2b)
summary(mnfit3c)

install.packages("effects")
library(effects) 
mnfit3d=allEffects(mnfit3c) 
plot(mnfit3d) 
plot(mnfit3d,style="stacked")


#####################
#   Adjacent Logit  #
#####################

# adjacent logits model 
# first option: function vglm 
alfit1=vglm(cbind(normal,malformed,dead)~dose,acat,mhd1) 
summary(alfit1)
alfit1p=vglm(cbind(normal,malformed,dead)~dose,acat(parallel=TRUE),mhd1) 
summary(alfit1p)


LRT.2p=2*(logLik(alfit1)-logLik(alfit1p))   #LIKELIHOOD RATIO TEST
c(LRT.2p,1-pchisq(LRT.2p,1)) #2 different slopes better #same with sas
#Look for Residual deviance for -2logL equal to SAS
#Look for Log-Likelihood to compute difference in 2models: 2*(LLsaturdated -LLourmodel)


#####################
# Continuation Ratio#
#####################

# continuation ratio model 
cr.fit1<-vglm(cbind(normal,malformed,dead)~dose,cratio(parallel=TRUE),mhd1) #more sensible #same slope
cr.fit2<-vglm(cbind(normal,malformed,dead)~dose,cratio,mhd1) #more sensible #different slope
LRT.2p=2*(logLik(cr.fit2)-logLik(cr.fit1))   #LIKELIHOOD RATIO TEST
c(LRT.2p,1-pchisq(LRT.2p,1)) #different slopes 
summary(cr.fit1)
summary(cr.fit2)

#####################
# Cumulative Links  #
#####################

# proportional odds model 
# first option: function vglm 
pofiti=vglm(cbind(normal,malformed,dead)~dose,propodds(reverse=F),mhd1) #same slope (same with cumulative parallel=FALSE)
pofit1=vglm(cbind(normal,malformed,dead)~dose,cumulative(reverse=F),mhd1) #different slope
pofit2=vglm(cbind(normal,malformed,dead)~dose,cumulative(reverse=F,parallel=TRUE),mhd1) #same slope
pofit1;summary(pofit1)
pofit2;summary(pofit2)

LRT.int=2*(logLik(pofit1)-logLik(pofit2))  #LIKELIHOOD RATIO TEST
c(LRT.int,1-pchisq(LRT.int,1)) #different slope better

# proportional odds model 
# second option: function polr 
library(MASS) 
pofit2b=polr(mental~life+ses,mhd2) 
summary(pofit2b) #Check AIC
pofit2c=allEffects(pofit2b)
plot(pofit2c) 
plot(pofit2c,style="stacked")

# test proportional odds assumption 
# other way to fit the PO model 
# pofit2=vglm(cbind(well,mild,moderate,impaired)~life+ses,cumulative(reverse=F,parallel=TRUE),mhd1) 
pofit2np=vglm(cbind(well,mild,moderate,impaired)~life+ses,cumulative(reverse=F),mhd1) 
summary(pofit2np)
LRT.po=2*(logLik(pofit2np)-logLik(pofit2))  
c(LRT.po,1-pchisq(LRT.po,4))

# other link function
summary(vglm(cbind(well,mild,moderate,impaired)~life+ses,cumulative(link="probit",parallel=TRUE),mhd1))


#############################################################
###################Chapter 2#################################
#############################################################

mhdori=read.table("EG.dat",sep=" ",header=T)
attach(mhdori)
normal=as.numeric(response==1) 
malformed=as.numeric(response==2) 
dead=as.numeric(response==3)
mhd1<-data.frame(normal,malformed,dead, dose) 
mhd2<-data.frame(responsen=factor(response), dose) 
detach(mhdori)
head(mhd1)
View(mhd1)
str(mhd2)
str(mhdori)
unique(mhdori$id)

mhdori2=read.csv("r_mice.csv",header=T)
str(mhdori2)   

mhdori3=read.csv("r_mice_orig.csv",header=T)
str(mhdori3)   



#####################
# Quasi Likelihood  #
#####################

# logistic regression ignoring clustering 
logitfit=glm(cbind(y,total-y)~dose, family=binomial(link="logit"),data=mhdori2) 
summary(logitfit)

# quasi-likelihood logistic regression (pearson)
quasifit=glm(cbind(y,total-y)~dose, family=quasibinomial(link="logit"),data=mhdori2) 
summary(quasifit)
#Residual Deviance = Deviance in SAS
#For Quasi Likelihood Williams, REFER TO SAS CODES


#####################
#       GEE         #
#####################

# generalized estimating equations 
install.packages("gee")
library(gee) 

library(gee) 
geefitind=gee(resp_normal~dose, family=binomial(link="logit"),id=id, corstr="independence",data=mhdori3) 
summary(geefitind) 
geefitexch=gee(resp_normal~dose, family=binomial(link="logit"),id=id, corstr="exchangeable",data=mhdori3) 
summary(geefitexch) #same with SAS


#####################
#  Full Likelihood  #
#####################


install.packages("aod")
library(aod) 
bbfit=betabin(cbind(y,total-y)~dose,~1,data=mhdori2) 
summary(bbfit) #FINAL FULL LIKELIHOOD RESULT


#####################
#  Random Effects   #
#####################
#FINAL RANDOM EFFECTS MODEL IN SAS!

# logistic-normal model 
library(lme4) #default is Laplacian approximation nAGQ=1 
glmmfit1 <- glmer(cbind(y,total-y)~dose + (1|id), family = binomial(link="logit"), data = mhdori2) 
summary(glmmfit1) 

# using Adaptive Gauss-Hermite Quadrature approximation to the log-likelihood 
glmmfit2 <- glmer(cbind(y,total-y)~dose + (1|id), family = binomial(link="logit"),nAGQ=500,data = mhdori2) 
summary(glmmfit2)

#alternative data structure 
gm1 <- glmer(cbind(y,total-y)~dose + (1|id), family = binomial(link="logit"), data = trmoore) 
summary(gm1)

# using penalized-quasi likelihood PQL approximation 
library(mgcv) 
glmmfit<-gamm(cbind(y,total-y)~dose + (1|id), family=binomial(link="logit"),random=list(litterid=~1),data=mhdori2)
summary(glmmfit$gam)


#############################################################
###################Chapter 3#################################
#############################################################

insomnia=read.csv("3_insomnia2.csv") 
head(insomnia) 

#####################
#       GEE         #
#####################

install.packages("multgee")
library(multgee) 
PO.GEE.fit=ordLORgee(response~dose, data=mhdori3,id=id,LORstr="independence") 
summary(PO.GEE.fit) 
PO.GEE.fit2=ordLORgee(response~dose, data=mhdori3,id=id,LORstr="uniform") 
summary(PO.GEE.fit2)
#FINAL MODEL IN SAS!


#####################
#  Random Effects   #
#####################
str(mhdori3)
install.packages("ordinal")
library(ordinal) 

PO.GLMM= clmm(as.factor(response)~dose+(1|id),data=mhdori3,link="logit",nAGQ=40) 
summary(PO.GLMM)
#FINAL MODEL IN SAS!
