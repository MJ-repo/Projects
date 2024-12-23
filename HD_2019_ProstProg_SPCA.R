##############################################################################################
# The prostat dataset from Efron 2010 and Hastie and Efron 2016                             ##
# The first 50 columns are the genetic activity measurements for the 50 control subjects,   ##
# while the last 52 columns represent the prostate cancer subjects.                         ##
##############################################################################################
rm(list=ls())
library(BVS)
library(Biobase)
library(multtest, verbose = FALSE)
par(mfrow=c(1,1))

prostmat <- read.csv("C:\\Users\\juach\\Documents\\UHasselt\\Analysis of High Dimensional Data\\Project\\prostmat.csv")
dim(prostmat)
#View(prostmat)
is.data.frame(prostmat)

x1<-prostmat[610,][1:50]
x2<-prostmat[610,][51:102]
t.test(x1,x2,var.equal=TRUE)
xx<-as.vector(prostmat[610,])
mean(xx)
p.transpose<-as.data.frame(t(prostmat))
(mean(p.transpose[,610]))

##########################################################
## Logistic regression model :gene by gene regression    #
##########################################################

# treatment group
classlabel<-c(rep(0,50),rep(1,52))
length(classlabel)
sum(classlabel)

length(prostmat[1,])
length(classlabel)
classlabel[1:50]
classlabel[51:102]

data1<-prostmat
ngenes<-length(data1[,1])

t.t<-p.val<-c(1:ngenes)
i=610
## t test ##
for(i in 1:ngenes)
{	
  xi<-as.numeric(data1[i,])
  modi<-glm(classlabel~xi,binomial(link = "logit"))
  t.t[i]<-summary(modi)$coefficients[2,3]
  p.val[i]<-summary(modi)$coefficients[2,4]
}
hist(p.val,nclass=50)
sum(p.val<0.05) #441
sum(p.val<0.1)  #772
sum(abs(t.t)>3.9)
plot(t.t,p.val)
abline(0.05,0,col=2)

index<-c(1:ngenes)
index1<-index[order(p.val)]
p.val.sort<-sort(p.val)[1]
kk<-index1[1]
xk<-as.numeric(data1[kk,])
modk<-glm(classlabel~xk,binomial(link = "logit"))
summary(modk)
modk$coefficients
plot(xk,classlabel, xlab="Gene Expression")
lines(sort(xk),modk$fit[order(xk)],lwd=3,col=2)
title(paste("Gene",kk))
xtable(summary(modk))
min(0.00000589,0.9/6033)
min(p.val)
sum(p.val<0.1/6033)


xx<-seq(from=-1,to=1,length=1000)
pred.x<-exp(modk$coefficients[1]+modk$coefficients[2]*xx)/(1+exp(modk$coefficients[1]+modk$coefficients[2]*xx))
plot(xx,pred.x,type="l",xlab="Gene Expression",ylab="P(Cancer)",ylim=c(0,1))
points(xk,classlabel)
abline(h=0.5,col=2)


########## Multiplicity ##########
library(Biobase)
library(multtest, verbose = FALSE)
rawp<-p.val
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),leg=c(4600,0.45),cex=0.8,lty=1,col=1:5,lwd=2)
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r
abline(0.05,0,col=1,lty=2)
abline(0.1,0,col=1,lty=2)



###################continue here###################

par(mfrow(c=1,1))
p.val.beta1<-p.val
t.val.beta1<-t.t
sort(p.val.beta1)[1:6]
sum(p.val.beta1 < 7.8478312e-05)
min(p.val.beta1)

index<-c(1:ngenes)
index1<-index[p.val.beta1 < 7.8478312e-05]
index1
plot(index1,p.val.beta1[p.val.beta1 <7.8478312e-05],type="l",ylab="p-value")


############################################################################################
##                                                                                        ##
##                                                                                        ##
## SPCA Top 5                                                                             ##
##                                                                                        ##
##                                                                                        ##
############################################################################################
K<-5
K<-10
cutoff<-sort(p.val.beta1)[K+1]
sum(p.val.beta1 < cutoff)
min(p.val.beta1)
index<-c(1:ngenes)
index1<-index[p.val.beta1 < cutoff]
index1

############################################################################################
###    the reduced matrix                                                                 ##
############################################################################################

data.R<-data1[index1,]
dim(data.R)

############################################################################################
###    PCA                                                                                ##
############################################################################################


pc.cr <- princomp(t(data.R),scores=TRUE,cor = TRUE)
loadings(pc.cr)
plot(pc.cr)
biplot(pc.cr)
summary(pc.cr)

library(xtable)
xtable(summary(pc.cr))

library(factoextra)
res.pca <- prcomp(t(data.R), scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)



############################################################################################
###    First PC                                                                           ##
############################################################################################

Ux1<-as.vector(pc.cr$scores[,1])
Ux2<-as.vector(pc.cr$scores[,2])
Ux3<-as.vector(pc.cr$scores[,3])
Ux4<-as.vector(pc.cr$scores[,4])
Ux5<-as.vector(pc.cr$scores[,5])
Ux6<-as.vector(pc.cr$scores[,6])
Ux7<-as.vector(pc.cr$scores[,7])
Ux8<-as.vector(pc.cr$scores[,8])
Ux9<-as.vector(pc.cr$scores[,9])
Ux10<-as.vector(pc.cr$scores[,10])
Ux11<-as.vector(pc.cr$scores[,11])
Ux12<-as.vector(pc.cr$scores[,12])
plot(Ux1,Ux2)
plot(Ux1,Ux2, col=c(rep(1,50), rep(2,52)))
############################################################################################
###    Plot: PC1 and response                                                             ##
############################################################################################


plot(Ux1,classlabel)
spca.mod <- glm(classlabel~Ux1,binomial(link = "logit"))
summary(spca.mod)
xtable(summary(spca.mod))

xx<-seq(from=-3,to=3,length=1000)
pred.x<-exp( spca.mod$coefficients[1] +spca.mod$coefficients[2]*xx)/(1+exp(spca.mod$coefficients[1] +spca.mod$coefficients[2]*xx))
plot(xx,pred.x,type="l",xlab="U(X)",ylab="P(Cancer)", col=2)
abline(h=0.5)
points(Ux1,classlabel)


plot(Ux1,classlabel)
summary(glm(classlabel~Ux1+Ux2,binomial(link = "logit")))
#Ux2 not significant

plot(Ux,classlabel)
summary(glm(classlabel~Ux1+Ux2+Ux3+Ux4,binomial(link = "logit")))
summary(glm(classlabel~Ux1+Ux2+Ux3+Ux4+Ux5+Ux6,binomial(link = "logit")))

############################################################################################
###   SPCA                                                                       ##
############################################################################################

data.Rtop5 <- t(data.R)
data.RUX <-as.data.frame(cbind(classlabel,Ux1,Ux2,Ux3,Ux4,Ux5))
data.RUX$Group[data.RUX$classlabel==1]="Cancer"
data.RUX$Group[data.RUX$classlabel==0]="Control"

library(MASS)
lda.fit = lda(Group~Ux1+Ux2, data=data.RUX)
lda.fit
plot(lda.fit, type = "both")

#Prior probabilities of groups:
#Cancer   Control 
#0.5098039 0.4901961 

#Group means:
#  Ux         Ux2
#Cancer  -0.9664451 -0.02350935
#Control  1.0051029  0.02444972

#Coefficients of linear discriminants:
#  LD1
#Ux  1.10221740
#Ux2 0.05027985
n = length(classlabel)
lda.pred = predict(lda.fit,data.RUX)
lda.class = lda.pred$class
table(lda.class,data.RUX$Group)
table(lda.class,data.RUX$Group)[1,1]
table(lda.class,data.RUX$Group)[1,2]
table(lda.class,data.RUX$Group)[2,1]
table(lda.class,data.RUX$Group)[2,2]
#MCE
lda.mce = (table(lda.class,data.RUX$Group)[1,2]+ table(lda.class,data.RUX$Group)[2,1])/n
lda.specificity = (table(lda.class,data.RUX$Group)[2,2])/((table(lda.class,data.RUX$Group)[2,2])+(table(lda.class,data.RUX$Group)[1,2]))
lda.sensitivity = (table(lda.class,data.RUX$Group)[1,1])/((table(lda.class,data.RUX$Group)[1,1])+(table(lda.class,data.RUX$Group)[2,1]))

#ROC
install.packages("pROC")
library(pROC)
classlabel
labels <- (ifelse(lda.pred$class=="Cancer", 1, 0))
roc.info<-roc(classlabel, labels,plot=T,legacy.axes=T)
roc(classlabel, labels,plot=T,legacy.axes=T, print.auc=T,col="#377eb8", percent=T)

roc.df <-data.frame(tpp = roc.info$sensitivities*100,
                    fpp = roc.info$specificities*100,
                    thresholds = roc.info$thresholds)
head(roc.df)

#Cross Validation
lda.cv <- lda(Group~Ux1+Ux2, data=data.RUX, CV=TRUE)
t(table(data.RUX$Group, lda.cv$class, dnn = c('Actual Group','Predicted Group')))


###########################################
n <- dim(data.RUX)[1]
cv.prediction <- c()

for (i in 1:n) {
  holdout <- data.RUX[-i,]
  
  holdout1 <- holdout[holdout$Group == "Cancer",][,2:6]
  holdout2 <- holdout[holdout$Group == "Control",][,2:6]
  
  holdout1.means <- apply(holdout1, 2, mean)
  holdout2.means <- apply(holdout2, 2, mean)
  
  n1 <- nrow(holdout1)
  n2 <- nrow(holdout2)
  
  w1 <- (n1 - 1) * var(holdout1)
  w2 <- (n2 - 1) * var(holdout2)
  
  sp1 <- 1 / (n1 + n2 - 2) * (w1 + w2)
  
  cutoff <- .5 * (holdout1.means - holdout2.means) %*% solve(sp1) %*% (holdout1.means + holdout2.means)
  
  ay <- (holdout1.means - holdout2.means) %*% solve(sp1) %*% as.numeric(data.RUX[i,2:6])
  group <- ifelse(ay > cutoff, 1, 0)
  cv.prediction <- append(cv.prediction, group)
}
t(table(data.RUX$Group, cv.prediction, dnn = c('Actual Group','Predicted Group')))

#############################################################

install.packages("superpc")
library(superpc)


x<-matrix(rnorm(1000*100),ncol=100)
v1<- svd(x[1:80,])$v[,1]

y<-2+5*v1+ .05*rnorm(100)

xtest<-x
ytest<-2+5*v1+ .05*rnorm(100)
censoring.status<- sample(c(rep(1,80),rep(0,20)))
censoring.status.test<- sample(c(rep(1,80),rep(0,20)))

featurenames <- paste("feature",as.character(1:1000),sep="")

data<-list(x=x,y=y, censoring.status= censoring.status, featurenames=featurenames)
data.test<-list(x=xtest,y=ytest, censoring.status=censoring.status.test, featurenames= featurenames)

train.obj<- superpc.train(data, type="survival")

cv.obj<-superpc.cv(train.obj, data)

superpc.plotcv(cv.obj)

lrtest.obj<-superpc.lrtest.curv(train.obj, data,data.test)

superpc.plot.lrtest(lrtest.obj)


fit.cts<- superpc.predict(train.obj, data, data.test, threshold=0.7, n.components=3, prediction.type="continuous")
fit.groups<- superpc.predict(train.obj, data, data.test, threshold=0.7, n.components=1, prediction.type="discrete")

superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.pred)

plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred), col=2:3, xlab="time", ylab="Prob survival")

fit.red<- superpc.predict.red(train.obj, data, data.test, threshold=0.7)

fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj,  data,  threshold=0.7)

superpc.plotred.lrtest(fit.redcv)
