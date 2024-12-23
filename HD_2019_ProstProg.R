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

x1<-prostmat[1,][1:50]
x2<-prostmat[1,][51:102]
t.test(x1,x2,var.equal=TRUE)


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

t.t<-p.val<-beta.t<-c(1:ngenes)

## t test ##
for(i in 1:ngenes)
{	
  xi<-as.numeric(data1[i,])
  modi<-glm(classlabel~xi,binomial(link = "logit"))
  t.t[i]<-summary(modi)$coefficients[2,3]
  beta.t[i]<-summary(modi)$coefficients[2,1]
  p.val[i]<-summary(modi)$coefficients[2,4]
}
hist(p.val,nclass=50)
sum(p.val<0.05) #441
sum(p.val<0.1)  #772

plot(t.t,p.val)
abline(0.05,0,col=2)
par(mfrow=c(1,2))
plot(t.t,-log(p.val), xlab="Test Statistic", ylab="-log(Pvalue)")
plot(beta.t,-log(p.val), xlab="Beta", ylab="-log(Pvalue)")
par(mfrow=c(1,1))
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


############################
#Usmallest pvalues
par(mfrow=c(1,1))
index<-c(1:ngenes)
index1<-index[order(p.val)]
small.index <- index1[1:10] #610 1720 3940  364  914 4546  332 4331 1089 3647
small.p <-p.val[small.index]
#Boxplot Per Gene (top 10)
par(mfrow=c(2,5))
for(i in 1:10){
  kk<-index1[i]
  xk<-as.numeric(data1[kk,])
  boxplot(split(xk,classlabel),names=c("Control","Cancer"))
  title(paste("Gene",kk)) #Gene 610
}

#Logistic Model Per Gene (top 10)
par(mfrow=c(2,5))
for(i in 1:10){
  kk<-index1[i]
  xk<-as.numeric(data1[kk,])
  modk<-glm(classlabel~xk,binomial(link = "logit"))
  summary(modk)
  modk$coefficients
  plot(xk,classlabel,xlab="Gene Expression")
  lines(sort(xk),modk$fit[order(xk)],lwd=3,col=2)
  title(paste("Gene", kk))
}
par(mfrow=c(1,1))
############################

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

#     rawp   
#0.05  441    1    1    1    0




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
### plot top gene                                                                        ###
############################################################################################


##################################



i<-610 #this must be gene 610
par(mfrow=c(1,1))
#fit.lm.g<-lm(y~as.numeric(data1[index1[i],])+as.factor(gr))
fit.lm.g<-glm(classlabel~as.numeric(data1[i,]),binomial(link = "logit"))
summary(fit.lm.g)
plot(as.numeric(data1[i,]),classlabel)
points(as.numeric(data1[i,]),fit.lm.g$fit,pch="+",col=2)


xi<-as.numeric(data1[i,])
modi<-glm(classlabel~xi,binomial(link = "logit"))
t.t[i]<-summary(modi)$coefficients[2,3]
p.val[i]<-summary(modi)$coefficients[2,4]


############################################################################################
### plot non significant  gene                                                           ###
############################################################################################


j<-2168
par(mfrow=c(1,1))
#fit.lm.g<-glm(classlabel~as.numeric(data1[index1[i],]),binomial(link = "logit"))
fit.lm.g<-glm(classlabel~as.numeric(data1[j,]),binomial(link = "logit"))
summary(fit.lm.g)
plot(as.numeric(data1[j,]),classlabel)
points(as.numeric(data1[j,]),fit.lm.g$fit,pch="+",col=2)




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
Ux1
Ux2<-as.vector(pc.cr$scores[,2])
Ux2
Ux3<-as.vector(pc.cr$scores[,3])
Ux3
Ux4<-as.vector(pc.cr$scores[,4])
Ux4
Ux5<-as.vector(pc.cr$scores[,5])
Ux6<-as.vector(pc.cr$scores[,6])
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
plot(xx,pred.x,type="l",xlab="gene expression",ylab="P(Cancer)", col=2)
points(Ux1,classlabel)
points(Ux1[classlabel==1],classlabel[classlabel==1])
# top Ux top 5
#    Null deviance: 141.363  on 101  degrees of freedom
#Residual deviance:  66.329  on 100  degrees of freedom
#AIC: 70.329

#top10
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   0.3664     0.4462   0.821    0.412    
#Ux1           3.8518     0.9611   4.008 6.13e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 141.363  on 101  degrees of freedom
#Residual deviance:  34.583  on 100  degrees of freedom
#AIC: 38.583


plot(Ux1,classlabel)
summary(glm(classlabel~Ux1+Ux2,binomial(link = "logit")))
#Ux2 not significant

plot(Ux,classlabel)
summary(glm(classlabel~Ux1+Ux2+Ux3+Ux4,binomial(link = "logit")))
summary(glm(classlabel~Ux1+Ux2+Ux3+Ux4+Ux5+Ux6,binomial(link = "logit")))

############################################################################################
###    Fisher LDA                                                                         ##
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


#####################################
# Leave one out CV                  #
#####################################
beetle1 <- data.RUX[data.RUX$Group == "Cancer",][,2:6]
beetle2 <- data.RUX[data.RUX$Group == "Control",][,2:6]

n1 <- nrow(beetle1)
n2 <- nrow(beetle2)

beetle1.means <- apply(beetle1, 2, mean)
beetle2.means <- apply(beetle2, 2, mean)

w1 <- (n1 - 1) * var(beetle1)
w2 <- (n2 - 1) * var(beetle2)

sp1 <- 1 / (n1 + n2 - 2) * (w1 + w2)
cutoff <- .5 * (beetle1.means - beetle2.means) %*% solve(sp1) %*% (beetle1.means + beetle2.means)
cutoff # -0.0463994

species.prediction <- apply(data.RUX[,2:6], 1, function(y) {
  z <- (beetle1.means - beetle2.means) %*% solve(sp1) %*% y # Calculate the discriminate function for the observation vector y
  ifelse(z > cutoff, 1,0)
})

t(table(data.RUX$classlabel, species.prediction, dnn = c('Actual Group','Predicted Group')))


n <- dim(data.RUX)[1]
1 / n

###########################################

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
xtable(t(table(data.RUX$Group, cv.prediction, dnn = c('Actual Group','Predicted Group'))))
#Actual Group
#Predicted Group Cancer Control
#0     10      45
#1     42       5

install.packages("ldaCMA")
library(ldaCMA)
### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run LDA
ldaresult <- ldaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(ldaresult)
ftable(ldaresult)
plot(ldaresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 10 genes
khanX <- as.matrix(khan[,2:11])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run LDA
ldaresult <- ldaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(ldaresult)
ftable(ldaresult)
plot(ldaresult)

############################################################################################
###    Tree based models                                                                  ##
############################################################################################

library(tree)
tree.fit <- tree(as.factor(classlabel)~Ux+Ux2)
summary(tree.fit)
tree.fit <- tree(as.factor(classlabel)~Ux+Ux2+Ux3+Ux4)
summary(tree.fit)
plot(tree.fit)
text(tree.fit, pretty=0)
cv.treefit <-cv.tree(tree.fit, FUN = prune.misclass)
cv.treefit

plot(cv.treefit$size, cv.treefit$dev, type= "b")
plot(cv.treefit$k, cv.treefit$dev, type= "b")


prune.treef <- prune.misclass(tree.fit, best=2)
plot(prune.treef)
text(prune.treef, pretty=0)

prune.treef

plot(Ux,Ux2)
plot(Ux,Ux2, col=c(rep(1,50), rep(2,52)))
abline(v=0.402409)


##Cross validation
library(caret)


# Define train control for k fold cross validation
train_control <- trainControl(method="cv", number=10)
# Fit Naive Bayes Model
model <- train(Species~., data=iris, trControl=train_control, method="nb")
# Summarise Results
print(model)

