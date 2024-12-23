##############################################################################################
# The prostat dataset from Efron 2010 and Hastie and Efron 2016                             ##
# The first 50 columns are the genetic activity measurements for the 50 control subjects,   ##
# while the last 52 columns represent the prostate cancer subjects.                         ##
##############################################################################################

library(BVS)
library(Biobase)
library(multtest, verbose = FALSE)
library(CMA)
library(plsgenomics)
library(randomForest)
library(MASS)
library(gbm)
library(class)
library(e1071)
par(mfrow=c(1,1))

rm(list=ls())
prostmat <- read.csv("C:\\Users\\juach\\Documents\\UHasselt\\Analysis of High Dimensional Data\\Project\\prostmat.csv")
dim(prostmat)

classlabel = c(rep("Control",50), rep("Cancer",52))
data.a <- as.data.frame(t(as.matrix(prostmat)))
data1 <- as.data.frame(cbind(classlabel, data.a))
str(data1)                        
head(data1[c(1:5),c(1:5)])

#Learning Sets
prostateY <- data1[, 1]
prostateX <- as.matrix(data1[, -1])

# Initials
numfeatures<-c(2,3,5,10,12,15,20,25,28,30,31,32,33,35,40,45,50,60)
#numfeatures<-seq(3,60,2)
#numfeatures<-c(2,3,10,20,30,40,50)
#numfeatures=c(5,10)
n.f<-length((numfeatures))
n.fold = 3
n.it = 1000
n.all = n.fold*n.it
mce<- matrix(data=NA, nrow=n.f, ncol=8)
spec<- matrix(data=NA, nrow=n.f, ncol=8)
sens<- matrix(data=NA, nrow=n.f, ncol=8)
auc<- matrix(data=NA, nrow=n.f, ncol=8)
eval.m<- matrix(data=NA, nrow=n.f, ncol=8)
se.m<- matrix(data=NA, nrow=n.f, ncol=8)
rf.auc<-matrix(data=NA, nrow=n.all, ncol=n.f)
w.name<-NA
par(mfrow=c(1,1))
#fiveCV5iter <- GenerateLearningsets(y = prostateY, method = "CV", fold = n.fold, niter = n.it, strat = TRUE)


# Looping
set.seed(123)
for (i in 1:n.f){
  num.top <- numfeatures[i]
  fiveCV5iter <- GenerateLearningsets(y = prostateY, method = "CV", fold = n.fold, niter = n.it, strat = TRUE)
  varsel_fiveCV <- GeneSelection(X = prostateX, y = prostateY, learningsets = fiveCV5iter, method = "limma")
  show(varsel_fiveCV)
  toplist(varsel_fiveCV, iter = 1)
  seliter <- numeric()
  for (j in 1:n.all) 
  {seliter <- c(seliter, toplist(varsel_fiveCV, iter = j,top = num.top, show = FALSE)$index)}
  sort(table(seliter), dec = TRUE)
  genesel_lim <- varsel_fiveCV
  
  
  #barplot
  len.sel=dim(table(seliter))-num.top
  Colo<-c(rep(2,num.top),rep(NA,len.sel))
  W.Color<-as.data.frame(cbind(as.data.frame(sort(table(seliter), dec = TRUE)),Colo))
  W.Color$seliter<-as.numeric(as.character(W.Color$seliter))
  cols = W.Color[order(W.Color$seliter),]
  barplot(table(seliter) ,col=cols$Colo, main=paste("Top",num.top,"Genes"), ylim=c(0,n.all), cex.axis=0.75)#xaxt="none"
  #axis(1, seq(1,6000),las=2, cex.axis=0.8, font=2)
  # try to put in a table top genes per run
  w.name<-cbind(w.name,W.Color[1:50,1:2])
  
  ### Classification Models ####
  # DLDA
  class_dlda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = dldaCMA, nbgene=num.top)
  a<-evaluation(class_dlda, measure = "misclassification", scheme="iter")
  eval.m[i,1]<-mean(a@score)
  se.m[i,1] <-sd(a@score)/sqrt(3000)
  
  # LDA
  class_lda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = ldaCMA, genesel = genesel_lim, nbgene = num.top)
  a<-evaluation(class_lda, measure = "misclassification", scheme="iter")
  eval.m[i,2]<-mean(a@score)
  se.m[i,2] <-sd(a@score)/sqrt(3000)
  
  #class_fda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = fdaCMA, genesel = genesel_lim, nbgene = num.top, comp = 2)
  #class_qda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = qdaCMA, genesel = genesel_lim, nbgene = 1)
  
  #Shrunken centroids discriminant analysis does not do variable selection, but hyperparameter tuning for the shinkage intensity.
  class_scda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = scdaCMA) #tuninglist = list(grids = list()))
  a<-evaluation(class_scda, measure = "misclassification", scheme="iter")
  eval.m[i,3]<-mean(a@score)
  se.m[i,3] <-sd(a@score)/sqrt(3000)
  
  
  # partial least squares (with two components)
  class_plsda <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = pls_ldaCMA) #tuninglist = list(grids = list()))
  a<-evaluation(class_plsda, measure = "misclassification", scheme="iter")
  eval.m[i,4]<-mean(a@score)
  se.m[i,4] <-sd(a@score)/sqrt(3000)
  
  
  #KNN
  class_knn <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, classifier = pknnCMA, genesel = genesel_lim, nbgene = num.top) #tuninglist =  list(grids = list()),
  a<-evaluation(class_knn, measure = "misclassification", scheme="iter")
  eval.m[i,5]<-mean(a@score)
  se.m[i,5] <-sd(a@score)/sqrt(3000)
  
  
  #svm
  class_svm <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter, genesel = genesel_lim, nbgene = num.top, classifier = svmCMA,  probability=TRUE) #tuninglist =  list(grids = list()), 
  a<-evaluation(class_svm, measure = "misclassification", scheme="iter")
  eval.m[i,6]<-mean(a@score)
  se.m[i,6] <-sd(a@score)/sqrt(3000)
  
  
  # LASSO
  class_lasso <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter,  genesel = genesel_lim, classifier = LassoCMA,nbgene=num.top) #tuninglist = list(grids = list())
  a<-evaluation(class_lasso, measure = "misclassification", scheme="iter")
  eval.m[i,7]<-mean(a@score)
  se.m[i,7] <-sd(a@score)/sqrt(3000)
  
  
  # Random Forest
  class_rf <- classification(X = prostateX, y = prostateY, learningsets = fiveCV5iter,  genesel = genesel_lim, classifier = rfCMA, nbgene=num.top) #tuninglist = list(grids = list()),
  a<-evaluation(class_rf, measure = "misclassification", scheme="iter")
  eval.m[i,8]<-mean(a@score)
  se.m[i,8] <-sd(a@score)/sqrt(3000)
  b<-evaluation(class_rf, measure = "auc", scheme="iter")
  rf.auc[,i]<- b@score
  
  
  #comparison
  #dalike <- list(class_dlda, class_lda, class_fda, class_qda, class_scda, class_plsda, class_knn, class_svm, class_lasso, class_rf)
  dalike <- list(class_dlda, class_lda, class_scda, class_plsda, class_knn, class_svm, class_lasso, class_rf)
  
  par(mfrow = c(2,2))
  golub.comp1 <- compare(dalike, plot = TRUE, measure = c("misclassification","auc", "sensitivity", "specificity" ))
  print(golub.comp1)
  
  # MCE, AUC, SENS, SPEC
  mce[i,]<-golub.comp1[[1]]
  auc[i,]<-golub.comp1[[2]]
  sens[i,]<-golub.comp1[[3]]
  spec[i,]<-golub.comp1[[4]]
  
  
  # ROC Plots
  schemes <- c("DLDA","LDA","SCDA","PLSDA","KNN","SVM","LASSO","RF")
  par(mfrow=c(2,2))
  result <- lapply(dalike, join)
  #for(i in seq(along = result)) plot(result[[i]], main = schemes[i])
  for(i in seq(along = result)) roc(result[[i]],main=schemes[i])
  
  par(mfrow=c(1,1))
}



numfeatures
rf.auc.dt<-as.data.frame(rf.auc[,-c(5,9,11,12,13,18)])
names(rf.auc.dt)
colnames(rf.auc.dt)<-numfeatures[-c(5,9,11,12,13,18)]
boxplot(rf.auc.dt, ylab="AUC", xlab="Number of Features" )

summary(rf.auc.dt)

library(ggplot2)
# Basic box plot
p <- ggplot(rf.auc.dt, aes(x=numfeatures)) +   geom_boxplot()
p


a<-evaluation(class_rf, measure = "misclassification", scheme="iter")
str(a)
eval.rf <-mean(a@score)
se.rf <-sd(a@score)/sqrt(1000)

r2<-as.data.frame(cbind(numfeatures,as.data.frame(r)))
colnames(r2)<-c("Top K","DLDA","LDA","SCDA","PLSDA","KNN","SVM","LASSO","RF")

schemes <- c("DLDA","LDA","SCDA","PLSDA","KNN","SVM","LASSO","RF")
par(mfrow=c(2,2))
result <- lapply(dalike, join)
#for(i in seq(along = result)) plot(result[[i]], main = schemes[i])
for(i in seq(along = result)) roc(result[[i]],main=schemes[i])

par(mfrow=c(1,1))
#MCE
mce.dt<-cbind(numfeatures,as.data.frame(mce))
str(mce.dt)
mce.dt<-mce.dt[-c(18),]
plot(mce.dt[,2]~mce.dt[,1], ylim=c(0,1),type="l", xlab="Top K Genes", ylab="MCE", col=1,lwd=2)
for (i in 3:dim(mce.dt)[2]){
  lines(mce.dt[,i]~mce.dt[,1],col=i-1,,lwd=2)
}
legend("topleft",legend = schemes, col=(1:8), cex = 0.7 , lty=1,lwd=2)

par(mfrow=c(1,2))
#Sensitivity
sens.dt<-cbind(numfeatures,as.data.frame(sens))
sens.dt<-sens.dt[-c(18),]
plot(sens.dt[,2]~sens.dt[,1], ylim=c(0,1),type="l", xlab="Top K Genes", ylab="Sensitivity",lwd=2)
for (i in 3:dim(sens.dt)[2]){
  lines(sens.dt[,i]~sens.dt[,1],col=i-1,lwd=2)
}
legend("bottomleft",legend = schemes, col=(1:8), cex = 0.7 , lty=1,lwd=2)


#Specificity
spec.dt<-cbind(numfeatures,as.data.frame(spec))
spec.dt<-spec.dt[-c(18),]
plot(spec.dt[,2]~spec.dt[,1], ylim=c(0,1),type="l", xlab="Top K Genes", ylab="Specificity",lwd=2)
for (i in 3:dim(spec.dt)[2]){
  lines(spec.dt[,i]~spec.dt[,1],col=i-1,lwd=2)
}
legend("bottomleft",legend = schemes, col=(1:8), cex = 0.7 , lty=1,lwd=2)


#AUC
auc.dt<-cbind(numfeatures,as.data.frame(auc))
auc.dt<-auc.dt[-c(18),]
plot(auc.dt[,2]~auc.dt[,1], ylim=c(0,1),type="l", xlab="Number of Features", ylab="AUC",lwd=2)
for (i in 3:dim(auc.dt)[2]){
  lines(auc.dt[,i]~auc.dt[,1],col=i-1,lwd=2)
}
legend("bottomleft",legend = schemes, col=(1:8), cex = 0.7 , lty=1,lwd=2)


dalike <- list(class_dlda, class_lda, class_scda, class_plsda, class_knn, class_svm, class_lasso, class_rf)
dalike2 <- c("class_dlda", "class_lda")#, class_scda, class_plsda, class_knn, class_svm, class_lasso, class_rf)
dalike2[1]
evaluation(dalike2[1])
a<-evaluation(class_dlda, measure = "misclassification", scheme="iter")
eval.m[1,1]<-mean(a@score)
se.m[1,1] <-sd(a@score)/sqrt(3000)


###############################
#Top 30 Genes using Random FOrest
top30<-c(610, 1720, 3940, 914, 364, 332, 3647, 579, 4331, 1089, 4546, 1068, 1077, 3991, 735, 3375, 4088, 702, 3665, 4073, 4316, 1113, 921, 2, 739, 4518, 2945, 694, 721, 3282)
length(top30)
str(data1)
library(ggplot2)
# Basic box plot

p <- ggplot(data1, aes(x=classlabel, y=data1[,top30[1]])) +   geom_boxplot()
p

boxplot(data1$top30[1])
head(data1[,top30[1]])
boxplot(data1$top30[1], data1$classlabel)

index1<-c(610, 1720, 3940, 914, 364, 332, 3647, 579, 4331, 1089, 4546, 1068, 1077, 3991, 735, 3375, 4088, 702, 3665, 4073, 4316, 1113, 921, 2, 739, 4518, 2945, 694, 721, 3282)
Index<-as.factor(c(rep(1,50), rep(2,52)))
levels(Index)<-c("Control", "Cancer")
par(mfrow=c(2,5))
i=1
for(i in 1:10){
  a1<-as.vector(prostmat[index1[i],1:50])
  a2<-as.vector(prostmat[index1[i],51:102])
  a12<-as.data.frame(cbind((Index), c(a1, a2)))
  a12$Index[a12$V1==1,]=="Control"
  a12$Index[a12$V1==2,]=="Cancer"
  summary(a12)
  str(a12)
  #title(paste("Gene",small.index[i]))
  boxplot(V2~Index, data=a12, main=paste("Gene",small.index[i]))
}


