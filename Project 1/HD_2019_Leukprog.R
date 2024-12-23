###################################################################################
# The leukemia dataset from Efron 2010                                            #
# 7128x72 matrix, columns labeled 0 and 1 for AML and ALL                         #
###################################################################################
rm(list=ls())

#########################################################
# packages                                              #
#########################################################



library(limma)
library(statmod)
library(AffyExpress)
library(vsn)
library(splines)
#limmaUsersGuide()
library(Biobase)
library(multtest, verbose = FALSE)


#Data Preparation
load("C:\\Users\\juach\\Documents\\UHasselt\\Analysis of High Dimensional Data\\Project\\Leukemia\\leukdata.RData")
dim(leukdata)
head(leukdata)
summary(leukdata)
str(leukdata)
leukdata[1,][1:47]
leukdata[1,][48:72]
#View(leukdata)


#####################################
#### ANALYSIS WITH t test        ####
#####################################

leukdata1 <-as.data.frame(leukdata)
classlabel<-c(rep(0,47),rep(1,25))
ngenes<-length(leukdata[,1])

r.i<-sd.i<-e.i<-t.t<-p.val<-c(1:ngenes)

## two-sample t test   ##
## two-sample t test   ##
## two-sample t test   ##
i=1
for(i in 1:ngenes)
{	
  cat(i)
  x1<-leukdata1[i,][1:47]
  x2<-leukdata1[i,][48:72]
  t.i<-t.test(x1,x2, alternative="two.sided")
  #t.w<-wilcox.test(x1,x2)  #cannot compute exact p-value with ties
  t.t[i]<-t.i$statistic
  p.val[i]<-t.i$p.value
  #t.tw[i]<-t.w$statistic
  #p.valw[i]<-t.w$p.value
  e.i[i]<-t.i$estimate[1]-t.i$estimate[2]
  r.i[i]<-t.i$estimate[1]/t.i$estimate[2]
  sd.i[i]<-e.i[i]/t.t[i]
}


#######################################
# figures for example in section 7.7 ##
#######################################
head(sd.i)
par(mfrow=c(1,1))
p.val.all<-p.val
e.i.all<-e.i
sum(abs(e.i)>2)
sum(abs(p.val)>0.05)
par(mfrow=c(1,2))
#FC
plot(e.i,xlab="Genes",ylab="Fold Change")
abline(2,0,lty=4)
abline(-2,0,lty=4)
plot(log(r.i),xlab="Genes",ylab="Fold Change")

#Pvalue vs FC
plot(e.i,p.val.all, xlab="Fold Change", ylab="Pvalue")
abline(h=0.05)
plot(e.i,-log(p.val.all), xlab="Fold Change", ylab="-Log(Pvalue)")
abline(h=-log(0.05),lty=4)
hist(p.val.all,col=0,nclass=50,main=" ",xlab="Pvalue")

#Ttest vs FC
plot(t.t,e.i,xlab="T-test Statistic",ylab="Fold Change")
abline(2,0,lty=4)
abline(-2,0,lty=4)
plot(e.i, t.t,ylab="T-test Statistic",xlab="Fold Change")
abline(v=2,lty=4)
abline(v=-2,lty=4)

#FC vs SE
plot(e.i,sd.i,xlab="Fold Change",ylab="Standard Error")
abline(v=2,lty=4)
abline(v=-2,lty=4)
#lines(c(2,2),c(0,1),lty=4)
#lines(c(-2,-2),c(0,1),lty=4)

#Ttest vs SE
plot(t.t,sd.i,xlab="T-test Statistic",ylab="Standard Error")
plot(sd.i, t.t,ylab="T-test Statistic",xlab="Standard Error",xlim=c(0,10))
 
#Ttest vs Pvalue
plot(t.t,p.val.all,xlab="T-test Statistic",ylab="Pvalue")
abline(h=0.05,lty=4)

#Ttest
plot(t.t,xlab="Genes",ylab="T-test Statistic")


#Number of mistakes
B <- 7128
ho <- c(1:7128)
mistakes <- c(1:7128)
for (i in 1:B){
  mistakes[i] <- 1-(1-0.05)^i
}

plot(ho[c(1:200)],mistakes[c(1:200)],xlab="No. of Null Hypothesis",ylab="P(at least one mistake)")



#Up regulated and down regulated
index<-c(1:ngenes)
index1<-index[order(e.i)]
sort(e.i)
down.index <- index1[1:5]
down.fc <-e.i[down.index]
#downregulated
summary(leukdata[index1[1],])
plot(leukdata[index1[1],], type="l", ylim = c(-3,3), xlab="Gene", ylab="Gene Expression")
for (i in 2:5){
  lines(leukdata[index1[i],],col=i) 
}
abline(v=47.5)

#upregulated
up.index <- index1[7123:7128]
up.fc <-e.i[up.index]
plot(leukdata[index1[7123],], type="l", ylim = c(-3,3), xlab="Gene", ylab="Gene Expression")
for (i in 7124:7128){
  lines(leukdata[index1[i],],col=i) 
}
abline(v=47.5)

ii<-index1[1]
pi<-p.val[ii]
head(sort(e.i))
tail(sort(e.i))


#######################################
# multiplicity                   7.7 ##
#######################################
par(mfrow=c(1,1))
library(Biobase)
library(multtest, verbose = FALSE)
rawp<-p.val

holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
#by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2])
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH"),leg=c(5800,0.5),lty=1,col=1:4,lwd=2,cex=0.8)
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2]),seq(0,1,0.05))$r
abline(h=0.05, lty=4)
rejectable<-mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2]),seq(0,1,0.05))$r
sum(allp[,4]<0.05)
#rawp                    
#0.05 2054  246  247 1249  617
library(xtable)
xtable(rejectable,digits=c(2,0,0,0,0,0))


#Upregulated or downregulated
small.index <- index1[1:10] #1882 6854 6041 3252 4847 2121  758 1745 4377 2020
small.t <-t.t[small.index]
sum(small.t<0)
sum(small.t>=0)

small.e <-e.i[small.index]
sum(small.e<0)
sum(small.e>=0)


#Usmallest pvalues
par(mfrow=c(1,1))
index<-c(1:ngenes)
index1<-index[order(p.val)]
small.index <- index1[1:10] #1882 6854 6041 3252 4847 2121  758 1745 4377 2020
small.p <-p.val[small.index]

summary(leukdata[index1[1],])
plot(leukdata[index1[1],], type="l", ylim = c(-3,3), xlab="Gene", ylab="Gene Expression")
for (i in 2:10){
  lines(leukdata[index1[i],],col=i) 
}

Index<-as.factor(c(rep(1,47), rep(2,25)))
levels(Index)<-c("Control", "Cancer")
par(mfrow=c(2,5))
for(i in 1:10){
a1<-as.vector(leukdata[index1[i],1:47])
a2<-as.vector(leukdata[index1[i],48:72])
a12<-as.data.frame(cbind((Index), c(a1, a2)))
a12$Index[a12$V1==1,]=="Control"
a12$Index[a12$V1==2,]=="Cancer"
summary(a12)
#title(paste("Gene",small.index[i]))
boxplot(V2~Index, data=a12, main=paste("Gene",small.index[i]),  col=c("skyblue3", "palevioletred"))
}

#combine all
allvalues <-as.data.frame(cbind(e.i,t.t,sd.i,p.val,allp))
str(allvalues)
colnames(allvalues)<-c("FC","TTest","SD","Pvalue","Rawpval","Bonferroni","Holm", "BH","BY")

index<-c(1:ngenes)
index1<-index[order(allvalues$BH)]
small.index <- index1[1:10] #1882 6854 6041 3252 4847 2121  758 1745 2020 4377


mask <- with(allvalues, abs(FC) < 2 & Pvalue < .01)
sum(mask==TRUE)
spike <- with(allvalues, abs(FC) > 2 & Pvalue < .01)
sum(spike==TRUE)
cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))

with(allvalues, plot(FC, -log10(Pvalue), cex=.8, pch=16,
               xlim=c(-1,1), ylim=c(0,10),
               xlab="difference in means",
               col=cols))
abline(h=2,v=c(-.2,.2), lty=2)


dim(leukdata)
names(as.data.frame(leukdata))
summary(leukdata)
par(mfrow=c(1,1))
#########################################################################################

## one-way ANOVA ##
## one-way ANOVA ##
## one-way ANOVA ##


data1<-leukdata
ngenes<-length(data1[,1])
p.val<-c(1:ngenes)
for(i in 1:ngenes)
{	
  cat(i)
  fit.aov<-aov(data1[i,]~as.factor(classlabel))
  p.val[i]<-anova(fit.aov)[1,5]
}

##########################################################
## MULTIPLICITY                                          #
##########################################################

library(Biobase)
library(multtest, verbose = FALSE)
rawp<-p.val
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),leg=c(5500,0.35),lty=1,col=1:5,lwd=2,cex=0.8)
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r
abline(0.05,0)

#rawp
#0.05 2028  264  269 1233  644

##########################################################
## LIMMA                                                 #
##########################################################
par(mfrow=c(1,1))
data1<-leukdata1
ngenes<-length(data1[,1])
gene.exp<-data1
dim(gene.exp)
head(gene.exp)
nn<-length(classlabel)
nn
library(limma)

design <- cbind(rep(1,72),abs(classlabel-1))
design

dim(gene.exp)
dim(design)
fit <- lmFit(gene.exp,design)
summary(fit)

fit1 <- eBayes(fit)
head(fit1$t)
head(fit1$p.value)

topTable <- topTable(fit1,number=nrow(gene.exp),coef=2)
?head
head(topTable, n=10) #3252  2288 1834  6854 1882  760 804 2354 1144 1685 
dim(topTable)
dim(topTable[topTable[,4]<=0.05,]) #2037
dim(topTable[topTable[,5]<=0.05,]) #1248
limma.bh <-as.vector(rownames((topTable[topTable[,5]<=0.05,])))
pvalbh<-index1[1:1249]
ttest.bh<-pvalbh

tst<-c(unique(limma.bh),unique(ttest.bh))
tst<-tst[duplicated(tst)]
tst[duplicated(tst)]
length(tst)

volcanoplot(fit1,coef=2, highlight=5)
help(topTable)

length(topTable[,4])
plot(t.t,-log(p.val), xlab="Test Statistic", ylab="-Log(PValue)")
points(topTable[,3],-log(topTable[,4]),col=2)
#how to add legend here


#check upregulated and downregulated
lim<-topTable[topTable[,5]<=0.05,]
lim<-topTable[1:1248,]
sum(lim[,1]<0)
sum(lim[,1]>=0)



###########


#Usmallest pvalues
par(mfrow=c(1,1))
index<-c(1:ngenes)
index1<-c(3252,  2288, 1834,  6854, 1882,  760, 804, 2354, 1144, 1685 )
#small.index <- index1[1:10] 
#small.p <-p.val[small.index]

Index<-as.factor(c(rep(1,47), rep(2,25)))
levels(Index)<-c("Control", "Cancer")
par(mfrow=c(2,5))
for(i in 1:10){
  a1<-as.vector(leukdata[index1[i],1:47])
  a2<-as.vector(leukdata[index1[i],48:72])
  a12<-as.data.frame(cbind((Index), c(a1, a2)))
  a12$Index[a12$V1==1,]=="Control"
  a12$Index[a12$V1==2,]=="Cancer"
  summary(a12)
  #title(paste("Gene",small.index[i]))
  boxplot(V2~Index, data=a12,  col=c("skyblue3", "palevioletred"), main=paste("Gene",index1[i]))
}
#######################################
# figures                            ##
#######################################


par(mfrow=c(1,1))
#plot(t.t,fit1$t[,2],xlab="t-test statistics",ylab="moderated t-test statistics") 
#abline(0,1)
#title("a")
plot(t.t,-log(p.val), xlab="T-test Statistic", ylab="-Log(PValue)")
points(topTable[,3],-log(topTable[,4]),col=8,pch="+")
#title("b")


par(mfrow=c(1,3))
plot(log(p.val),log(fit1$p.value[,2]),pch=".",xlab="log(raw p value): t-test",ylab="log(raw p value): limma")
abline(log(0.05),0)
lines(c(log(0.05),log(0.05)),c(-50,1))

plot(rank(p.val),rank(fit1$p.value[,2]),xlab="rank of raw p value: t-test",ylab="rank of raw p value: limma")


bh.1<-mt.rawp2adjp(p.val, proc=c("BH"))
summary(bh.1)
head(bh.1)
bh.1$adjp[,2]

p.val2<-as.vector(fit1$p.value[,2])
bh.2<-mt.rawp2adjp(p.val2, proc=c("BH"))
plot(log(bh.1$adjp[,2]),log(bh.2$adjp[,2]),xlab="log(adjusted p value): t-test",ylab="log(adjested p value): limma")
abline(log(0.05),0)
lines(c(log(0.05),log(0.05)),c(-50,1))



