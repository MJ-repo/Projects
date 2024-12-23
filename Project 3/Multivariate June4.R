###################################################
###############Data Preparation####################
###################################################
rm(list=ls())
load("CanadianWeather.rda")
da<-CanadianWeather[[1]]


#Average Daily Rainfall (mm/day) for 365 days, 35 Canadian Cities
da<-da[,,"Precipitation.mm"] # precipitation data
head(da)
colnames(da)

str(da)

#City, Region, Province, Latitude, Longitude (35x5 matrix)
MetaData<-data.frame(city=colnames(da), region=CanadianWeather$region, province=CanadianWeather$province, coord=CanadianWeather$coordinates)
head(MetaData)
str(MetaData)

library(ldr)
# Using the polynomial basis functions up to degree 18
days<-1:365
days<-(days-min(days))/(diff(range(days))) # rescaling to [0,1]
phi<-poly(days,degree=18)#polynomial
phi<-bf(days,case="fourier",degree=18)#polynomial
dim(phi)

for (i in 1:dim(da)[2]){
  # estimation of the theta parameters for Vancouver
  m.i<-lm(da[,i]~phi)
  summary(m.i)
  
  # plot of fitted function
  #plot(1:365,da[,i],main=colnames(da)[i],"(m=10)", xlab="day", ylab="precipitation (mm/day)")
  #lines(1:365,m.i$fitted.values,type="l", col=2)
  
  if (i==1){
    dfrow <- data.frame(City=colnames(da)[i],t(m.i$coefficients))  
    dfcol <- data.frame(m.i$coefficients)  
    #    dfrsqd <- data.frame(City=colnames(da)[i],R2=summary(m.i)$r.squared*100)
  }else {
    dfrow<-rbind(dfrow,data.frame(City=colnames(da)[i],t(m.i$coefficients))) 
    dfcol<-cbind(dfcol,data.frame(m.i$coefficients)) 
    #    dfrsqd <-rbind(dfrsqd,data.frame(City=colnames(da)[i],summary(m.i)$r.squared*100))
  }
}
colnames(dfcol)=colnames(da)

#smoothing
str(phi2)
theta<-dfrow[,-1]#Using the data bigtheta (dfrow) (35x19)
ybar<-as.data.frame(t(apply(theta,2,mean))) #mean of original theta 1x(18+1) matrix
phi2<-cbind(replicate(365,1), phi) #Add a vector of ones to the first col of the phi matrix 
Ybar<-as.matrix(ybar)%*%t(phi2)  #Create the vector of average values to make the back transformation

t<-as.data.frame(t(apply(theta,2,mean))) #mean of original theta 1x(18+1) matrix
x_1<-matrix(1,nrow=365,ncol=1)
x_orig<-as.data.frame(cbind(x_1,as.matrix(phi))) #Xphi regressors with intercept 365x(18+1) matrix
str(x_orig)
tX<-as.matrix(t)%*%t(as.matrix(x_orig)) #newmean data to replace ybar for smooth curve 1x365
str(tX)

n<-nrow(theta) 
H<-diag(n)-1/n*matrix(1,ncol=n,nrow=n) 
theta[,]<-H%*%as.matrix(theta)  #column centering
theta.svd<-svd(theta) 
k<-2 
Uk<-theta.svd$u[,1:k] 
Dk<-diag(theta.svd$d[1:k]) 
Vk<-theta.svd$v[,1:k] 
thetak<-Uk%*%Dk%*%t(Vk)

#MDS Plot
Zk<-Uk%*%Dk 
rownames(Zk)<-dfrow[,1]
abbr=abbreviate(rownames(Zk))

plot(Zk,type="n",xlab="Z1",ylab="Z2", xlim=c(-50,25)) 
text(Zk,abbr,cex=0.7,col=as.integer(MetaData$region))
legend("topright", legend=unique(MetaData$region),pch=1,col=unique(MetaData$region),
       cex=0.7)


###########################################################
############Backtransformation and Biplots ################
###########################################################

str(Ybar)
str(V1.max+Ybar)
View(V1.max)
V1.max<- (max(Zk[,1])) %*% t(Vk[,1]) %*% t(phi2) #Create the vector to be plot when Zk is max 
V1.min<- (min(Zk[,1])) %*% t(Vk[,1]) %*% t(phi2) #Create the vector to be plot when Zk is min 
plot(days, Ybar,type='l', col=1, ylim=c(0,11)) 
lines(days, (V1.max)+Ybar, col=2) #plot the max values adding the average 
lines(days, (V1.min)+Ybar, col=3) #plot the min values adding the average


V2.max<- (max(Zk[,2])) %*% t(Vk[,2]) %*% t(phi2)#Create the vector to be plot when Zk is max 
V2.min<- (min(Zk[,2])) %*% t(Vk[,2]) %*% t(phi2)#Create the vector to be plot when Zk is min 
plot(days, Ybar,type='l', col=1, ylim=c(-1,8)) 
lines(days, (V2.max)+Ybar, col=2) 
lines(days, (V2.min)+Ybar, col=3)
str(v1.max)


#Y average per day 35x365 matrix
data2<-as.data.frame(t(da))
str(data2)
ybar<-as.data.frame(t(apply(data2,2,mean)))
str(ybar)
for (i in 2:dim(data2)[1]){
  ybar<-rbind(ybar,as.data.frame(t(apply(data2,2,mean))))
}
dim(ybar)
View(ybar)
ybarm<-as.matrix(ybar)
str(ybarm)

str(phi)    #365x18 matrix
str(theta)  #35x19 matrix
str(tX)
#using the new tX matrix
tXnew = tX
for (i in 2:dim(data2)[1]){
  tXnew<-rbind(tXnew,tX)
}
str(tXnew)
View(tXnew)
tXnewm<-as.matrix(tXnew)
str(tXnewm)
#i need to get the matrix X (m+1)xp!!!!
#Y=THETA*transpose(X)
#Apply biplot in Y!!!
x_1<-matrix(1,nrow=365,ncol=1)
x_orig<-as.data.frame(cbind(x_1,as.matrix(phi)))
str(x_orig)

#backtransform!!!! yhat = thetacentered*transpose(X) + tXnew (replacing Ybar(old)) 
y_new<-as.matrix(theta)%*%t(as.matrix(x_orig)) + tXnew #+ybarm
str(y_new)
head(y_new)
#same as:::
#backtransform with truncated SVD k=2: yhat = Zk*transpose(Vk)*transpose(X) + tXnew (replacing Ybar(old)) 
y_new2<-Zk%*%t(Vk)%*%t(as.matrix(x_orig))+ tXnew #+ybarm
str(y_new)
head(y_new)



########################
########Z1 Plot#########

View(Zk)
str(Zk) #35x2
str(Vk) #19x2

y_z1<-Zk[,1]%*%t(Vk[,1])%*%t(as.matrix(x_orig))+tXnewm#original
y_z1<-Zk[,1]%*%t(Vk[,1])%*%t(as.matrix(x_orig))+tXnewm#new

str(x_orig)#365x19
str(tXnewm)#35x365


    #new<-t(as.matrix(x_orig))%*%t(tXnewm)
  #y_z1<-Zk[,1]%*%t(Vk[,1])%*%new
View(y_z1)
#str(y_z1)
#head(y_z1)
#precip(y_z1)

yhat<-as.data.frame(y_z1)
str(yhat)
rownames(yhat)<-colnames(da)
#head(yhat)  

yhatm<-as.data.frame(t(yhat))
#str(yhatm)

plot(yhatm[,1],type="n",ylim=c(-2,15),ylab="Precipitation(in mm)", xlab="Days", main="Precipitation Level (based on Z1 Scores)")
for (i in c(23,29)) {
  lines(yhatm[,i],col=i,lwd=5)
}
yave<-as.data.frame(t(tX))
#str(yave)
lines(yave,col=1,lwd=5)
legend("topright", lty=1, col=c(23,29,1), legend=c("Highest Z1Score","Lowest Z1Score","Average Precipitation"))

########################
########Z2 Plot#########
y_z2<-Zk[,2]%*%t(Vk[,2])%*%t(as.matrix(x_orig))+tXnewm
#str(y_z2)
#head(y_z2)
#precip(y_z2)
yhat<-as.data.frame(y_z2)
#str(yhat)
rownames(yhat)<-colnames(da)
#head(yhat)  
yhatm<-as.data.frame(t(yhat))
#str(yhatm)


plot(yhatm[,1],type="n",ylim=c(-2,15),ylab="Precipitation(in mm)", xlab="Days", main="Precipitation Level (based on Z2 Scores)")
for (i in c(27,29)) {
  lines(yhatm[,i],col=i,lwd=5)
  #text(rownames(Zk)[i],col=i) not working!!!!
}
yave<-as.data.frame(t(tX))
#str(yave)
lines(yave,col=1,lwd=5)
legend("topright", lty=1, col=c(29,27,1), legend=c("Highest Z2Score","Lowest Z2Score","Average Precipitation"))


ggplot(data)+
  geom_smooth(aes(x=Days, y=Precipitation, group = City, color = City), 
              method = 'lm', formula = y ~ poly(x,12), size = 1) +
  geom_vline(xintercept = c(91.25,91.25*2, 91.25*3, 91.25*4), colour = "grey85", alpha = 0.5)+
  annotate("text", x = c(91.25-91.25*1/2,91.25*2-91.25*1/2, 91.25*3-91.25*1/2, 91.25*4-91.25*1/2)
           , y = 11, label = c('Quater 1','Quater 2', 'Quater 3', 'Quater 4'), colour = 'grey70', size = 4)+
  theme(
    text = element_text(size = 10),
    panel.background = element_blank(),
    legend.position = 'top'
  )

mos <- seq( as.Date("2019-01-01"), as.Date("2019-12-31"), by="+1 day",format="%b")
fmos<- format(mos, "%b")
View(fmos)


######################################

data2<-as.data.frame(t(da))
str(data2)
ybar<-as.data.frame(t(apply(data2,2,mean)))
str(ybar)
for (i in 2:dim(data2)[1]){
  ybar<-rbind(ybar,as.data.frame(t(apply(data2,2,mean))))
}
dim(ybar)
View(ybar)
ybarm<-as.matrix(ybar)
str(ybarm)

str(phi)    #365x18 matrix
str(theta)  #35x19 matrix

#using the new tX matrix
tXnew = tX
for (i in 2:dim(data2)[1]){
  tXnew<-rbind(tXnew,tX)
}
str(tXnew)
View(tXnew)
tXnewm<-as.matrix(tXnew)
str(tXnewm)
#i need to get the matrix X (m+1)xp!!!!
#Y=THETA*transpose(X)
#Apply biplot in Y!!!
x_1<-matrix(1,nrow=365,ncol=1)
x_orig<-as.data.frame(cbind(x_1,as.matrix(phi)))
str(x_orig)

#backtransform!!!! yhat = thetacentered*transpose(X) + tXnew (replacing Ybar(old)) 
y_new<-as.matrix(theta)%*%t(as.matrix(x_orig)) + tXnew #+ybarm
str(y_new)
head(y_new)
#same as:::
#backtransform with truncated SVD k=2: yhat = Zk*transpose(Vk)*transpose(X) + tXnew (replacing Ybar(old)) 
y_new2<-Zk%*%t(Vk)%*%t(as.matrix(x_orig))+ tXnew #+ybarm
str(y_new)
head(y_new)



y_z1<-Zk[,1]%*%t(Vk[,1])%*%t(as.matrix(x_orig))+tXnewm
#str(y_z1)
#head(y_z1)
#precip(y_z1)

yhat<-as.data.frame(y_z1)
#str(yhat)
rownames(yhat)<-colnames(da)
#head(yhat)  

yhatm<-as.data.frame(t(yhat))
#str(yhatm)
mos <- seq( as.Date("2019-01-01"), as.Date("2019-12-31"), by="+1 day",format="%b")
fmos<- format(mos, "%b")

plot(yhatm[,1],type="n",ylim=c(-2,15),ylab="Precipitation(in mm)", xlab="Days", main="Precipitation Level (based on Z1 Scores)")
for (i in c(23,29)) {
  lines(yhatm[,i],col=i,lwd=5)
}
yave<-as.data.frame(t(tX))
#str(yave)
lines(yave,col=1,lwd=5)
legend("topright", lty=1, col=c(23,29,1), legend=c("Highest Z1Score","Lowest Z1Score","Average Precipitation"))


#trial
View(yhatz1)
str(y_z1)
str(fmos)
y_z1<-as.data.frame((max(Zk[,1])%*%t(Vk[,1])%*%t(as.matrix(x_orig))+tX))
yhatz1<-cbind(fmos,yhatm)

plot(yhatz1[,1],yhatz1[,2],type="n",ylim=c(-2,15),ylab="Precipitation(in mm)", xlab="Days", main="Precipitation Level (based on Z1 Scores)")
for (i in c(23,29)) {
  lines(yhatm[,i],col=i,lwd=5)
}
yave<-as.data.frame(t(tX))
#str(yave)
lines(yave,col=1,lwd=5)
legend("topright", lty=1, col=c(23,29,1), legend=c("Highest Z1Score","Lowest Z1Score","Average Precipitation"))

t<-as.data.frame(t(apply(theta,2,mean))) #mean of original theta 1x(18+1) matrix
x_orig<-as.data.frame(cbind((matrix(1,nrow=365,ncol=1)),as.matrix(phi))) #Xphi with intercept 365x(18+1) matrix
tX<-as.matrix(t)%*%t(as.matrix(x_orig)) #newmean data to replace ybar for smooth curve
Z1max<-(max(Zk[,1])) %*% t(Vk[,1]) %*% t(x_orig)+tX #Create the vector to be plot when Zk is max 
Z1min<-(min(Zk[,1])) %*% t(Vk[,1]) %*% t(x_orig)+tX #Create the vector to be plot when Zk is min 
days<-1:365
Z1max<- (max(Zk[,1])) %*% t(Vk[,1]) %*% t(x_orig)+tX #Create the vector to be plot when Zk is max 
Z1min<- (min(Zk[,1])) %*% t(Vk[,1]) %*% t(x_orig)+tX #Create the vector to be plot when Zk is min 
plot(days, tX,type='l', col=1, ylim=c(-2,15),ylab="Precipitation(in mm)", xlab="Days", main="Precipitation Level (based on Z1 Scores)") 
lines(days, (Z1max), col=6) #plot the max values adding the average 
lines(days, (Z1min), col=4) #plot the min values adding the average
legend("topright", lty=1, col=c(6,4,1), legend=c("Highest Z1Score","Lowest Z1Score","Average Precipitation"))

