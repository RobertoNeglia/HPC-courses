#### Ex2


library(MASS)

dataset = read.table("orthopaedics.txt", header = T)

group = factor(dataset$norm_abnorm, levels=c('NO', 'AB'))   
# mettendo i levels sono sicura che l'ordine è quello che penso io

g = 2 # number of groups 


i1 <- which(group=='NO')       # rows of group1
i2 <- which(group=='AB')
 # tante righe quanti gruppi

n1 <- length(i1)
n2 <- length(i2)

n <- n1+n2

dataset_num <- dataset[,1:2]
load('mcshapiro.test.RData')
data_group1 = dataset_num[i1,]
data_group2 = dataset_num[i2,]

P = c( mcshapiro.test(data_group1)$pvalue, 
       mcshapiro.test(data_group2)$pvalue)
P

# multivariate
S  <-  cov(dataset_num)
S1 <-  cov(dataset_num[i1,])   # S nel gruppo i = 1..g
S2 <-  cov(dataset_num[i2,])
x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))


round(S1,digits=1)
round(S2,digits=1)


prior.c = c(0.65, 0.35)
qda.m <- qda(dataset_num, group, prior = prior.c)
qda.m




## Covariances 
x_new = dataset_num # vettore in Rp

C1_pred <- predict(qda.m, x_new) 
C1_pred$class

new_group_NO = which(C1_pred$class == 'NO')
new_group_AB = which(C1_pred$class == 'AB')
S1 <-  cov(dataset_num[new_group_NO,])   # S nel gruppo i = 1..g
S2 <-  cov(dataset_num[new_group_AB,])
S1
S2

x11()
plot(dataset_num, main='Partition', xlab='V1', ylab='V2', pch=20)
points(dataset_num[i1,], col='red', pch=20)   
# se i dati sono in dataframe separati, ad esempio il  gruppo uno è nel dataset g1
# allora semplicemente: dataset_num[i1,] = g1
points(dataset_num[i2,], col='blue', pch=20)
legend("topright", legend=levels(group), fill=c('red','blue'), cex=.7)

points(qda.m$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)
# the x are the centroids

# i create some fine grid in order to plot the regions that we have obtained 
x  <- seq(min(dataset_num[,1]), max(dataset_num[,1]), length=200)
y  <- seq(min(dataset_num[,2]), max(dataset_num[,2]), length=200)
xy <- expand.grid(incidence=x, tilt=y)
# Sepal.Length = NOME (non colonna) variabile 1 
# Sepal.Width = nome variabile 2

z  <- predict(qda.m, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - z[,2] # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - z[,1] # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    


# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
# the function contour with level=0 plots the line 
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)


# b)
qdaCV <- qda(dataset_num, group,prior=c(0.65, 0.35), CV=TRUE)  # specify the argument CV
#this perform the leave one out cross val

table(class.true=group, class.assignedCV=qdaCV$class)

# calcolo AER: leggi i valori dalla tabella
# ad esempio con due classi True False
# P12 = P(assegnare F quando è T)
# P21 = P(assegnare T quando è F)
# pt = prior di T
# pf = prior di F
# AER = P12/totale(T) * pt + P21/totale(F) * pf
# T=N ,F=A

8/80 * 0.65 + 32/70 * 0.35 


# c)
x_new = c(60,0)  # vettore in Rp

C1_pred <- predict(qda.m, x_new) 
C1_pred




# d) 
library(e1071)


y <- rep(0,n)
y[which(group=='AB')] <- 1  # 1 = AB

# dat <- data.frame(x=dataset_num[,c(2,1)], y=as.factor (y)) # the name x and y are mandatory

# dat <- data.frame(x=dataset_num, y=as.factor (y)) # the name x and y are mandatory

dat <- data.frame(x=dataset_num, y= group)

svmfit <- svm(y~., data=dat , kernel ='linear', cost =0.1, scale =FALSE )

summary(svmfit)

x11()
par(mfrow=c(1,2))
plot(svmfit , dat , col =c('salmon', 'light blue'), pch=19, asp=1)


# Prediction
testdat = data.frame(x.incidence = 60 , x.tilt = 0 )
ypred <- predict(svmfit,testdat)
ypred


















