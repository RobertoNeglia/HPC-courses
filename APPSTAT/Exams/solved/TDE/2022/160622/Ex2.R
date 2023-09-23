### Ex1

library(MASS)

dataset = read.table("musiccountry.txt",header=T)

attach(dataset)

group = factor(dataset$release.country, levels=c('Germany', 'US'))   # Species è il nome della colonna del grouping nel dataset
# mettendo i levels sono sicura che l'ordine è quello che penso io

g = 2 # number of groups 


i1 <- which(group=='Germany')       # rows of group1
i2 <- which(group=='US')

n1 <- length(i1)
n2 <- length(i2)

n <- n1+n2

detach(dataset)

dataset_num <- dataset[,1:2]

x11()
plot(dataset_num, main='Data', xlab='V1', ylab='V2', pch=19)
points(dataset_num[i1,], col='red', pch=19)
points(dataset_num[i2,], col='green', pch=19)

legend("topright", legend=levels(group), fill=c('red','blue'))



# gaussianity
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
   # g righe

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)
#round(S3,digits=1)

# to compare visually
x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
#image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))




# Set i prior:
prior.c = c(0.1 , 0.9 )
qda.m <- qda(dataset_num, group, prior = prior.c)
qda.m

x11()
plot(dataset_num, main='Partition', xlab='V1', ylab='V2', pch=20)
points(dataset_num[i1,], col='red', pch=20)   
# se i dati sono in dataframe separati, ad esempio il  gruppo uno è nel dataset g1
# allora semplicemente: dataset_num[i1,] = g1
points(dataset_num[i2,], col='blue', pch=20)
legend("topright", legend=levels(group), fill=c('red','blue'), cex=.7)

#      qui è con qda!
points(qda.m$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)
# the x are the centroids

# i create some fine grid in order to plot the regions that we have obtained 
x  <- seq(min(dataset_num[,1]), max(dataset_num[,1]), length=200)
y  <- seq(min(dataset_num[,2]), max(dataset_num[,2]), length=200)

#                 nome_var1   nome_var2
xy <- expand.grid( price   =x ,  average.length  = y)


z  <- predict(qda.m, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - z[,2] # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - z[,1] # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    


# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
# the function contour with level=0 plots the line 
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)




# means and covariances 
qda.m$means



# 
qdaCV <- qda(dataset_num, group,prior=c(0.1 , 0.9), CV=TRUE)  # specify the argument CV
#this perform the leave one out cross val

table(class.true=group, class.assignedCV=qdaCV$class)

# F= germ
# T= us
4/152*0.9 + 9/36*0.1

# Predire una nuova osservazione
x_new = c( 50, 3.5)  # vettore in Rp

C1_pred <- predict(qda.m, x_new) 
C1_pred


## SVM

library(e1071)
y <- rep(0,n)
y[which(group=='US')] <- 1

dat <- data.frame(x=dataset_num[,c(2,1)], y=as.factor (y)) # the name x and y are mandatory


set.seed(1)
tune.out <- tune(svm,y~.,data=dat ,kernel = 'linear',
                 ranges =list(cost=c(0.001 , 0.01, 0.1, 1,10,100) ) , tunecontrol = tune.control(cross = 10), scale = F)
summary(tune.out)

# Extract the best model from the result of tune
bestmod <- tune.out$best.model
summary(bestmod)
x11()
plot(bestmod , dat, col =c('salmon', 'light blue'), pch=19, asp=1)

testdat = data.frame(x.price = 50, x.average.length = 3.5)
ypred <- predict(bestmod,testdat)
ypred
# CORREGGI SU TEMPLATE







