## Ex2

dataset <- read.table('activity.txt', header = T)  # se serve aggiungi il nome delle colonne col.names=c('',..)
dataset

load("mcshapiro.test.RData")

group = factor(dataset$activity, levels=c('walking', 'sitting', 'laying'))   # Species è il nome della colonna del grouping nel dataset
# mettendo i levels sono sicura che l'ordine è quello che penso io

g = 3 # number of groups 

# a) hyp
attach(dataset)

i1 <- which(activity=='walking')   # rows of group A: group è il nome della colonna del grouping nel dataset
i2 <- which(activity=='sitting') 
i3 <- which(activity=='laying') 

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n <- n1+n2+n3

detach(dataset)

dataset_num <- dataset[,1:2]

# gaussianity
data_group1 = dataset_num[i1,]
data_group2 = dataset_num[i2,]
data_group3 = dataset_num[i3,]

P = c( mcshapiro.test(data_group1)$pvalue, 
       mcshapiro.test(data_group2)$pvalue,
       mcshapiro.test(data_group3)$pvalue)
g = 3

# homogeneity
# multivariate
S  <-  cov(dataset_num)
S1 <-  cov(dataset_num[i1,])   # S nel gruppo i = 1..g
S2 <-  cov(dataset_num[i2,])
S3 <-  cov(dataset_num[i3,])   # g righe

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)
round(S3,digits=1)
# to compare visually
x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))

library(MASS)
# QDA

# se voglio settare i prior:
prior.c = c(3/24, 12/24, 9/24)
qda.m <- qda(dataset_num, group, prior = prior.c)
qda.m



# PLOT

x11()
plot(dataset_num, main='Partition', xlab='V1', ylab='V2', pch=20)
points(dataset_num[i1,], col='red', pch=20)   
# se i dati sono in dataframe separati, ad esempio il  gruppo uno è nel dataset g1
# allora semplicemente: dataset_num[i1,] = g1
points(dataset_num[i2,], col='green', pch=20)
points(dataset_num[i3,], col='blue', pch=20)
legend("topright", legend=levels(group), fill=c('red','green','blue'), cex=.7)

points(qda.m$means, pch=4,col=c('red','green','blue') , lwd=2, cex=1.5)
# the x are the centroids

# i create some fine grid in order to plot the regions that we have obtained 
x  <- seq(min(dataset_num[,1]), max(dataset_num[,1]), length=200)
y  <- seq(min(dataset_num[,2]), max(dataset_num[,2]), length=200)

#                nome Var1   nome Var2
xy <- expand.grid( accel = x , gyro =y)

z  <- predict(qda.m, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2], z[,3])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1], z[,3])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
z3 <- z[,3] - pmax(z[,1], z[,2])  # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
# the function contour with level=0 plots the line 
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)




# b) APER
qda.training <- predict(qda.m, dataset_num)  #prediciton on my original data
names(qda.training)

qda.training$class   # assigned classes
group     # true labels
table(class.true=group, class.assigned=qda.training$class)

3/150 * 0.125 + 5/150*0.5 + 10/150* 0.375 

# c)
x_new = c(0.45,0.52)  # vettore in Rp

C1_pred <- predict(qda.m, x_new) 
C1_pred


# d)
k = 5 
x  <- seq(min(dataset_num[,1]), max(dataset_num[,1]), length=200)
y  <- seq(min(dataset_num[,2]), max(dataset_num[,2]), length=200)
xy <- expand.grid(accel=x, gyro=y)

library(class)
dataset.knn <- knn(train = dataset_num, test = xy, cl = group, k = k)


z  <- as.numeric(dataset.knn)

x11()
plot(dataset_num, main='activity', xlab='accel', ylab='gyro', pch=20)
points(dataset_num[i1,], col=2, pch=20)
points(dataset_num[i2,], col=3, pch=20)
points(dataset_num[i3,], col=4, pch=20)
legend("topright", legend=levels(group), fill=c(2,3,4))
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

# come calcolo l'error rate?



#general command
knn.m <- knn(train = dataset_num, test = dataset_num, cl = group, k = 5, prob=T )


#AER (by cross validation) (actual error rate)
misc <- table(class.true = group, class.assigned = knn.m)
# secondo me leggi da qi i misclassificati
misc




xy <- expand.grid(x= 0.4, y= 0.22)
knn.m <- knn(train = dataset_num, test = xy, cl = group, k = 5, prob=T)
knn.m











