### Ex1

library(MASS)

dataset = read.table("occupancy.txt",header=T)

group = factor(dataset$X, levels=c('0', '1'))   # Species è il nome della colonna del grouping nel dataset
# mettendo i levels sono sicura che l'ordine è quello che penso io

g = 2 # number of groups 


i1 <- which(group=='0')       # rows of group1
i2 <- which(group=='1')
   # tante righe quanti gruppi

n1 <- length(i1)
n2 <- length(i2)

n <- n1+n2

dataset_num <- dataset[,1:2]

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


# to compare visually
x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))


# se voglio settare i prior:
prior.c = c(15/24, 9/24)
qda.m <- qda(dataset_num, group, prior = prior.c)


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
xy <- expand.grid(Humidity=x, CO2=y)


z  <- predict(qda.m, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - z[,2] # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - z[,1] # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    


# Plot the contour line of level (levels=0) of z1, z2, z3: 
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j 
# where j realizes the max.
# the function contour with level=0 plots the line 
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)






# 1) Compute the APER 
qda.training <- predict(qda.m, dataset_num)  #prediciton on my original data
names(qda.training)

qda.training$class   # assigned classes
group     # true labels
misc = table(class.true=group, class.assigned=qda.training$class)



G <- 2 #number of group
APER <- 0
for(g in 1:G)
  APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior.c[g] 

APER

APER = 5/60*(9/24) + 1/40 * (15/24)
# 0 = F, 1=V
APER

x_new = c(26,9)  # vettore in Rp

C1_pred <- predict(qda.m, x_new) 
C1_pred


# d) 
library(class)

k = 5
x  <- seq(min(dataset_num[,1]), max(dataset_num[,1]), length=200)
y  <- seq(min(dataset_num[,2]), max(dataset_num[,2]), length=200)
xy <- expand.grid(Humidity=x, CO2=y)
# Sepal.Length è il nome della prima variabile in dataset
# se invece ho più file dove ho diversi gruppi vedi su lda come unirli

## KNN
dataset.knn <- knn(train = dataset_num, test = xy, cl = group, k = k)

z  <- as.numeric(dataset.knn)

x11()
plot(dataset_num, main='Occupancy', xlab='Humidity', ylab='CO2', pch=20)
points(dataset_num[i1,], col=2, pch=20)
points(dataset_num[i2,], col=3, pch=20)

legend("topright", legend=levels(group), fill=c(2,3,4))
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)


dim(dataset)

# domanda error rate

# knn sul training
k = 5

knn.m <- knn(train = dataset_num, test = dataset_num, cl = group, k = k , prob=T)


# Table
n = dim(dataset)[1]
misc <- table(class.true = group, class.assigned = knn.m)

# APER somma i misclassificati e dividi per n
APER = 0
for(i in 1:g)
  APER = APER + sum(misc[i,-i])     # no priors
APER = APER / n 

APER




