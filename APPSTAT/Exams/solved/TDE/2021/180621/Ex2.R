### Ex2

dataset_tot <- read.table('beans.txt', header=T)

dataset = dataset_tot[, -9]



n <- dim(dataset)[1]
p <- dim(dataset)[2]

# Boxplot centered in the means
x11()
boxplot(scale(x=dataset,center = T, scale=F), las=2, col='gold')

# scale
dataset <- scale(dataset)  # the default command is with center=T and scale= T 
dataset <- data.frame(dataset)

x11()
boxplot(dataset, las=2, col='gold')


# PCA
pc.dataset <- princomp(dataset, scores=T)
pc.dataset
summary(pc.dataset)


# Explained variance
max1 = max(pc.dataset$scores)   # va bene così?
max2 = max(sapply(dataset,sd)^2)
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.dataset, las=2, main='Principal components', ylim=c(0,max1))
abline(h=1, col='blue')
barplot(sapply(dataset,sd)^2, las=2, main='Original Variables', ylim=c(0,max2), ylab='Variances')
plot(cumsum(pc.dataset$sd^2)/sum(pc.dataset$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')  #0.8 is a good treshold
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(dataset),labels=1:ncol(dataset),las=2)




# b) 
# plot only the first k PCs
load <- pc.dataset$loadings
k = 2
x11()
par(mar = c(1,k,0,2), mfrow = c(k,1))
for(i in 1:k) barplot(load[,i], ylim = c(-1, 1))




scores.dataset <- pc.dataset$scores
# scores.dataset

x11()
plot(scores.dataset[which(dataset_tot$Type == 'cannellini'),1:2], col= 'red', xlim=c(-5,5))
abline(h=0, v=0, lty=2, col='grey')
points(scores.dataset[which(dataset_tot$Type == 'adzuki'),1:2], col='green', xlim=c(-5,5))
points(scores.dataset[which(dataset_tot$Type == 'black-eyed'),1:2], col='blue', xlim=c(-5,5))
abline(h=0, v=0, lty=2, col='grey')


# c)
data_cann = data.frame(scores.dataset[which(dataset_tot$Type == 'cannellini'),1:2])


n <- 50
p <- 2
# test of Gaussianity
load("mcshapiro.test.RData")
mcshapiro.test(data_cann)$pvalue
x.mean <- colMeans(data_cann)
x.cov <- cov(data_cann)


# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

alpha = 0.05
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values) 

library(car)
x11()
plot(data_cann, asp = 1)
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)










