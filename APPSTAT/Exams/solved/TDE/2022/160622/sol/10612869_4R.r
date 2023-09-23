# The file listening.txt contains daily measurements of music listening of 129 different songs along 2021. Consider
# a functional data analysis approach where, for each song, the measurements provided are considered as discrete
# sampling of underlying smooth functions.
# a) Perform a smoothing of each datum through cubic smoothing splines (order = 3), with a smoothing parameter
# λ = 102 . Specify the choice for the nodes of the splines. Provide a plot of the smoothed data and report the
# first 3 coefficients obtained for the first song.
# b) Perform a functional principal component analysis of the smoothed data obtained at point (a). Report the
# variance explained along the first 5 functional principal components and the screeplot.
# c) Propose a possible dimensionality reduction for the data and justify your choice. Plot the retained principal
# components.
# d) Plot the retained principal components as perturbation of the mean, and interpret them.
# e) Provide a plot of the scores along the first two principal components and comment the results.

library(fda)



listening <- read.table('listening.txt', header=T)

# a) Perform a smoothing of each datum through cubic smoothing splines (order = 3), with a smoothing parameter
# λ = 102 . Specify the choice for the nodes of the splines. Provide a plot of the smoothed data and report the
# first 3 coefficients obtained for the first song.

data_W <- t(as.matrix(listening))

time <- 1:365
basis <- create.bspline.basis(rangeval=c(0,365),norder = 3,nbasis=100)
data_W.fd <- Data2fd(y = data_W, argvals = time, basisobj = basis, lambda=102)
plot.fd(data_W.fd)

#first song
data_W.fd$coefs[1:3,1]


# b) Perform a functional principal component analysis of the smoothed data obtained at point (a). Report the
# variance explained along the first 5 functional principal components and the screeplot.

pca_W <- pca.fd(data_W.fd,nharm=5,centerfns=TRUE)

#variance along the first 5 fpc
pca_W$values[1:5] 

#screeplot
plot(cumsum(pca_W$values)/sum(pca_W$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

# c) Propose a possible dimensionality reduction for the data and justify your choice. Plot the retained principal
# components.
cumsum(pca_W$values)/sum(pca_W$values)
layout(cbind(1,2,3))
plot(pca_W$harmonics[1,],col=2,ylab='FPC1',ylim=c(-0.1,0.1))
plot(pca_W$harmonics[2,],col=3,ylab='FPC2',ylim=c(-0.1,0.1))
plot(pca_W$harmonics[3,],col=4,ylab='FPC3',ylim=c(-0.1,0.1))


# d) Plot the retained principal components as perturbation of the mean, and interpret them.
par(mfrow=c(3,1))
plot(pca_W, nx=100, pointplot=TRUE, harm=c(1:3), expand=0, cycle=FALSE)


# e) Provide a plot of the scores along the first two principal components and comment the results.
plot(pca_W$scores[,1],pca_W$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
abline(v=0)
abline(h=0)
