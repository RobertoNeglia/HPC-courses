#### Ex4

library(fda)
library(fields)


### Visualize data

dataset = read.table("wind.txt", header = T)

dataset = t(dataset)

x11()
matplot(dataset,type='l',main='data',xlab='Abscissa',ylab='Y')


abscissa = 1:24
n_functions = 30

### b-spline basis
nbasis = 12  # number of basis
m = 3
degree = m-1

basis <- create.bspline.basis(rangeval=range(abscissa),nbasis= nbasis, norder = m)

#                   y = dataset va bene se nel dataset ho solo le valutazioni delle funzioni
#                   se ho anche l'ascissa penso di doverla togliere
data.fd <- Data2fd(y = dataset,argvals = abscissa,basisobj = basis)
x11()
plot.fd(data.fd, main="B-splines")

data.fd$coefs

pca_result <- pca.fd(data.fd,nharm=5,centerfns=TRUE)


# scree - plot
x11()
plot(pca_result$values[1:n_functions],xlab='j',ylab='Eigenvalues')



# Explained variability
x11()
plot(cumsum(pca_result$values)[1:n_functions]/sum(pca_result$values),xlab='j',ylab='CPV',ylim=c(0.5,1))

pca_result$varprop

x11()
layout(cbind(1,2))
plot(pca_result$harmonics[1,],col=1,ylab='FPC1')  # se serve setta ylim=c(.., ..)
plot(pca_result$harmonics[2,],col=2,ylab='FPC2')
# plot(pca_result$harmonics[3,],col=3,ylab='FPC3')  -> non va 

# Plot of the FPCs as perturbation of the mean
n_PC = 3  # number of PC i want to plot
x11()
par(mfrow=c(1,n_PC))   
#                                           harm = c(indici delle PC che voglio plottare)                  
plot(pca_result, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)

# Scatter plot of the scores
x11()
par(mfrow=c(1,2))
plot(pca_result$scores[,1],pca_result$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
plot(pca_result$scores[1,1],pca_result$scores[1,2], col= 'red', lwd=2)
plot(pca_result$scores[,1],pca_result$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")

text(pca_result$scores[,1],pca_result$scores[,2],dimnames(dataset)[[2]], cex=1)

x11()
par(mfrow=c(1,2))
plot(pca_result$scores[,1],pca_result$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca_result$scores[1,1],pca_result$scores[1,2], col= 'red', lwd=2)

plot(pca_result$scores[,1],pca_result$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_result$scores[,1],pca_result$scores[,2],dimnames(dataset)[[2]], cex=1)




















