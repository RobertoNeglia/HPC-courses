### Ex4 



library(fda)
library(fields)


### Visualize data

dataset = read.table("traffic.txt", header = T)


dataset = t(dataset)

x11()
matplot(dataset,type='l',main='data',xlab='Abscissa',ylab='Y')


abscissa = 1:24
n_functions = dim(dataset)[2]


### b-spline basis
m <- 4         # spline order
degree <- m-1    # spline degree = order -1
nbasis = 15   # number of basis

basis <- create.bspline.basis(rangeval=range(abscissa),nbasis= nbasis, norder = m)
#                                                                      se non è specificato l'ordine togli questo input (default)

#                   y = dataset va bene se nel dataset ho solo le valutazioni delle funzioni
#                   se ho anche l'ascissa penso di doverla togliere
data.fd <- Data2fd(y = dataset,argvals = abscissa,basisobj = basis)
x11()
plot.fd(data.fd, main="B-splines")


## Extract the coeff
data.fd$coefs


# b)
pca_result <- pca.fd(data.fd,nharm=5,centerfns=TRUE)

# variance explained by each FPC
pca_result$varprop

# Scree-plot
# Eigenvalues: only the first N-1 are non-null
x11()
plot(pca_result$values[1:n_functions],xlab='j',ylab='Eigenvalues')

x11()
layout(cbind(1,3))
plot(pca_result$harmonics[1,],col=1,ylab='FPC1')  # se serve setta ylim=c(.., ..)
plot(pca_result$harmonics[2,],col=2,ylab='FPC2')


# c) 
x11()
plot(cumsum(pca_result$values)[1:n_functions]/sum(pca_result$values),xlab='j',ylab='CPV',ylim=c(0.8,1))



# Plot of the FPCs as perturbation of the mean
n_PC = 2  # number of PC i want to plot
x11()
par(mfrow=c(1,n_PC))   
#                                           harm = c(indici delle PC che voglio plottare)                  
plot(pca_result, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)
# + : andamento se ho valori positivi di quella FPC
# - : andamento se ho valori positivi di quella FPC



# d)
x = c(1:30)
x11()
plot(x, pca_result$scores[,1],xlab="Day",ylab="Scores FPC1",lwd=2)












