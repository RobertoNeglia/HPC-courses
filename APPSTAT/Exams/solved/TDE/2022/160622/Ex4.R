### Ex4

library(fda)
library(fields)


### Visualize data

dataset = read.table("listening.txt", header = T)

dataset = t(dataset)

x11()
matplot(dataset,type='l',main='data',xlab='Abscissa',ylab='Y')

abscissa <- 1:365
n_functions = dim(dataset)[2]

m = 3
nbasis = 100
lambda = 100

basis <- create.bspline.basis(rangeval=range(abscissa),norder = m,nbasis= nbasis)
data.fd <- Data2fd(y = dataset, argvals = abscissa, basisobj = basis, lambda= lambda)

data.fd$coefs

x11()
plot.fd(data.fd, main="B-splines")


## -- WHY????

### FPCA
#                            nharm = number of PC to compute 
pca_result <- pca.fd(data.fd,nharm=5,centerfns=TRUE)

# variance explained by each FPC
pca_result$varprop

# Scree-plot
# Eigenvalues: only the first N-1 are non-null
x11()
plot(pca_result$values[1:n_functions],xlab='j',ylab='Eigenvalues')

# Explained variability
x11()
plot(cumsum(pca_result$values)[1:n_functions]/sum(pca_result$values),xlab='j',ylab='CPV',ylim=c(0.7,1))


n_PC = 3  # number of PC i want to plot
x11()
par(mfrow=c(1,n_PC))   
plot(pca_result$harmonics[1,],col=1,ylab='FPC1')  # se serve setta ylim=c(.., ..)
plot(pca_result$harmonics[2,],col=2,ylab='FPC2')
plot(pca_result$harmonics[3,],col=3,ylab='FPC3')

n_PC = 3  # number of PC i want to plot
x11()
par(mfrow=c(1,n_PC))   
#                                           harm = c(indici delle PC che voglio plottare)                  
plot(pca_result, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)




# Scatter plot of the scores
x11()

plot(pca_result$scores[,1],pca_result$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
# se voglio evidenziare il punto degli score associato alla funzione i
# points(pca_result$scores[i,1],pca_result$scores[i,2], col= 'red', lwd=2)
















