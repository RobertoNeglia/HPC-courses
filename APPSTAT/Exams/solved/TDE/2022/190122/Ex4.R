### Ex4

library(fda)
library(fields)


### Visualize data

dataset = read.table("temperature.txt", header = T)

dataset = dataset[, -366]

dataset = t(dataset)


abscissa = 1:365
N = length(abscissa)
n_functions = dim(dataset)[2]

x11()
matplot(dataset,type='l',main='data',xlab='Abscissa',ylab='Y')


nbasis = 21    #number of basis

basis <- create.fourier.basis(rangeval= range(abscissa),nbasis = nbasis)

#                   y = dataset va bene se nel dataset ho solo le valutazioni delle funzioni
#                   se ho anche l'ascissa penso di doverla togliere
data.fd <- Data2fd(y = dataset ,argvals = abscissa ,basisobj = basis)
x11()
plot.fd(data.fd)


data.fd$coefs


# b) 
pca_result <- pca.fd(data.fd,nharm=5,centerfns=TRUE)

# variance explained by each FPC
pca_result$varprop

# Scree-plot
# Eigenvalues: only the first N-1 are non-null
x11()
plot(pca_result$values[1:n_functions],xlab='j',ylab='Eigenvalues', xlim=c(0,50))

# Explained variability
x11()
plot(cumsum(pca_result$values)[1:n_functions]/sum(pca_result$values),xlab='j',ylab='CPV',ylim=c(0.8,1), xlim=c(0,50))



x11()
layout(cbind(1,2))
plot(pca_result$harmonics[1,],col=1,ylab='FPC1')  # se serve setta ylim=c(.., ..)
plot(pca_result$harmonics[2,],col=2,ylab='FPC2')


# Plot of the FPCs as perturbation of the mean
n_PC = 3  # number of PC i want to plot
x11()
par(mfrow=c(1,n_PC))   
#                                           harm = c(indici delle PC che voglio plottare)                  
plot(pca_result, nx=100, pointplot=TRUE, harm=c(1,2, 3), expand=0, cycle=FALSE)



# scatterplot
DATA =read.table("temperature.txt", header = T)
country = DATA[, 366]
italy = which(country == 'Italy')
germ = which(country == 'Germany')

x11()
plot(pca_result$scores[italy,1],pca_result$scores[italy,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, col='red', xlim=c(-100,100), ylim= c(-80,80))
# se voglio evidenziare il punto degli score associato alla funzione i
# points(pca_result$scores[i,1],pca_result$scores[i,2], col= 'red', lwd=2)
points(pca_result$scores[germ,1],pca_result$scores[germ,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, col='blue')
legend('topright', c('Italy', 'Germany'), pch = 1,col=c('red', 'blue'), cex = 1)












