### Ex4

library(fda)
library(KernSmooth)

dataset = read.table("power.txt",header=T)

Xobs0 <- dataset$power           # valore
abscissa <-  1:365    # ascissa
N <- length(abscissa) # number of locations of observations

# generalized cross-validation
nbasis <- 6:50
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.fourier.basis(range(abscissa), nbasis[i])
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis)$gcv
}
x11()
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis_opt = nbasis[which.min(gcv)]
abline(v=nbasis_opt,col='red')

nbasis = 12
# create the basis
basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis)

x11()
plot(basis)

# Evaluate the basis on the grid of abscissa
Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)

# predictions in the observation location
Xsp0 <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="red",lwd=2)


# b)
rappincX1 <- (Xobs0[3:N]-Xobs0[1:(N-2)])/(abscissa[3:N]-abscissa[1:(N-2)])

x11()
par(mfrow=c(1,3))
plot(abscissa,Xobs0,xlab="t",ylab="observed data",type="l")


# First derivative
Xsp1 <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative

x11()
plot(abscissa[2:(N-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1,type='l',col="orange",lwd=3)


# Oversmoo = underfit
nbasis_over = 3

basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis_over)


# Evaluate the basis on the grid of abscissa
Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)

# predictions in the observation location
Xsp0 <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="red",lwd=2)



# Undersmooth

nbasis_over = 40

basis <- create.fourier.basis(rangeval=range(abscissa), nbasis=nbasis_over)


# Evaluate the basis on the grid of abscissa
Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)

# predictions in the observation location
Xsp0 <- eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="red",lwd=2)
























