## Ex4

library(fda)
library(KernSmooth)

## Visualize the data

dataset = read.table("tide.txt",header=T)

Xobs0 <- dataset$level           # valore
abscissa <- dataset$time         # ascissa
N <- length(abscissa) # number of locations of observations


# Set parameters
m <- 4          # spline order
degree <- m-1    # spline degree = order -1

nbasis = ...     # numero di basi


# We can choose the correct number of basis via GCV
nbasis_range <- 6:30
gcv <- numeric(length(nbasis_range))
for (i in 1:length(nbasis_range)){
  basis_temp <- create.bspline.basis(range(abscissa), nbasis_range[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis_temp)$gcv
}
x11()
par(mfrow=c(1,1))
plot(nbasis_range,gcv)
n_basis_opt = nbasis_range[which.min(gcv)]
abline(v = n_basis_opt, col = 2)
n_basis_opt
# nbasis too low -> Oversmoothing = underfitting
# nbasis too high -> Undersmoothing = overfitting 

# choose 
nbasis <- n_basis_opt  

# Create the basis
basis <- create.bspline.basis(rangeval = range(abscissa), nbasis=nbasis, norder=m)
# If breaks are not provided, equally spaced knots are created
names(basis)

x11()
plot(basis)
# The basis functions have compact support

# Evaluate the basis on the grid of abscissa
basis_eval <- eval.basis(abscissa, basis)   
dim(basis_eval) # number of data x number of basis
head(basis_eval)
# each column is the evaluation of the basis function in the location where we have the data
# we see that basis 6,7,8,9 are 0 in the first nodes and this is because of the support

# Fit via LS
est_coef = lsfit(basis_eval, Xobs0, intercept=FALSE)$coef
est_coef
# remove the intercept since it is already considered 

# predictions in the observation location
Xsp0 <- basis_eval %*% est_coef  


x11()
par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
abline(v=basis$params)
#




# c)


# First derivative
rappincX1 <- (Xobs0[3:N]-Xobs0[1:(N-2)])/(abscissa[3:N]-abscissa[1:(N-2)])

Xsp1 <- eval.fd(abscissa, Xsp$fd, Lfd=1) # first derivative

x11()
plot(abscissa[2:(N-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1,type='l',col="orange",lwd=3)


# First derivative (argument Lfdobj=1)

# we compute the derivative of the basis and since our function is a lin combination 
# of the basis and we have computed the coeff its derivative will be a 
#lin combination of the derivative of the basis with same  coeff
basis_eval_1<- eval.basis(abscissa, basis, Lfdobj=1)  # =1 first derivative
Xsp1 <- basis_eval_1 %*% est_coef

# to obtain the second derivative (argument Lfdobj=2)
basis_eval_2<- eval.basis(abscissa, basis, Lfdobj=2)
Xsp2 <- basis_eval_2 %*% est_coef

x11(width = 14)
par(mfrow=c(1,2))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
# points(abscissa,X0_vera ,type="l",col="orange",lwd=3)  # T
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
legend("topleft", legend = c("noisy data","true curve","estimated curve"), col = c("black", "orange","blue"), lwd = c(1,3,2))
plot(abscissa[2:(N-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
# points(abscissa_vera,X1_vera,type='l',col="orange",lwd=3)   # T
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
# plot(abscissa[2:(N-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
# points(abscissa_vera,X2_vera,type='l',col="orange",lwd=3)   # T

x11()
plot(abscissa[2:(N-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
# points(abscissa_vera,X1_vera,type='l',col="orange",lwd=3)   # T
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)


# b)
S <- basis_eval%*%solve(t(basis_eval)%*%basis_eval)%*%t(basis_eval) #projection operator = smoothing operator
df = sum(diag(S))  
sigmahat <- sqrt(sum((Xsp0-Xobs0)^2)/(N-df)) #estimate of sigma
lb <- Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub <- Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))

x11()
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")

















