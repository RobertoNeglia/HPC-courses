###--------------------###
### LAB 4 (20/03/2023) ###
###--------------------###

### TOPIC:
### Multivariate normality tests
### Box-Cox transformations
### Tests and confidence regions for the mean of a multivariate Gaussian

library(MVN) # Can need separate installation of "gsl" software library
library(car)
library(mvtnorm)
library(mvnormtest)

quartz.options(height=6.5, width=11, reset=FALSE)

#_______________________________________________________________________________
### Test for multivariate normality
# We generate a sample of n=150 observation from bivariate Gaussian

mu <- c(1,2)
mu
sig <- matrix(c(1,1,1,2), 2)
sig

n   <-  150

set.seed(20230320)
X <- rmvnorm(n, mu, sig)

### mvnormtest package ###

# Multivariate version of the Shapiro-Wilk test
# To be used with sample sizes *between 3 and 5000*
# Takes a matrix where each observation is a column
mshapiro.test(t(X)) # Don't forget to transpose your matrix before !!

# Univariate Q-Q plots
quartz()
par(mfrow=c(1,2))
qqnorm(X[,1], main='QQplot of X.1',xlab='theoretical quantiles', ylab='sample quantiles')
qqline(X[,1])
qqnorm(X[,2], main='QQplot of X.2',xlab='theoretical quantiles', ylab='sample quantiles')
qqline(X[,2])

### MVN package ###

# Henze-Zirkler’s test
result <- mvn(data = X, mvnTest = "hz") # Test performed by default
result$multivariateNormality

# Royston’s test (Shapiro-Wilk test under the hood)
result <- mvn(data = X, mvnTest = "royston")
result$multivariateNormality

# with Q-Q plot of the squared mahalanobis distance over chi-square
quartz()
result <- mvn(data = X, mvnTest = "royston", multivariatePlot = "qq")
result$multivariateNormality

# Example: test of normality for the dataset stiff
stiff <- read.table('stiff.dat')
stiff

quartz()
plot(stiff, asp=1, pch=19)

# Normality of the components
quartz()
par(mfcol=c(2,4))

for(i in 1:4)
{
  hist(stiff[,i], prob=T, main=paste('Histogram of V', i, sep=''), xlab=paste('V', i, sep=''))
  lines(900:2800, dnorm(900:2800,mean(stiff[,i]),sd(stiff[,i])), col='blue', lty=2)
  qqnorm(stiff[,i], main=paste('QQplot of V', i, sep=''))
  qqline(stiff[,i])
  print(shapiro.test(stiff[,i])$p)
}

quartz()
result <- mvn(data = stiff, multivariatePlot = "qq", covariance = FALSE)
result$multivariateNormality

### The data don't seem Gaussian. What can we do?
### Identify clusters 
### Identify (and possibly remove) outliers
### Transform the data (e.g., Box-Cox transformations)
### Work without the Gaussian assumption (e.g., permutation tests)

### Let's try to identify and remove outliers:
M <- colMeans(stiff)
S <- cov(stiff)

d2 <- matrix(mahalanobis(stiff, M, S))

stiff_wo_outliers <- stiff[which(d2 <= 8), ]

result <- mvn(data = stiff_wo_outliers)
result$multivariateNormality

quartz()
plot(stiff_wo_outliers, asp=1, pch=19)

# Normality of the components
quartz()
par(mfcol=c(2,4))

for(i in 1:4)
{
  hist(stiff_wo_outliers[,i], prob=T, main=paste('Histogram of V', i, sep=''), xlab=paste('V', i, sep=''))
  lines(900:2800, dnorm(900:2800,mean(stiff_wo_outliers[,i]),sd(stiff_wo_outliers[,i])), col='blue', lty=2)
  qqnorm(stiff_wo_outliers[,i], main=paste('QQplot of V', i, sep=''))
  qqline(stiff_wo_outliers[,i])
  print(shapiro.test(stiff_wo_outliers[,i])$p)
}

#_______________________________________________________________________________
### Box-Cox transformations

###########################################################################
### Univariate Box-Cox transformation (Johnson-Wichern Chap.4.8)
# The Box-Cox transformations are based on the family of power transformation
# (for x>0)
# x_lambda = (x^lambda-1)/lambda if lambda!=0
#            ln(x)               if lambda==0
# that is continuous in lambda for x>0.
# The parameter lambda is determined as the maximum likelihood solution
# under the Gaussian assumption (see J-W)

# Let's plot the transformations for some values of lambda

box_cox <- function(x,lambda)
{
  if(lambda!=0)
    return((x^lambda-1)/lambda)
  return(log(x))
}

x<-seq(0, 25, by=0.01)

quartz()
plot(x, box_cox(x,0), col='gold', type='l', xlab='x', ylab=expression(x[lambda]),
     ylim=c(-20,30), xlim=c(-5,25), asp=1)
title(main='Box-Cox transformations')
curve(box_cox(x,-1),from=0,to=25,col='blue',add=TRUE)
curve(box_cox(x,1),from=0,to=25,col='red',add=TRUE)
curve(box_cox(x,2),from=0,to=25,col='springgreen',add=TRUE)
points(1,0,pch=19,col='black')
abline(v=0,lty=2,col='grey')
legend('topright',c(expression(paste(lambda,'=-1')),expression(paste(lambda,'=0')),
                    expression(paste(lambda,'=1')),expression(paste(lambda,'=2'))),
       col=c('blue','gold','red','springgreen'),lty=c(1,1,1,1,1))

xx=seq(0.01, 20.01, .05)

par(cex=.5)
points(xx, rep(0,length(xx)), col='grey', pch=19)
points(-.5+rep(0,length(xx)), box_cox(xx, -1), col='blue', pch=19)
points(-.5+rep(-.5,length(xx)), log(xx), col='gold', pch=19)
points(-.5+rep(-1,length(xx)), box_cox(xx,1), col='red', pch=19)
points(-.5+rep(-1.5,length(xx)), box_cox(xx,2), col='springgreen', pch=19)
points(1, 0, pch=19, col='black')

# For lambda<1: observations <1 are "spread", observations >1 are "shrinked"
# For lambda>1: observations <1 are "shrinked", observations >1 are "spread"

graphics.off()

##### A trivial example
n <- 100

set.seed(20032023)
x <- rnorm(n)^2

quartz()
hist(x, col='grey', prob=T, xlab='x', main ='Histogram of X')
points(x, rep(0,n), pch=19)
# we would need to spread the obs with small values and shrink the obs with large values
# => we expect a small lambda (lambda<1)

# Univariate Box-Cox transformation
# We compute the optimal lambda of the univariate Box-Cox transformation
# (command powerTransform of library car)
lambda.x <- powerTransform(x) 
lambda.x
# lambda<1: observations <1 are "spread", observations >1 are "shrinked"

# Transformed sample with the optimal lambda (command bcPower of library car)
bc.x <- bcPower(x, lambda.x$lambda)      # it transforms the data of the first argument through the Box-Cox 
                                         # transformation with lambda given as second argument
hist(bc.x, col='grey', prob=T, main='Histogram of BC(X)', xlab='BC(x)')
points(bc.x, rep(0,n), pch=19)

shapiro.test(x)
shapiro.test(bc.x)

### Multivariate Box-Cox transformation (Johnson-Wichern Cap.4.8)
# Similar to univariate transformation, but jointly on all the variables
rm(x)

##### Simulated Example
b <- read.table('data_sim.txt')
head(b)
dim(b)

attach(b)

quartz()
plot(b, pch=19,main='Data', xlim=c(1,7), ylim=c(-20,800))
points(x, rep(-20,dim(b)[1]), col='red', pch=19)
points(rep(1,dim(b)[1]), y, col='blue', pch=19)

quartz()
par(mfrow=c(2,2))

hist(x, prob=T,col='grey85')
lines(0:1000/100, dnorm(0:1000/100, mean(x), sd(x)), col='blue', lty=2)

hist(y, prob=T,col='grey85')
lines(0:1000, dnorm(0:1000, mean(y), sd(y)), col='blue', lty=2)

qqnorm(x, main='QQplot of x')
qqline(x)

qqnorm(y, main='QQplot of y')
qqline(y)
# For y: left tail too light, right tail too heavy

shapiro.test(x)
shapiro.test(y)

graphics.off()

mshapiro.test(t(b))

### Bivariate Box-Cox transformation
# Compute the optimal lambda of a bivariate Box-Cox transformation
# (command powerTransform, with a multivariate input)
lambda <- powerTransform(cbind(x,y))
lambda
# lambda[1] close to 1 

# lambda[2] close to 0, it is almost a log transform; 
# lambda[2]<1:  obs <1 "spread", obs >1 "shrinked" 
# => it adds weight to the left tail, lightens the right tail

# Compute the transformed data with optimal lambda (of the bivariate transf.)
# (command bcPower)
BC.x <- bcPower(x, lambda$lambda[1])
BC.y <- bcPower(y, lambda$lambda[2])

xx <- seq(0, 7, by=0.01)

quartz()
par(mfrow=c(1,3))

plot(xx, box_cox(x=xx, lambda=lambda$lambda[1]), col='red', lty=1, type='l',
     xlab=expression(x), ylab=expression(x[lambda]), ylim=c(-5,10), asp=1)
title(main='Box-Cox transformation')
lines(xx, box_cox(x=xx, lambda=lambda$lambda[2]), col='blue')
points(1, 0, pch=19, col='black')
abline(a=-1, b=1, lty=2, col='grey')
legend('bottomright', c(expression(lambda[x]),expression(lambda[y]),
                        expression(paste(lambda,'=1'))),
       col=c('blue','red','grey'),lty=c(1,1,1))

plot(b, pch=19, main='Data', xlim=c(1,7))
points(rep(1,200), y, pch=19, col='blue') # projection on y
points(x, rep(-20,200), pch=19, col='red')  # projection on x

plot(BC.x, BC.y, pch=19, main='Bivariate BC', xlim=c(0,5))
points(rep(0,200), BC.y, pch=19, col='blue') # projection on y
points(BC.x, rep(0.7,200), pch=19, col='red')  # projection on x

graphics.off()

# Let's formulate an hypothesis of transformation:
# since we get lambda[1]~1 and lambda[2]~0, we could reasonably consider:
hyp.x <- x
hyp.y <- log(y)

# Remark: from the application-oriented viewpoint, explaining a log
# transform could be simpler than introducing Box-Cox transformations

quartz()
par(mfrow=c(3,3))
plot(b, pch=19, main='Data', xlim=c(1,7))
points(rep(1,200), y, pch=19, col='blue')   # projection on y
points(x, rep(-20,200), pch=19, col='red')  # projection on x

plot(BC.x, BC.y, pch=19, main='BC Bivariate', xlim=c(0,7), ylim=c(0.5,7.5))
points(rep(0,200), BC.y, pch=19, col='blue')   # projection on y
points(BC.x, rep(0.5,200), pch=19, col='red')  # projection on x

plot(hyp.x, hyp.y, pch=19, main='According to hyp.', xlim=c(0,7), ylim=c(0.5,7.5))
points(rep(0,200), hyp.y, pch=19, col='blue')  # projection on y
points(hyp.x, rep(0.5,200), pch=19, col='red') # projection on x

qqnorm(x, main="x", col='red')
qqnorm(BC.x, main="BC.x", col='red')
qqnorm(hyp.x, main="hyp.x", col='red')

qqnorm(y, main="y", col='blue')
qqnorm(BC.y, main="BC.y", col='blue')
qqnorm(hyp.y, main="hyp.y", col='blue')

dev.off()

# Univariate Shapiro-Wilk (H0: univariate gaussianity along a given direction)
shapiro.test(x)$p
shapiro.test(BC.x)$p
shapiro.test(hyp.x)$p

shapiro.test(y)$p
shapiro.test(BC.y)$p
shapiro.test(hyp.y)$p

# MC-Shapiro test (H0: multivariate gaussianity)
mshapiro.test(rbind(x,y))$p
mshapiro.test(rbind(BC.x,BC.y))$p
mshapiro.test(rbind(hyp.x,hyp.y))$p

# High p-value for the transformed data: with the obtained transformation we
# have no evidence of non-normality, i.e., we can assume that the transformed
# data are Gaussian

detach(b)

graphics.off()

### Box-Cox transformation on the dataset stiff

stiff <- read.table('stiff.dat')
head(stiff)
dim(stiff)

lambda.mult <- powerTransform(stiff)    
lambda.mult

BC.x <- bcPower(stiff[,1], lambda.mult$lambda[1]) 
BC.y <- bcPower(stiff[,2], lambda.mult$lambda[2])
BC.z <- bcPower(stiff[,3], lambda.mult$lambda[3]) 
BC.w <- bcPower(stiff[,4], lambda.mult$lambda[4])

# Plot of the original variables
attach(stiff)
quartz()
par(mfrow=c(2,4))
xx<-seq(1320, 3000, length=100)
hist(V1, prob=T, breaks=8, col='grey85')
lines(xx, dnorm(xx,mean(V1), sd(V1)), col='blue', lty=2)
yy<-seq(1150, 2800, length=100)
hist(V2, prob=T,breaks=8,col='grey85')
lines(yy, dnorm(yy,mean(V2), sd(V2)), col='blue', lty=2)
zz<-seq(1000, 2420, length=100)
hist(V3, prob=T, breaks=8, col='grey85')
lines(zz, dnorm(zz,mean(V3), sd(V3)), col='blue', lty=2)
ww<-seq(1100, 2600, length=100)
hist(V4, prob=T, breaks=8, col='grey85')
lines(ww, dnorm(ww, mean(V4), sd(V4)), col='blue', lty=2)

qqnorm(V1, main='QQplot of V1')
qqline(V1)

qqnorm(V2, main='QQplot of V2')
qqline(V2)

qqnorm(V3, main='QQplot of V3')
qqline(V3)

qqnorm(V4, main='QQplot of V4')
qqline(V4)

detach(stiff)

# Plot of the transformed variables
quartz()
par(mfrow=c(2,4))
xx<-seq(10, 12, length=100)
hist(BC.x, prob=T, breaks=8, col='grey85')
lines(xx, dnorm(xx, mean(BC.x), sd(BC.x)), col='blue', lty=2)
yy<-seq(3, 3.2, length=100)
hist(BC.y, prob=T, breaks=8, col='grey85')
lines(yy, dnorm(yy, mean(BC.y), sd(BC.y)), col='blue', lty=2)
zz<-seq(12, 15, length=100)
hist(BC.z, prob=T, breaks=8, col='grey85')
lines(zz, dnorm(zz, mean(BC.z), sd(BC.z)), col='blue', lty=2)
ww<-seq(250, 500, length=100)
hist(BC.w, prob=T, breaks=8, col='grey85')
lines(ww, dnorm(ww, mean(BC.w), sd(BC.w)), col='blue', lty=2)

qqnorm(BC.x, main='QQplot of BC.x')
qqline(BC.x)

qqnorm(BC.y, main='QQplot of BC.y')
qqline(BC.y)

qqnorm(BC.z, main='QQplot of BC.z')
qqline(BC.z)

qqnorm(BC.w, main='QQplot of BC.w')
qqline(BC.w)

# One transformed variable at a time
shapiro.test(BC.x)
shapiro.test(BC.y)
shapiro.test(BC.z)
shapiro.test(BC.w)

# All together
mshapiro.test(rbind(BC.x, BC.y, BC.z, BC.w))

# Note: the higher the dimensionality, the more difficult to 
# recover the gaussianity

# In this case, Box-Cox transformations are not enough. We saw that removing
# outliers indeed solve the problems of non-Gaussianity.
# Alternative approaches:
# - Asymptotics
# - Non-parametrics (e.g., permutation tests)

graphics.off()



#_______________________________________________________________________________
### Tests and confidence regions for the mean of a multivariate Gaussian

### Example with simulated data from a bivariate Gaussian

mu <- c(1, 0)
sig <- matrix(c(1, 1, 1, 2), nrow=2)

set.seed(123)
x <- rmvnorm(n=30, mean=mu, sigma=sig)
x <- data.frame(X.1=x[,1], X.2=x[,2])

quartz()
plot(x, asp = 1, pch = 19)

dev.off()

n <- dim(x)[1]
p <- dim(x)[2]

x.mean   <- sapply(x, mean)
x.cov    <- cov(x)
x.invcov <- solve(x.cov)

#_______________________________________________________________________________
### Test for the mean of a multivariate Gaussian 

################################################################################
# Premise: general rule to perform a test
# 1)  Formulate the test (and test the Gaussian assumption, if needed)
# 2)  Compute the test statistics 
# 3a) Having set the level of the test, verify whether the test statistics 
#     belongs to the region of rejection (i.e., if there is statistical  
#     evidence to reject H0)
# 3b) Compute the p-value of the test
################################################################################

#_______________________________________________________________________________
### Test on the mean of level alpha=1%
### H0: mu == mu0 vs H1: mu != mu0
### with mu0 = c(1, 0)
###-----------------------------------
mshapiro.test(t(x))

alpha <- 0.01
mu0 <- c(1,0)

# T2 Statistics
x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
# Test: 
x.T2 < cfr.fisher   # no statistical evidence to reject H0 at level alpha
# Rejection region: {x.T2 > cfr.fisher}
# (we reject for large values of the T2 statistics)

# Compute the p-value 
P <- 1 - pf(x.T2*(n-p) / ((n-1)*p), p, n-p)
P

quartz()
xx <- seq(0,40,by=0.05)
plot(xx, df(xx*(n-p)/((n-1)*p),p,n-p), type="l", lwd=2, main='Density F(p,n-p)',
     xlab='x*(n-p)/((n-1)*p)', ylab='Density')
abline(h=0, v=x.T2*(n-p)/((n-1)*p), col=c('grey','red'), lwd=2, lty=c(2,1))
# The P-value is high because the test statistics is central with respect to 
# its distribution under H0
# => we cannot reject for any reasonable level (we would reject for a level alpha>81.7%)

dev.off()

# Region of rejection (centered in mu0)
quartz()
plot(x, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)

# We add a red point in correspondence of the sample mean
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)

# Confidence region (centered in x.mean)
# { m \in R^2 s.t. n * (x.mean-m)' %*% (x.cov)^-1 %*% (x.mean-m) < cfr.fisher }

ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)

# Remark: the radius and the shape of the ellipse are the same, but the center changes:
# - Rejection region: the center is the mean mu0 under H0 (blue ellipse)
# - Confidence region: the center is the sample mean (red ellipse)

# Which relation between the two ellipses?
# - If the rejection region does NOT contain the sample mean (i.e., we
#   are in the acceptance region), then we cannot reject H0 
#   (i.e., if the sample mean falls within the ellipse we accept H0)
# - If the mean under H0 (mu0) is contained in the confidence region
#   of level 1-alpha, then we do not reject H0 at level alpha
# => the confidence region of level 1-alpha contains all the mu0
#    that we would accept at level alpha

# Note: by definition, the confidence region of level 1-alpha
# produces ellipsoidal regions that contain the true mean
# 100(1-alpha)% of the times. If H0 is true (i.e., mu0 is 
# the true mean), those ellipsoidal regions will contain mu0 
# 100(1-alpha)% of the times

#_______________________________________________________________________________
### Asymptotic test on the mean
### H0: mu == mu0 vs H1: mu != mu0
### with mu0 = c(1,0)
###---------------------------------------------

# Note: we don't need to verify the Gaussian assumption!
# Warning: we are going to use an asymptotic test, but we only have n = 30 data!

mu0   <- c(1,0)

x.T2A   <- n * (x.mean-mu0) %*%  x.invcov  %*% (x.mean-mu0)
cfr.chisq <- qchisq(1-alpha, p)
x.T2A < cfr.chisq # no statistical evidence to reject H0 at level alpha

# Compute the p-value
PA <- 1 - pchisq(x.T2A, p)
PA

quartz()
par(mfrow=c(1,2))
plot(x, asp = 1,main='Comparison rejection regions')
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 1,
        center.pch = 4, center.cex=1.5, lwd=2)
ellipse(mu0, x.cov/n, sqrt(cfr.chisq), col = 'lightblue', lty = 1,
        center.pch = 4, center.cex=1.5, lwd=2)
points(mu0[1], mu0[2], pch = 4, cex = 1.5, lwd = 2, col ='lightblue')
legend('topleft', c('Exact', 'Asymptotic'), col=c('blue', 'lightblue'),
       lty=c(1), lwd=2)

plot(x, asp = 1,main='Comparison of confidence regions')
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 1, center.pch = 4,
        center.cex=1.5, lwd=2)
ellipse(x.mean, x.cov/n, sqrt(cfr.chisq), col = 'orange', lty = 1, center.pch = 4,
        center.cex=1.5, lwd=2)
points(x.mean[1], x.mean[2], pch = 4, cex = 1.5, lwd = 2, col ='orange')
legend('topleft', c('Exact', 'Asymptotic'), col=c('red','orange'),
       lty=c(1), lwd=2)

graphics.off()

#_______________________________________________________________________________
###  We change the null Hypothesis  
###-------------------------------
mu0 <- c(1.5,-0.5)
mu0

mu

# Remark: now H0 is false (we will want to reject H0)

### Test of level 1%
### H0: mu == mu0=c(1.5,-0.5) vs H1: mu != mu0=c(1.5,-0.5) 
###--------------------------------------------------

### Ellipsoidal region
x.T2 <- n * t(x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
x.T2 < cfr.fisher

# Compute the p-value
P <- 1 - pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P

quartz()
plot(x, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col = 'red', cex=1.5)

dev.off()

# How would we communicate the results?
# For instance, we could provide confidence intervals along interesting
# directions (if they don't contain the mean under H0, it means that
# we reject the associated univariate test along that direction)

# Let's try with simultaneous T2 confidence intervals on the coordinate directions:
# Recall: these are projections of the ellipsoidal confidence region

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2
# Both the intervals contain the mean under H0
# (i.e., mu0 is contained in the rectangular region determined by
# the projection of the ellipsoid along the coordinate directions)

# Remark: this is not in contrast with the previous findings
# Rejecting the global T2-test means that we reject H0 along at least one
# direction, not necessarily along the coordinate direction

quartz()
par(mfrow=c(1,1))

# Plot of the acceptance region and the sample mean
plot(x, asp = 1,main='Confidence and rejection regions')
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
points(x.mean[1], x.mean[2], pch = 16, col = 'red', cex=1.5)

# Plot of the confidence region and the region defined by the cartesian product
# of the sim. confidence intervals in direction of x1 and x2
ellipse(x.mean, shape=x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, center.pch = 16)
rect(T2[1,1], T2[2,1], T2[1,3], T2[2,3], border='red', lwd=2)

# Let's try with Bonferroni intervals
k <- p # number of intervals I want to compute (set in advance)
cfr.t <- qt(1-alpha/(2*k), n-1)
Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf
# Both the intervals contain the mean under H0
# (i.e., mu0 is contained in the rectangular region determined by
# the Bonferroni intervals along the coordinate directions)

# we add the Bonferroni intervals to the plot
rect(Bf[1,1], Bf[2,1], Bf[1,3], Bf[2,3], border='orange', lwd=2)
legend('topleft', c('Rej. Reg.', 'Conf. Reg','T2-sim', 'Bonferroni'),
       col=c('blue','red','red','orange'), lty=c(2,2,1,1), lwd=2)

# Remark: if we wanted to compute additional Bonferroni intervals
# along other directions, we would need to re-compute all the Bonferroni
# intervals with another correction k

graphics.off()

#_______________________________________________________________________________
### Example: analysis of a real dataset (dataset stiff)

stiff <- read.table('stiff.dat')
head(stiff)
dim(stiff)

n <- dim(stiff)[1]
p <- dim(stiff)[2]

quartz()
plot(stiff,pch=19)

dev.off()

### Normality test
mvn(stiff)$multivariateNormality
mvn(stiff_wo_outliers)$multivariateNormality

n <- dim(stiff_wo_outliers)[1]
p <- dim(stiff_wo_outliers)[2]

### Test for the mean of level 5%
###--------------------------------
# 1) Formulate the test:
#    H0: mu == mu0 vs H1: mu != mu0
#    with mu0 = c(1850, 1750, 1500, 1700)

mu0   <- c(1850, 1750, 1500, 1700)
alpha <- 0.05

# 2) Compute the test statistics
x.mean   <- colMeans(stiff_wo_outliers)
x.cov    <- cov(stiff_wo_outliers)
x.invcov <- solve(x.cov)

x.T2     <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 

# 3a) Verify if the test statistics belongs to the rejection region
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
x.T2 < cfr.fisher # we accept H0 at 5%

# 3b) Compute the p-value
P <- 1 - pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P

### Confidence region for the mean of level 95%
###------------------------------------------------------------
# 1)  Identify the type of region of interest:
#     CR for the mean (ellipsoidal region) 
#     {m \in R^4 s.t. n * (x.mean-m)' %*% x.invcov %*% (x.mean-m) < cfr.fisher}
# 2)  Characterize the region: compute the center, direction of 
#     the principal axes, length of the axes

# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values)

# Question: how to plot the confidence region? 
# (I don't know how to plot in R^4!)
# I can work with 'representations' of the confidence regions
# (e.g., projections along particular directions)

# We plot the projections of the ellipsoid in some directions 
# of interest (e.g. the X1,...,Xp coordinates)
# => We plot the simultaneous T2 confidence intervals in each
#    direction of interest (with global coverage alpha)

### Simultaneous T2 intervals on the components of the mean
### with global level 95%
###----------------------------------------------------

T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2

quartz()

# Plot of the simultaneous T2 intervals in the direction of X1,...,X4
matplot(1:4,1:4, pch='',ylim=range(stiff_wo_outliers), xlab='Variables', ylab='T2 for a component', 
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:4) segments(i, T2[i,1], i, T2[i,3], lwd=3, col=i)
points(1:4, T2[,2], pch=16, col=1:4)

# Is mu0 inside the rectangular region?
# We add it to the plot
points(1:4, mu0, lwd=3, col='orange')

# Yes, it is, because it is inside all the T2-intervals,

### Bonferroni intervals on the components of the mean
### with global level 95%
###----------------------------------------------------
k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf

# Let's do a plot
quartz()
matplot(1:4, 1:4, pch='', ylim=range(stiff_wo_outliers), xlab='Variables',
        ylab='Confidence intervals along a component', main='Confidence intervals')

for(i in 1:4) segments(i, T2[i,1], i, T2[i,3], lwd=2, col='grey35', lty=3)
points(1:4, T2[,1], pch='-', col='grey35')
points(1:4, T2[,3], pch='-', col='grey35')

for(i in 1:4) segments(i, Bf[i,1], i, Bf[i,3], lwd=2, col=i)
points(1:4, Bf[,2], pch=16, col=1:4)
points(1:4, Bf[,1], pch='-', col=1:4)
points(1:4, Bf[,3], pch='-', col=1:4)

# Is mu0 inside the Bonferroni confidence region?
# we add it to the plot
points(1:4, mu0, lwd=3, col='orange')

# Yes, it is, because it belongs to all the intervals along the components

graphics.off()
