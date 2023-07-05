#_______________________________________________________________________________
### Further material
### Simulation to compute the actual level and the power of the T2 test

library(mvtnorm)

### 1) Exact test (T2) & H0 TRUE
#-----------------------------------------
mu0 <- c(1,0)
alpha <- 0.01

mu <- c(1,0)
sig <- matrix(c(1,1,1,2),nrow=2)
n=100
p=2

T2.stat <- NULL
# We generate 10000 samples under H0 and
# count how many times we reject H0 (true)
# [this is an estimate of the actual level]
for(i in 1:10000)
{
  x <- rmvnorm(n, mean=mu, sigma=sig)
  x.mean   <- colMeans(x)
  x.cov    <- cov(x)
  x.invcov <- solve(x.cov)
  x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
  T2.stat <- c(T2.stat, x.T2)
}
# The vector T2.stat collects 10000 independent realizations
# of the T2 statistics, computed under H0
# => we know its distribution!

x11()
layout(cbind(c(1,1,1),c(2,3,4)), widths=c(2,1), heights=c(1,1))
plot(T2.stat, xlab='Simulation', ylab='T2 Statistics',xlim=c(0,10500))
abline(h = ((n-1)*p/(n-p))*qf(1-alpha,p,n-p), col='red', lwd=2)  # criterion to rejecto H0
points(which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p)),T2.stat[which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))], col='orange')
points(rep(10500,10000),T2.stat, xlab='Simulations', ylab='T2 statistics',pch=19)
points(rep(10500, length(which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p)))),T2.stat[which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))], pch=19, col='orange')
title(main='T2 Statistics (Exact Test, H0 true)')
legend('topleft',c('Accept H0', 'Reject H0', 'cfr.fisher'), col=c('black','orange','red'),pch=c(19,19,-1),lty=c(-1,-1,1),lwd=c(1,1,2),cex=.8)

plot(0,0, ylab='Simulations', xlab='',ylim=c(0,10000),main='Number of times we reject H0')
count=0
for(i in 1:10000)
{
  if(T2.stat[i]>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))
  {
    count=count+1
    points(0,count,col='orange',pch=19)
  }
}
abline(h=10000*alpha,lty=2,col='grey')
legend('topleft', c('(Empirical level)*10000', '(Nominal level)*10000'), col=c('orange','grey'), lty=c(-1,2),pch=c(19,-1),cex=.8)

# we compute the proportion of times on the total number
# of simulations in which I reject H0:
# this is an empirical estimate of the level of the test
# (recall: H0 is true!)
sum(T2.stat > ((n-1)*p/(n-p))*qf(1-alpha,p,n-p))/10000

# we test the distribution of T2 under H0:
# (it should be proportional to a Fisher)
hist(T2.stat, prob=TRUE,col='grey85',ylim=c(0,1),xlab='T2.stat',ylab='Density',main='Histogram of T2.stat')
xx <- seq(0,40,by=0.05)
lines(xx,df(xx*(n-p)/((n-1)*p),p,n-p),type="l",lwd=2)
abline(v = ((n-1)*p/(n-p))*qf(1-alpha,p,n-p), col='red', lwd=2)

quantiles.H0 <- ((n-1)*p/(n-p))*qf((1:10000 - 0.5)/10000,p,n-p)
qqplot(quantiles.H0, T2.stat, asp=1,main='QQ-plot')
abline(0,1)

### 2) Asymptotic test (Chi-squared) & H0 TRUE
#----------------------------------------------------
# I repeat the same experiment but with an
# asymptotic test
T2.stat <- NULL
for(i in 1:10000)
{
  x <- rmvnorm(n, mean=mu, sigma=sig)
  x.mean   <- colMeans(x)
  x.cov    <- cov(x)
  x.invcov <- solve(x.cov)
  x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
  T2.stat <- c(T2.stat, x.T2)
}

x11()
layout(cbind(c(1,1,1),c(2,3,4)), widths=c(2,1), heights=c(1,1))
plot(T2.stat, xlab='Simulation', ylab='T2 statistics',xlim=c(0,10500))
abline(h = qchisq(1-alpha,p), col='red', lwd=2)  # criterion to reject
points(which(T2.stat>qchisq(1-alpha,p)),T2.stat[which(T2.stat>qchisq(1-alpha,p))], col='orange')
points(rep(10500,10000),T2.stat, xlab='Simulations', ylab='T2 Statistics',pch=19)
points(rep(10500, length(which(T2.stat>qchisq(1-alpha,p)))),T2.stat[which(T2.stat>qchisq(1-alpha,p))], pch=19, col='orange')
title(main='T2 Statistics (Asymptotic test, H0 true)')
legend('topleft',c('Accept H0', 'Reject H0', 'cfr.chisq'), col=c('black','orange','red'),pch=c(19,19,-1),lty=c(-1,-1,1),lwd=c(1,1,2),cex=.8)

plot(0,0, ylab='Simulation', xlab='',ylim=c(0,10000),main='Number of times we reject H0')
count=0
for(i in 1:10000)
{
  if(T2.stat[i]>qchisq(1-alpha,p))
  {
    count=count+1
    points(0,count,col='orange',pch=19)
  }
}
abline(h=10000*alpha,lty=2,col='grey')
legend('topleft', c('(Empirical level)*10000', '(Nominal level)*10000'), col=c('orange','grey'), lty=c(-1,2),pch=c(19,-1),cex=.8)

# we compute the proportion of times on the total number
# of simulations in which I reject H0
# this is an empirical estimate of the level of the test
# (recall: H0 is true!)
sum(T2.stat > qchisq(1-alpha,p))/10000

# we test the distribution of T2 under H0:
# they should be (approximately) draws from a chi-squared
hist(T2.stat, prob=TRUE,col='grey85',ylim=c(0,.5),xlab='T2.stat',ylab='Density',main='Histogram of T2.stat')
xx <- seq(0,40,by=0.05)
lines(xx,dchisq(xx,p),type="l",lwd=2)
abline(v = qchisq(1-alpha,p), col='red', lwd=2)

quantili.H0 <- qchisq((1:10000 - 0.5)/10000,p)
qqplot(quantili.H0, T2.stat, asp=1)
abline(0,1)

graphics.off()

### 3) H0 FALSE
#-----------------------------------------
# We keep the same test but change the 
# distribution according to which we 
# generate the data 
# (=> H0 is false, we now aim to reject it)

mu <- c(2,1)

T2.stat <- NULL
for(i in 1:10000)
{
  x <- rmvnorm(n, mean=mu, sigma=sig)
  x.mean   <- colMeans(x)
  x.cov    <- cov(x)
  x.invcov <- solve(x.cov)
  x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
  T2.stat <- c(T2.stat, x.T2)
}
# The vector T2.stat collects 10000 independent realizations
# of the T2 statistics, computed under H0 (but now H0 is false!)

x11()
layout(cbind(c(1,1,1),c(2,3,4)), widths=c(2,1), heights=c(1,1))
plot(T2.stat, xlab='Simulation', ylab='T2 Statistics',xlim=c(0,10500))
abline(h = ((n-1)*p/(n-p))*qf(1-alpha,p,n-p), col='red', lwd=2)  # criterio di rifiuto
points(which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p)),T2.stat[which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))], col='orange')
points(rep(10500,10000),T2.stat, xlab='Simulazioni', ylab='Statistica T2',pch=19)
points(rep(10500, length(which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p)))),T2.stat[which(T2.stat>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))], pch=19, col='orange')
title(main='T2 statistics (Exact test, H0 false)')
legend('topleft',c('Accept H0', 'Reject H0', 'cfr.fisher'), col=c('black','orange','red'),pch=c(19,19,-1),lty=c(-1,-1,1),lwd=c(1,1,2),cex=.8)

plot(0,0, ylab='Simulation', xlab='',ylim=c(0,10000),main='Number of times we reject H0')
count=0
for(i in 1:10000)
{
  if(T2.stat[i]>((n-1)*p/(n-p))*qf(1-alpha,p,n-p))
  {
    count=count+1
    points(0,count,col='orange',pch=19)
  }
}

# we compute the proportion of times on the total number
# of simulations in which I reject H0
# this is an empirical estimate of the POWER of the test
# for mu = c(2,1) 
# (recall: H0 is false!)
sum(T2.stat > ((n-1)*p/(n-p))*qf(1-alpha,p,n-p))/10000

# we test the distribution of the T2 stat.:
# now H0 is false, we expect that they are NOT draws
# from a Fisher
hist(T2.stat, prob=TRUE,col='grey85',ylim=c(0,.5),xlab='T2.stat',ylab='Density',main='Histogram of T2.stat')
xx <- seq(0,40,by=0.05)
lines(xx,df(xx*(n-p)/((n-1)*p),p,n-p),type="l",lwd=2)
abline(v = ((n-1)*p/(n-p))*qf(1-alpha,p,n-p), col='red', lwd=2)

quantiles.H0 <- ((n-1)*p/(n-p))*qf((1:10000 - 0.5)/10000,p,n-p)
qqplot(quantiles.H0, T2.stat, asp=1,main='QQ-plot')
abline(0,1)

graphics.off()
# Remark: the power depends on the position of the theoretical
#         mean with respect to the mean under H0
