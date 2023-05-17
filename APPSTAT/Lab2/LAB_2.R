###--------------------###
###        LAB 2       ###
###--------------------###

### TOPICS:
### Probability density functions, Cumulative distribution functions, Quantiles
### Random number generation
### QQplots

#_______________________________________________________________________________
### Probability density functions, Cumulative distribution functions, Quantiles

### Example: Gaussian distribution
# Probability density function (pdf)
dnorm(0)                  # density function at 0 for a distribution N(0,1)
dnorm(0, mean=1, sd=2)    # density function at 0 for a distribution N(1,4)

# Cumulative distribution function (cdf)
pnorm(0)       # P(Z<0), with Z ~ N(0,1)
pnorm(0, 1, 2) # P(X<0), with X ~ N(1,4)

# Quantiles (inverse of cdf)
qnorm(0.95)        #  = z s.t.  P(Z<z)=0.95, with Z ~ N(0,1)
qnorm(0.95, 1, 2)  #  = z s.t.  P(Z<z)=0.95, with Z ~ N(1,4)

# Random generation
set.seed(8321)
rnorm(10)       # X_1,..,X_10 ~ N(0,1) i.i.d.
rnorm(10, 1, 2) # X_1,..,X_10 ~ N(1,4) i.i.d.

#########################################################
##### Commands to generate random numbers, obtain pdfs 
##### cdf and its inverse for the most popular models:
#
# Commands rnorm(),  dnorm(),  pnorm(),  qnorm(), 
#          rexp(),   dexp(),   pexp(),   qexp(),           
#          runif(),  dunif(),  punif(),  qunif(),      
#          rbinom(), dbinom(), pbinom(), qbinom(),      
#          rpois(),  dpois(),  ppois(),  qpois(),      
#          rgamma(), dgamma(), pgamma(), qgamma(),      
#########################################################

### Other examples of distributions
quartz() # x11() on PC
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, byrow=T))

s <- seq(-2, 2, by=0.01)

plot(s, dunif(s, 0, 1), main='Uniform Unif(0,1)', type='l', ylim=c(0, 1))
plot(s, dexp(s, 1), main='Exponential Exp(1)', type='l', ylim=c(0, 1))
plot(s, dnorm(s, 0, 1), main='Gaussian N(0,1)', type='l', ylim=c(0, 1))

plot(s, punif(s, 0, 1), main='Uniform Unif(0,1)', type='l', ylim=c(0, 1))
plot(s, pexp(s, 1), main='Exponential Exp(1)', type='l', ylim=c(0, 1))
plot(s, pnorm(s, 0, 1), main='Gaussian N(0,1)', type='l', ylim=c(0, 1))

w <- seq(0, 1, by=0.01)

plot(w, qunif(w, 0, 1), main='Uniform Unif(0,1)', type='l')
plot(w, qexp(w, 1),    main='Exponential Exp(1)', type='l')
plot(w, qnorm(w, 0, 1), main='Gaussian N(0,1)', type='l')

#_______________________________________________________________________________
### Generation of random numbers

set.seed(8321)
x <- runif(n=1000, min=0, max=1)
y <- rexp(n=1000, rate=1)
z <- rnorm(n=1000, mean=0, sd=1)

quartz()
par(mfrow=c(2, 3))
plot(x, main='Uniform Unif(0,1)')
plot(y, main='Exponential Exp(1)')
plot(z, main='Gaussian N(0,1)')

hist(x, main='', col='grey', xlab='x', prob=T)
lines(seq(-0.2, 1.2, length=100), dunif(seq(-0.2, 1.2, length=100)), col='blue', lty=2, lwd=2)
box()
hist(y, main='', col='grey', xlab='x', prob=T, ylim=c(0,1))
lines(seq(-1, 9, length=100), dexp(seq(-1, 9, length=100)), col='blue', lty=2, lwd=2)
box()
hist(z, main='', col='grey', xlab='x', prob=T, ylim=c(0,.45))
lines(seq(-4, 4, length=100), dnorm(seq(-4, 4, length=100)), col='blue', lty=2, lwd=2)
box()

graphics.off()

#_______________________________________________________________________________
### QQplot
# The QQplot (command qqplot()) can be used to verify if a sample comes
# from a given distribution

# Remark: the QQplot can be plotted for any distribution - not necessarily
# Gaussian 

quartz()
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, byrow=T))

plot(x, main='Uniform Unif(0,1)')
plot(y, main='Exponential Exp(1)')
plot(z, main='Gaussian N(0,1)')

hist(x, main='', col='grey', xlab='x', prob=T)
lines(seq(-0.2, 1.2, length=100), dunif(seq(-0.2, 1.2, length=100)), col='blue',
      lty=2, lwd=2)
box()
hist(y, main='', col='grey', xlab='x', prob=T, ylim=c(0,1))
lines(seq(-1, 9, length=100), dexp(seq(-1, 9, length=100)), col='blue', lty=2, lwd=2)
box()
hist(z, main='', col='grey', xlab='x', prob=T, ylim=c(0,.45))
lines(seq(-4, 4, length=100), dnorm(seq(-4, 4, length=100)), col='blue', lty=2, lwd=2)
box()

qqplot(qunif((1:1000/1000-0.5/1000)), x, col='red', xlab='Theoretical quantile',
       ylab='Sample Quantile', asp=1)
abline(0, 1, col='blue')
qqplot(qexp((1:1000/1000-0.5/1000)), y, col='red', xlab='Theoretical quantile',
       ylab='Sample Quantile', asp=1)
abline(0, 1, col='blue')
qqplot(qnorm((1:1000/1000-0.5/1000)), z, col='red', xlab='Theoretical quantile',
       ylab='Sample Quantile', asp=1)
abline(0, 1, col='blue')

# The QQplot is mostly used to qualitatively verify if a sample
# comes from a Gaussian distribution

quartz()
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, byrow=T))

x <- runif(n=1000, min=0, max=1)
y <- rexp(n=1000, rate=1)
z <- rnorm(n=1000, mean=0, sd=1)

qqplot(qnorm((1:1000/1000-1/2000)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, main='Unif(0,1)')
qqplot(qnorm((1:1000/1000-1/2000)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, main='Exp(1)')
qqplot(qnorm((1:1000/1000-1/2000)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, main='N(0,1)')

x <- runif(n=100, min=0, max=1)
y <- rexp(n=100, rate=1)
z <- rnorm(n=100, mean=0, sd=1)

qqplot(qnorm((1:100/100-1/200)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)
qqplot(qnorm((1:100/100-1/200)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)
qqplot(qnorm((1:100/100-1/200)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)

x <- runif(n=10, min=0, max=1)
y <- rexp(n=10, rate=1)
z <- rnorm(n=10, mean=0, sd=1)

qqplot(qnorm((1:10/10-1/20)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)
qqplot(qnorm((1:10/10-1/20)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)
qqplot(qnorm((1:10/10-1/20)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1)

# These plots can be obtained by using the command qqnorm() that automatically
# computes the theoretical quantiles of the Gaussian distribution.
# Further, R creates automatically the line with the command qqline(),
# that plots the straight line through the first and third quartile

# What happens when the sample comes from a non-standard Gaussian? 

quartz()
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, byrow=T))

x <- rnorm(n=1000, mean=0, sd=1)
y <- rnorm(n=1000, mean=2, sd=1)
z <- rnorm(n=1000, mean=0, sd=2)

qqplot(qnorm((1:1000/1000-1/2000)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2, main='N(0,1)')
abline(0, 1)
qqline(x, col='red')
qqplot(qnorm((1:1000/1000-1/2000)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2, main='N(2,1)')
abline(0,1)
qqline(y, col='red')
qqplot(qnorm((1:1000/1000-1/2000)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2, main='N(0,4)')
abline(0,1)
qqline(z, col='red')

x <- rnorm(n=100, mean=0, sd=1)
y <- rnorm(n=100, mean=2, sd=1)
z <- rnorm(n=100, mean=0, sd=2)

qqplot(qnorm((1:100/100-1/200)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0, 1)
qqline(x, col='red')
qqplot(qnorm((1:100/100-1/200)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0,1)
qqline(y, col='red')
qqplot(qnorm((1:100/100-1/200)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0,1)
qqline(z, col='red')

x <- rnorm(n=10, mean=0, sd=1)
y <- rnorm(n=10, mean=2, sd=1)
z <- rnorm(n=10, mean=0, sd=2)

qqplot(qnorm((1:10/10-1/20)), x, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0,1)
qqline(x, col='red')
qqplot(qnorm((1:10/10-1/20)), y, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0,1)
qqline(y, col='red')
qqplot(qnorm((1:10/10-1/20)), z, col='red', xlab='Theoretical quantile N(0,1)',
       ylab='Sample quantile', asp=1, ylim=c(-5,5)*2)
abline(0, 1)
qqline(z, col='red')

# If the data are Gaussian, the slope of the qqline is an estimate of 
# the standard deviation, the intercept is an estimate of the mean
# [this comes from the observation that: if x is the quantile of order
#  alpha of N(mu,sigma^2), than z=(x-mu)/sigma is the quantile of order
#  alpha of Z~N(0,1), i.e., x=mu+sigma*z].

### Test of Gaussianity: Shapiro-Wilks test
# H0: X ~ N     vs    H1=H0^c
shapiro.test(x)
shapiro.test(y)
shapiro.test(z)
