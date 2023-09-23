### Ex1

dataset <- read.table('chicca.txt', header = T)

n <- dim(dataset)[1]
p <- dim(dataset)[2]



load("mcshapiro.test.RData")
mcshapiro.test(dataset)$pvalue

mu0   <- c(0, 90)
alpha <- 0.01


# 2) Compute the test statistics
x.mean   <- colMeans(dataset)
x.cov    <- cov(dataset)
x.invcov <- solve(x.cov)

x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# T statistics is the Mahalanobis distance between the sample mean and mu0

# 3a) Verify if the test statistics belongs to the acceptance/rejection region
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
x.T2 < cfr.fisher # Q: accettiamo?
# reject

# 3b) Compute the p-value
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P


# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values) 
# Warning: Conf Reg => x.cov/n 

# Question: how to plot the confidence region? 
# (I don't know how to plot in R^p with p>3!)

# SE SIAMO NEL CASO P=2 POSSO PLOTTARE L'ELLISSE della regione di rifiuto e CR

# Region of rejection (centered in mu0)
# NB: la regione di rifiuto è il complementare dell'ellisse centrata in mu0! 
library(car)
x11()
plot(dataset, asp = 1)
ellipse(mu0, shape=x.cov/n, sqrt(cfr.fisher), col = 'blue', lty = 2, center.pch = 16)
# We add a red point in correspondence of the sample mean
points(x.mean[1], x.mean[2], pch = 16, col ='red', cex = 1.5)


# c)
k <- p
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt(diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt(diag(x.cov)/n))
Bf


# d) 

DATA = dataset
DATA$delay[which(DATA$delay < 0)] = 0

x.mean2   <- colMeans(DATA)
x.cov2    <- cov(DATA)
x.invcov2 <- solve(x.cov)

mu0 = 90
a = c(1,1)
alpha = 0.1
ax.T2       <- sqrt(n) *((t(a)%*%x.mean2)-mu0)/sqrt(t(a)%*%x.cov2%*%a)

cfr.t = qt(1 - alpha/2, n-1)

abs(ax.T2) < cfr.t
# reject

P = 2*(1-pt(ax.T2, cfr.t))
P




























