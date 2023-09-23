## Ex1


dataset <- read.table('pollution.txt', header = T)
head(dataset)
dim(dataset)

n <- dim(dataset)[1]
p <- dim(dataset)[2]


# test of Gaussianity
load("mcshapiro.test.RData")
mcshapiro.test(dataset)$pvalue

mu0   <- c(50,50)
alpha <- 0.05


# 2) Compute the test statistics
x.mean   <- colMeans(dataset)
x.cov    <- cov(dataset)
x.invcov <- solve(x.cov)

x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# T statistics is the Mahalanobis distance between the sample mean and mu0

# 3a) Verify if the test statistics belongs to the acceptance/rejection region
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
x.T2 < cfr.fisher # Q: accettiamo?

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


library(car)



x11()
plot(dataset, asp = 1)
ellipse(x.mean, x.cov/n, sqrt(cfr.fisher), col = 'red', lty = 2, lwd=2, center.cex=1)


T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2


















