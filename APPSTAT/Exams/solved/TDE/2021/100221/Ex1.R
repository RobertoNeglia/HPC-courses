### Ex1

dataset = read.table("shopping.txt", header = T)

n <- dim(dataset)[1]
p <- dim(dataset)[2]


# test of Gaussianity
load("mcshapiro.test.RData")
mcshapiro.test(dataset)$pvalue

x.mean   <- colMeans(dataset)
x.cov    <- cov(dataset)
x.invcov <- solve(x.cov)
# T statistics is the Mahalanobis distance between the sample mean and mu0
alpha = 0.05

# 3a) Verify if the test statistics belongs to the acceptance/rejection region
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)

# Center:
x.mean

# Directions of the principal axes:
eigen(x.cov/n)$vectors

# Length of the semi-axes of the ellipse:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(x.cov/n)$values) 



# b)
T2 <- cbind(inf = x.mean - sqrt(cfr.fisher*diag(x.cov)/n),
            center = x.mean, 
            sup = x.mean + sqrt(cfr.fisher*diag(x.cov)/n))
T2


a = c(0,1,1)

T2_sum = cbind(inf = t(a)%*%x.mean - sqrt(cfr.fisher*t(a)%*%(x.cov)%*%a/n),
               center = t(a)%*%x.mean, 
               sup = t(a)%*%x.mean + sqrt(cfr.fisher*t(a)%*%(x.cov)%*%a/n))

T2_sum


# c) 
a = c(-0.2,1 ,1)

mu0 = 0


ax.T2       <- sqrt(n) * (t(a)%*%x.mean-mu0)/ (sqrt(t(a)%*%(x.cov)%*%a))  

cfr.t = qt(1-alpha, n-1)

ax.T2 < cfr.t
# reject 
P = 1 - pt(ax.T2 , n-1)
P


ax.T2 < -cfr.t
# accept H0

P = pt(ax.T2 , n-1)
P

# pvalue = 0.99 -> accept 
























