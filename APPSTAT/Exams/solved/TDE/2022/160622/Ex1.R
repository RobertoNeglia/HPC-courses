### Ex1


library(car)
# read data from a table

data_1 <- read.table('discomaniac.txt', header = T)

data_2 <- read.table('lipsticks.txt', header = T)

D = data.frame( price =data_1$price-data_2$price, media_cond=data_1$media.condition - data_2$media.condition)

# Test the Gaussian assumption (on D!)
load("mcshapiro.test.RData")  # 0.744
mcshapiro.test(D)

n <- dim(D)[1]  
p <- dim(D)[2]  

D.mean   <- sapply(D,mean)   # sample mean
D.cov    <- cov(D)    # sample covariance
D.invcov <- solve(D.cov)

alpha   <- .05
delta.0 <- c(0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)  # test statistic
D.T2

cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher

D.T2 < cfr.fisher   # Q: accettiamo?

# we compute the p-value
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P



# Ellipsoidal confidence region with confidence level 95%
x11()
plot(D, asp=1, pch=1, main='Dataset of the Differences')
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2)

points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='grey35')


## Bonf
k <- 4  # 2
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.price <- c( D.mean[1]-cfr.t*sqrt(D.cov[1,1]/n) , D.mean[1], D.mean[1]+cfr.t*sqrt(D.cov[1,1]/n) )
IC.BF.cond  <- c( D.mean[2]-cfr.t*sqrt(D.cov[2,2]/n) , D.mean[2], D.mean[2]+cfr.t*sqrt(D.cov[2,2]/n) )

Bf <- rbind(IC.BF.price, IC.BF.cond)
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf

Bf_var <- cbind(inf=diag(D.cov)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
                center=diag(D.cov),
                sup=diag(D.cov)*(n-1) / qchisq(alpha/(2*k), n-1))
Bf_var
















