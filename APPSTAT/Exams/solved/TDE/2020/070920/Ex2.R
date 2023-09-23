### Ex2


pop1 = read.table("candle.txt", header = T)
pop2 = read.table("sunshine.txt", header = T)


n1 <- dim(pop1)[1] 
n2 <- dim(pop2)[1] 
p  <- dim(pop1)[2] 
g = 2

# we compute the sample mean, covariance matrices and the matrix Spooled
# we are assuming that the true covariance matrix is the same 

pop1.mean <- sapply(pop1,mean)
pop2.mean <- sapply(pop2,mean)
pop1.cov  <-  cov(pop1)
pop2.cov  <-  cov(pop2)
Sp      <- ((n1-1)*pop1.cov + (n2-1)*pop2.cov)/(n1+n2-2)

load("mcshapiro.test.RData") 
mcshapiro.test(pop1)$pvalue
mcshapiro.test(pop2)$pvalue


# to compare visually
x11()
par(mfrow=c(1,g))
image(pop1.cov, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(pop1.cov,pop2.cov), (0:100)/100, na.rm=TRUE))
image(pop2.cov, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(pop1.cov,pop2.cov), (0:100)/100, na.rm=TRUE))

# b)
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .05
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

To <- n1*n2/(n1+n2) * (pop1.mean-pop2.mean-delta.0) %*% Spinv %*% (pop1.mean-pop2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
To < cfr.fisher # TQ : accept?

P <- 1 - pf(To/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P

# c)
alpha <- 0.05
k = p
IC <- cbind(pop1.mean-pop2.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k*2), n1+n2-2),
            pop1.mean-pop2.mean,
            pop1.mean-pop2.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k*2), n1+n2-2))
IC


# d)

# 
diff1 = data.frame(pop1[,1]-pop1[,2])
diff2 = data.frame(pop2[,1]-pop2[,2])

pop1.mean <- sapply(diff1,mean)
pop2.mean <- sapply(diff2,mean)
pop1.cov  <-  cov(diff1)
pop2.cov  <-  cov(diff2)
Sp      <- ((n1-1)*pop1.cov + (n2-1)*pop2.cov)/(n1+n2-2)

alpha   <- .05
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

To <- n1*n2/(n1+n2) * (pop1.mean-pop2.mean-delta.0) %*% Spinv %*% (pop1.mean-pop2.mean-delta.0)

# H0 : mu1 -mu2 <= 0
# H1 : mu1 -mu2 > 0
# reject if T0 is too high
cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
To < cfr.fisher # TQ : accept?























