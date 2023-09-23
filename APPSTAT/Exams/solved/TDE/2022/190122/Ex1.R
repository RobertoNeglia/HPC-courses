### Ex1

pop1 <- read.table("acoruna.txt", header = T)
pop2 <- read.table("pontevedra.txt", header = T)



n1 <- dim(pop1)[1]
n2 <- dim(pop2)[1]
p <- dim(pop1)[2]
g <- 2 # number of populations

# we compute the sample mean, covariance matrices and the matrix Spooled
# we are assuming that the true covariance matrix is the same

pop1.mean <- sapply(pop1, mean)
pop2.mean <- sapply(pop2, mean)
pop1.cov <- cov(pop1)
pop2.cov <- cov(pop2)
Sp <- ((n1 - 1) * pop1.cov + (n2 - 1) * pop2.cov) / (n1 + n2 - 2)
# we compare the matrices
list(S1 = pop1.cov, S2 = pop2.cov, Spooled = Sp)

# gaussianity
load("mcshapiro.test.RData")
mcshapiro.test(pop1)$pvalue
mcshapiro.test(pop2)$pvalue


# to compare visually
x11()
par(mfrow = c(1, g))
image(pop1.cov, col = heat.colors(100), main = "Cov. S1", asp = 1, axes = FALSE, breaks = quantile(rbind(pop1.cov, pop2.cov), (0:100) / 100, na.rm = TRUE))
image(pop2.cov, col = heat.colors(100), main = "Cov. S2", asp = 1, axes = FALSE, breaks = quantile(rbind(pop1.cov, pop2.cov), (0:100) / 100, na.rm = TRUE))


alpha <- .01
delta.0 <- c(0, 0)
Spinv <- solve(Sp)

To <- n1 * n2 / (n1 + n2) * (pop1.mean - pop2.mean - delta.0) %*% Spinv %*% (pop1.mean - pop2.mean - delta.0)

cfr.fisher <- (p * (n1 + n2 - 2) / (n1 + n2 - 1 - p)) * qf(1 - alpha, p, n1 + n2 - 1 - p)
To < cfr.fisher # TQ : accept?

P <- 1 - pf(To / (p * (n1 + n2 - 2) / (n1 + n2 - 1 - p)), p, n1 + n2 - 1 - p)
P


# b)
alpha <- 0.01
k <- p
IC <- cbind(
    pop1.mean - pop2.mean - sqrt(diag(Sp) * (1 / n1 + 1 / n2)) * qt(1 - alpha / (k * 2), n1 + n2 - 2),
    pop1.mean - pop2.mean,
    pop1.mean - pop2.mean + sqrt(diag(Sp) * (1 / n1 + 1 / n2)) * qt(1 - alpha / (k * 2), n1 + n2 - 2)
)
IC

# c)

pop1 <- 0.5 * pop1[, 1] + 0.5 * pop1[, 2]
pop2 <- 0.5 * pop2[, 1] + 0.5 * pop2[, 2]

n1 <- 30
n2 <- 30
p <- 1
g <- 2 # number of populations

pop1.mean <- mean(pop1)
pop2.mean <- mean(pop2)
pop1.cov <- var(pop1)
pop2.cov <- var(pop2)
Sp <- ((n1 - 1) * pop1.cov + (n2 - 1) * pop2.cov) / (n1 + n2 - 2)
# we compare the matrices


alpha <- .01
delta.0 <- 0

a <- c(..)
mu0 <- ..
ax.T2 <- sqrt(n) * ((t(a) %*% x.mean) - mu0) / sqrt(t(a) %*% x.cov %*% a)

cfr.t <- qt(1 - alpha / 2, n - 1) # caso bilaterale
cfr.t <- qt(1 - alpha, n - 1) # caso unilaterale

abs(ax.T2) < cfr.t # Q: accettiamo?   <- caso bilaterale: H0: = mu0
ax.T2 > -cfr.t # Q: accettiamo?   <- caso bilaterale: H0: > mu0
ax.T2 < cfr.t # Q: accettiamo?   <- caso bilaterale: H0: < mu0

P <- 2 * (1 - pt(ax.T2, cft.t)) # test bilaterale
P <- pt(ax.T2, n - 1) # test unilaterale H0: > mu0
P <- 1 - pt(ax.T2, n - 1) # test unilaterale H0: < mu0

# H1:  mu.1 - mu.2 > 0
