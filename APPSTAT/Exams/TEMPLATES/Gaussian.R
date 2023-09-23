X <- read.table("~/GitHub/Applied-Statistics-Exam/Exams of previous years/2018/2018-09-13/IAMG.txt") # example file p=3
# X <- read.table("~/GitHub/Applied-Statistics-Exam/Exams of previous years/2017/2017-07-18/castle.txt") #example file p=2

# a) Conf. region (ellipse) for mean of multivariate gaussian (95%)
n <- dim(X)[1]
p <- dim(X)[2]
sample.mean <- sapply(X, mean)
S <- cov(X)
invS <- solve(S)
alpha <- 0.05
qfish <- (n - 1) * p / (n - p) * qf(1 - alpha, p, n - p)
# The analytical expression is (x-sample.mean)%*%invS%*%(x-sample.mean) <= qfish/n
# notice that we divide qfish by n one extra time, this comes form the variance which is S/n
center <- sample.mean
radius <- sqrt(qfish / n)
axes <- eigen(S)$vectors
lenghts <- sqrt(eigen(S)$values) * radius
if (p == 2) {
  library(car)
  x11()
  plot(X, asp = 1, pch = 1, main = "Dataset")
  car::ellipse(
    center = sample.mean,
    shape = S,
    radius = sqrt(qfish * (1 / n)),
    lwd = 2,
    col = "blue"
  ) # 95% of means, replace 1/n with 1/n+1 to get 95% of values
  points(45.733, 7.333, pch = 19, col = "red")
}

# b) T2 confidence intervals for the mean (fisher)
n <- dim(X)[1]
p <- dim(X)[2]
sample.mean <- sapply(X, mean)
S <- cov(X)
invS <- solve(S)
alpha <- 0.05
qfish <- (n - 1) * p / (n - p) * qf(1 - alpha, p, n - p)
conf.int.T2 <- cbind(
  inf = sample.mean - sqrt(qfish) * sqrt(diag(S) / n),
  mean = sample.mean,
  sup = sample.mean + sqrt(qfish) * sqrt(diag(S) / n)
)

# c) Bonferroni confidence intervals for the mean (t-student)
n <- dim(X)[1]
p <- dim(X)[2]
sample.mean <- sapply(X, mean)
S <- cov(X)
invS <- solve(S)
alpha <- 0.05
qT <- qt(1 - alpha / 2 / p, n - 1)
conf.int.B <- cbind(
  inf = sample.mean - qT * sqrt(diag(S) / n),
  mean = sample.mean,
  sup = sample.mean + qT * sqrt(diag(S) / n)
)

# d) Hypothesis test | e.g (p=3) H0: 0.1*X$1 - X$3 = 0, bilateral, alpha = 0.05
n <- dim(X)[1]
alpha <- 0.05
R <- cbind(0.1, 0, 1)
sample.mean <- sapply(X, mean)
mu0 <- 0
S <- cov(X)
qT <- qt(1 - alpha / 2, n - 1)
T.stat <- sqrt(n * (R %*% sample.mean - mu0) %*% solve(R %*% S %*% t(R)) %*% (R %*% sample.mean - mu0))
T.stat > qT # if true = reject H0
pvalue <- 1 - pt(T.stat, n - 1)

# e) Confidence intervals for the variance (p=3)
n <- dim(X)[1]
alpha <- 0.05
SS <- diag(cov(X))
chi_up <- qchisq(1 - alpha / 2, n - 1)
chi_down <- qchisq(alpha / 2, n - 1)
CI2 <- cbind(
  inf = SS * (n - 1) / chi_up,
  center = SS,
  sup = SS * (n - 1) / chi_down
)
CI2

# f) Test on the mean of a gaussian
n <- dim(D)[1]
p <- dim(D)[2]
mu0 <- rep(0, p)
sample.mean <- sapply(D, mean)
S <- cov(D)
invS <- solve(S)
alpha <- 0.05
qfish <- (n - 1) * p / (n - p) * qf(1 - alpha, p, n - p)
T_0 <- n * (sample.mean - mu0) %*% invS %*% (sample.mean - mu0)
reject <- T_0 > qfish
p.value <- 1 - pf(T_0 * (n - p) / (p * (n - 1)), p, n - p)
print(paste("Reject H0:", reject, ", p-value:", p.value))
