load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2017/2017-07-03/bento.txt")
D <- X[c(1, 2, 3, 4)] - X[c(5, 6, 7, 8)]
# Or if the datasets are in two different files:
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-01-19/acoruna.txt")
# Y <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-01-19/pontevedra.txt")
# D <- X-Y
boxplot(D)

# a) Check assumptions: gaussianity of difference dataset
mcshapiro.test(D)$p

# b) Test on the mean of a gaussian with H_0: mu = c(0,0,0,0)
n <- dim(D)[1]
p <- dim(D)[2]
mu0 <- rep(0, p)
sample.mean <- sapply(D, mean)
S <- cov(D)
invS <- solve(S)
alpha <- 0.01
qfish <- (n - 1) * p / (n - p) * qf(1 - alpha, p, n - p)
T_0 <- n * (sample.mean - mu0) %*% invS %*% (sample.mean - mu0)
reject <- T_0 > qfish
p.value <- 1 - pf(T_0 * (n - p) / (p * (n - 1)), p, n - p)
print(paste("Reject H0:", reject, ", p-value:", p.value))

# c) T2 intervals to check if there is INCREASE in mean (fisher)
n <- dim(D)[1]
p <- dim(D)[2]
mu0 <- rep(0, p)
sample.mean <- sapply(D, mean)
S <- cov(D)
alpha <- 0.05
qfish <- (n - 1) * p / (n - p) * qf(1 - alpha, p, n - p)
rejection.region.T2 <- cbind("-inf", mu0 + sqrt(qfish) * sqrt(diag(S) / n))
print(
    cbind(
        rejection.region.T2,
        sample.mean,
        ifelse(
            sample.mean > as.numeric(rejection.region[, 2]),
            "Increase",
            "No statistical increase"
        )
    )
)

# d) Bonferroni intervals to check if there is INCREASE in mean (t-student)
n <- dim(D)[1]
p <- dim(D)[2]
k <- p
mu0 <- rep(0, p)
sample.mean <- sapply(D, mean)
S <- cov(D)
alpha <- 0.05
# we don't divide alpha by 2 since we wnat to check the increase (unilateral)
qT <- qt(1 - alpha / k, n - 1)
rejection.region.b <- cbind("-inf", mu0 + qT * sqrt(diag(S) / n))
print(
    cbind(
        rejection.region.b,
        sample.mean,
        ifelse(
            sample.mean > as.numeric(rejection.region[, 2]),
            "Increase",
            "No statistical increase"
        )
    )
)
