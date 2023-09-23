library(MVN)

## Test for the mean of two indepenedent Gaussian populations
setwd("~/HPC/APPSTAT/Exams/2022-01-19/E1")
city1 <- read.table("acoruna.txt", header = TRUE)
city2 <- read.table("pontevedra.txt", header = TRUE)

dim(city1)
dim(city2)

head(city1)
head(city2)

n1 <- dim(city1)[1] # n1 = 30
n2 <- dim(city2)[1] # n2 = 30
p <- dim(city1)[2] # p = 2

# Multivariate normality
mvn(city1)$multivariateNormality$`p value`
mvn(city2)$multivariateNormality$`p value`

# Means and covariances
city1.mean <- colMeans(city1)
city1.mean
city2.mean <- colMeans(city2)
city2.mean

city1.cov <- cov(city1)
city1.cov
city2.cov <- cov(city2)
city2.cov



# Spooled covariance
Spooled <- ((n1 - 1) * city1.cov + (n2 - 1) * city2.cov) / (n1 + n2 - 2)
Spooled

# Statistical test of 99% confidence to verify if the means differ
# Test H0: mu1 = mu2 vs H1: mu1 != mu2
# i.e.
# Test H0: mu1 - mu2 = c(0, 0) vs H1: mu1 - mu2 != c(0, 0)
alpha <- 0.01
diff <- c(0, 0)
Spooled.inv <- solve(Spooled)

T2 <- n1 * n2 / (n1 + n2) * (city1.mean - city2.mean - diff) %*% Spooled.inv %*% (city1.mean - city2.mean - diff)
T2

F.test <- (p * (n1 + n2 - 2) / (n1 + n2 - p - 1)) * qf(1 - alpha, p, n1 + n2 - p - 1)
T2 < F.test
# FALSE: there is statistical evidence reject the null hypothesis that the means are equal

# Compute the p-value
p.value <- 1 - pf(T2 / (p * (n1 + n2 - 2) / (n1 + n2 - p - 1)), p, n1 + n2 - p - 1)
p.value
# indeed the p-value is very small, so we reject the null hypothesis

# question b)
# Bonferroni confidence intervals for the difference of the means.
# I have to compute 2 comparisons: one for the first component and one
# for the second component
k <- 2

qT <- qt(1 - alpha / (2 * k), n1 + n2 - 2)
qT

BC.I <- cbind(
    lower = (city1.mean - city2.mean) - qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2)),
    mean = city1.mean - city2.mean,
    upper = (city1.mean - city2.mean) + qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2))
)
BC.I

svg("bonferroni.svg", width = 5, height = 5)
plot(c(1, 2), range(BC.I), pch = "")
lines(x = c(1, 1), y = c(BC.I[1, 1], BC.I[1, 3]), col = "red")
lines(x = c(2, 2), y = c(BC.I[2, 1], BC.I[2, 3]), col = "blue")
abline(h = 0, col = "gray")
points(city1.mean - city2.mean, col = c("red", "blue"), pch = 16)
dev.off()

# question c)
avg.eval <- cbind(city1 = rowMeans(city1), city2 = rowMeans(city2))
head(avg.eval)
# check for normality (it should be normal because it is the average of two
# normal variables)
mvn(avg.eval)$multivariateNormality$`p value`

avg.mean <- colMeans(avg.eval)
avg.mean
avg.mean[1] - avg.mean[2]

avg.cov <- cov(avg.eval)
avg.cov

Sp <- sqrt(avg.cov[1, 1]^2 * (n1 - 1) + avg.cov[2, 2]^2 * (n2 - 1)) /
    sqrt(n1 + n2 - 2)
Sp
# Test H0: mu1 = mu2 vs H1: mu1 != mu2
# i.e.
# Test H0: mu1 - mu2 = 0 vs H1: mu1 - mu2 != 0
alpha <- 0.01

qT <- qt(1 - alpha / 2, n1 + n2 - 2)
qT

t <- (avg.mean[1] - avg.mean[2]) / (Sp * sqrt(1 / n1 + 1 / n2))
t

t > qT

p.value <- 2 * pt(-abs(t), n1 + n2 - 2)
p.value

# p-value is smaller than any significance level, so we reject the null hypothesis

CI <- c(
    lower = avg.mean[1] - avg.mean[2] - qT * Sp * sqrt(1 / n1 + 1 / n2),
    upper = avg.mean[1] - avg.mean[2] + qT * Sp * sqrt(1 / n1 + 1 / n2)
)
CI

svg("CI-average-diff.svg", width = 5, height = 5)
plot(c(0.9, 1.1), c(-0.01, 0.01), pch = "", xlab = "", ylab = "")
lines(x = c(1, 1), y = c(CI[1], CI[2]), col = "red")
abline(h = 0, col = "gray")
dev.off()
