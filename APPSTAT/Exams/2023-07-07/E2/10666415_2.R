setwd("~/shared-folder/HPC/APPSTAT/Exams/2023-07-07/E2")
rm(list = ls())

consumption <- read.table("consumption.txt", header = TRUE)
dim(consumption)
head(consumption)

cons.mean <- colMeans(consumption)
cons.mean
library(MVN)
mvn(consumption)$multivariateNormality$`p value`

svg("matplot.svg", width = 6, height = 6)
matplot(t(consumption), type = "l", lty = 1)
dev.off()

n <- dim(consumption)[1]
p <- dim(consumption)[2]

M <- sapply(consumption, mean)
S <- cov(consumption)

S

C <- matrix(c(
    -1, 1, 0, 0, 0, 0,
    0, 0, -1, 1, 0, 0,
    0, 0, 0, 0, -1, 1
), nrow = 3, ncol = 6, byrow = TRUE)
C

alpha <- 0.05
delta.0 <- c(0, 0, 0)

Md <- C %*% M
Md

Sd <- C %*% S %*% t(C)
Sd

T2 <- n * t(Md - delta.0) %*% solve(Sd) %*% (Md - delta.0)

qF <- ((p - 1) * (n - 1) / (n - (p - 1))) * qf(1 - alpha, p - 1, n - p + 1)

T2 < qF
# Output is FALSE -> we reject H0 at 5%

P <- 1 - pf(T2 * (n - (p - 1)) / ((p - 1) * (n - 1)), p - 1, n - p + 1)
P
