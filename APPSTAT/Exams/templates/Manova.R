load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2017/2017-07-18/horsecolic.txt")

# a) Manova
X.values <- X[c(1, 2, 3, 4)]
X.groups <- X[[5]]
man <- manova(as.matrix(X.values) ~ X.groups)
summary(man)
summary.aov(man) # to see non simultaneous anova tests
# Check assumption (normality, same variance between groups)
mcshapiro.test(X.values[which(X.groups == "Yes"), ])$p
mcshapiro.test(X.values[which(X.groups == "No"), ])$p
cov(X.values[which(X.groups == "Yes"), ])
cov(X.values[which(X.groups == "No"), ])
bartlett.test(X.values, X.groups) # ok if pvalue big

# b) Bonferroni conf. int. for the difference of means
group1 <- X.values[which(X.groups == "Yes"), ]
group2 <- X.values[which(X.groups == "No"), ]
g <- length(unique(X.groups))
p <- dim(X.values)[2]
n1 <- dim(group1)[1]
n2 <- dim(group2)[1]
n <- n1 + n2
k <- p * g * (g - 1) / 2
man <- manova(as.matrix(X.values) ~ X.groups)
SSres <- summary.manova(man)$SS$Residuals
Spooled <- SSres / (n - g)
mean.g1 <- sapply(group1, mean)
mean.g2 <- sapply(group2, mean)
alpha <- 0.01
qT <- qt(1 - alpha / 2 / k, n - g)
conf.int.B <- cbind(
    inf = mean.g1 - mean.g2 - qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2)),
    mean.dif = mean.g1 - mean.g2,
    sup = mean.g1 - mean.g2 + qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2))
)
conf.int.B
