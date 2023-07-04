load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2017/2017-07-03/geisha.txt")
#  X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-02-09/streaming.txt", header=TRUE)
if (dim(X)[2] == 2) {
  x11()
  plot(X, main = "Dataset")
}

# a) Hierarchical clustering
D <- dist(X) # distance matrix
clust <- hclust(D, method = "single") # linkages: single,complete, average
x11()
plot(clust, hang = -0.1, xlab = "", labels = F, cex = 0.6, sub = "") # dendrogram
groups <- cutree(clust, 3)
# centers:
if (length(table(groups)) == 2) {
  centers <- rbind(
    sapply(X[which(groups == 1), ], mean),
    sapply(X[which(groups == 2), ], mean)
  )
} else {
  centers <- rbind(
    sapply(X[which(groups == 1), ], mean),
    sapply(X[which(groups == 2), ], mean),
    sapply(X[which(groups == 3), ], mean)
  )
}
centers
if (dim(X)[2] == 2) {
  x11()
  plot(X, col = groups + 1, main = "Clustering") # red=group 1, green=group 2
  points(centers, pch = 19)
}
# sizes:
sizes <- table(groups)
sizes
# cophenetic coefficient:
cor(cophenetic(clust), D)

# b) K-means
result.k <- kmeans(X, centers = 2) # centers: fixed number of clusters
groups.k <- result.k$cluster # labels of clusters
centers.k <- result.k$centers # centers of the clusters
# result.k$totss        # tot. sum of squares
# result.k$withinss     # sum of squares within clusters
# result.k$tot.withinss # sum(sum of squares within cluster)
# result.k$betweenss    # sum of squares between clusters
sizes.k <- result.k$size # dimension of the clusters
if (dim(X)[2] == 2) {
  x11()
  plot(X, col = groups.k + 1) # red=group 1, green=group 2
  points(centers.k, pch = 19)
}
if (dim(X)[2] == 3) {
  open3d()
  plot3d(X, size = 3, col = groups.k + 1, aspect = F)
  points3d(centers.k, size = 10)
}

# c) Choosing k in K-means (elbow method)
b <- NULL
w <- NULL
max.k <- 5
for (k in 1:max.k) {
  result.k <- kmeans(X, k)
  w <- c(w, sum(result.k$wit))
  b <- c(b, result.k$bet)
}
x11()
matplot(1:max.k, w / (w + b), pch = "", xlab = "clusters", ylab = "within/tot", main = "Choice of k", ylim = c(0, 1))
lines(1:max.k, w / (w + b), type = "b", lwd = 2)

# d) Bonferroni conf. int. for the difference of mean of groups
# Assumptions (gaussianity and same variance)
load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
mcshapiro.test(X[which(groups == 1), ])$p
mcshapiro.test(X[which(groups == 2), ])$p
bartlett.test(X[which(groups == 1), ], X[which(groups == 2), ])

man <- manova(as.matrix(X) ~ factor(groups))
summary.manova(man, test = "Wilks")
n1 <- table(groups)[1]
n2 <- table(groups)[2]
n <- n1 + n2
p <- dim(X)[2]
g <- 2
k <- p * g * (g - 1) / 2
alpha <- 0.05
qT <- qt(1 - alpha / 2 / k, n - g)
m1 <- sapply(X[which(groups == 1), ], mean)
m2 <- sapply(X[which(groups == 2), ], mean)
S <- summary.manova(man)$SS$Residuals
Spooled <- S / (n - g)
conf.int.diff12 <- cbind(
  inf = m1 - m2 - qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2)),
  dif_mean = m1 - m2,
  sup = m1 - m2 + qT * sqrt(diag(Spooled) * (1 / n1 + 1 / n2))
)
conf.int.diff12

# e) Bonferroni conf. int. for the mean of one of the groups
G <- X[which(groups == 1), ]
load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
mcshapiro.test(G)$p # Assumption: gaussianity
n <- dim(G)[1]
p <- dim(G)[2]
sample.mean <- sapply(G, mean)
S <- cov(G)
invS <- solve(S)
alpha <- 0.05
qT <- qt(1 - alpha / 2 / p, n - 1)
conf.int.B <- cbind(
  inf = sample.mean - qT * sqrt(diag(S) / n),
  mean = sample.mean,
  sup = sample.mean + qT * sqrt(diag(S) / n)
)
conf.int.B
