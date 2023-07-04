load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2018/2018-06-28/Precolombian.txt") # c=3
X.values <- X[c(1, 2)]
X.classes <- factor(X[[3]])
x11()
plot(X.values, col = ifelse(X.classes == levels(as.factor(X.classes))[1], "red", ifelse(X.classes == levels(as.factor(X.classes))[2], "blue", "green2")), main = "True classification")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue", "green2"), bty = "n")

# a) Assumptions:  gaussianity (qda/lda) and same variance of the groups (lda only)
mcshapiro.test(X.values[which(X.classes == levels(X.classes)[1]), ])$p
mcshapiro.test(X.values[which(X.classes == levels(X.classes)[2]), ])$p
mcshapiro.test(X.values[which(X.classes == levels(X.classes)[2]), ])$p
bartlett.test(X.values[which(X.classes == levels(X.classes)[1]), ], X.values[which(X.classes == levels(X.classes)[2]), ])
bartlett.test(X.values[which(X.classes == levels(X.classes)[2]), ], X.values[which(X.classes == levels(X.classes)[3]), ])
cov(X.values[which(X.classes == levels(X.classes)[1]), ])
cov(X.values[which(X.classes == levels(X.classes)[2]), ])
cov(X.values[which(X.classes == levels(X.classes)[3]), ])

# b) LDA
library(MASS) # for lda/qda
lda.mod <- lda(as.matrix(X.values), X.classes) # remember to add prior if needed
lda.mod
classification <- predict(lda.mod)$class
# Plot:
x11()
plot(X.values, col = ifelse(classification == levels(as.factor(X.classes))[1], "red", ifelse(classification == levels(as.factor(X.classes))[2], "blue", "green2")), main = "lda")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue", "green2"), bty = "n")
# Separation lines:
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = 200)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)
z.q <- predict(lda.mod, xy)$post
z1.q <- z.q[, 1] - pmax(z.q[, 2], z.q[, 3])
z2.q <- z.q[, 2] - pmax(z.q[, 3], z.q[, 1])
z3.q <- z.q[, 3] - pmax(z.q[, 2], z.q[, 1])
contour(x, y, matrix(z1.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z2.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z3.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)

# c) QDA
library(MASS) # for lda/qda
qda.mod <- qda(as.matrix(X.values), X.classes) # remember to add prior if needed
qda.mod
classification <- predict(qda.mod)$class
# Plot:
x11()
plot(X.values, col = ifelse(classification == levels(as.factor(X.classes))[1], "red", ifelse(classification == levels(as.factor(X.classes))[2], "blue", "green2")), main = "qda")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue", "green2"), bty = "n")
# Separation lines:
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = 200)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)
z.q <- predict(qda.mod, xy)$post
z1.q <- z.q[, 1] - pmax(z.q[, 2], z.q[, 3])
z2.q <- z.q[, 2] - pmax(z.q[, 3], z.q[, 1])
z3.q <- z.q[, 3] - pmax(z.q[, 2], z.q[, 1])
contour(x, y, matrix(z1.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z2.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z3.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)

# d) APER
misclass <- sum(classification != X.classes)
conf.mat <- table(true.class = X.classes, classifier = classification)
conf.mat
prior <- qda.mod$prior
APER <- (conf.mat[1, 2] + conf.mat[1, 3]) / (conf.mat[1, 2] + conf.mat[1, 1] + conf.mat[1, 3]) * prior[1] +
  (conf.mat[2, 1] + conf.mat[2, 3]) / (conf.mat[2, 1] + conf.mat[2, 2] + conf.mat[2, 3]) * prior[2] +
  (conf.mat[3, 1] + conf.mat[3, 2]) / (conf.mat[3, 1] + conf.mat[3, 2] + conf.mat[3, 3]) * prior[3]
APER

# e) AER through leave-one-out cross-validation
errors_CV <- 0
n <- dim(X)[1]
for (i in 1:n) {
  modCV.i <- qda(X.values[-i, ], X.classes[-i], prior = prior)
  errors_CV <- errors_CV + as.numeric(predict(modCV.i, X.values[i, ])$class != X.classes[i])
}
errors_CV
AERCV <- sum(errors_CV) / n
AERCV

# f)Probability of a new sample being classified as class i
i <- 1
Prob <- conf.mat[i, i] / (sum(conf.mat[i, ])) * prior[i] + t(conf.mat[-i, i] / (table(X.classes)[-i])) %*% prior[-i]
Prob

# g) Predict the class of a new unit
x0 <- data.frame(x1 = 50, x2 = 3.5)
names(x0) <- names(X.values)
predict(lda.mod, x0)

# h) k-nearest-neighbors
library(class) # for knn
k <- 1
x11()
plot(X.values, col = ifelse(X.classes == levels(as.factor(X.classes))[1], "red", ifelse(X.classes == levels(as.factor(X.classes))[2], "blue", "green2")), main = "True classification")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue", "green2"), bty = "n")
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = 200)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), length = 200)
xy <- expand.grid(x, y)
names(xy) <- names(X.values)
knn.mod <- knn(train = X.values, test = xy, cl = X.classes, k = k)
z <- as.numeric(knn.mod)
contour(x, y, matrix(z, 200), levels = c(1.5, 2.5, 3.5), drawlabels = F, add = T)
