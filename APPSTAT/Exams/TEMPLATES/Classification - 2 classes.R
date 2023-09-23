load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2017/2017-07-18/horsecolic.txt")
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-06-16/Exercise 2/musicCountry.txt",header = TRUE)
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-01-19/fish.txt")

X.values <- X[c(1, 2)]
X.classes <- factor(X[[3]])
x11()
plot(X.values, col = ifelse(X.classes == levels(as.factor(X.classes))[1], "red", "blue"), main = "True classification")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue"), bty = "n")

# a) Assumptions:  gaussianity (qda/lda) and same variance of the groups (lda only)
load("~/GitHub/Applied-Statistics-Exam/mcshapiro.test.RData")
mcshapiro.test(X.values[which(X.classes == levels(X.classes)[1]), ])$p
mcshapiro.test(X.values[which(X.classes == levels(X.classes)[2]), ])$p
bartlett.test(X.values[which(X.classes == levels(X.classes)[1]), ], X.values[which(X.classes == levels(X.classes)[2]), ])
cov(X.values[which(X.classes == levels(X.classes)[1]), ])
cov(X.values[which(X.classes == levels(X.classes)[2]), ])

# b) LDA
library(MASS) # for lda/qda
lda.mod <- lda(as.matrix(X.values), X.classes) # remember to add prior if needed
lda.mod
classification <- predict(lda.mod)$class
# Plot:
x11()
plot(X.values, col = ifelse(classification == levels(as.factor(X.classes))[1], "red", "blue"), main = "LDA")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue"), bty = "n")
points(X.values[which(classification == levels(as.factor(X.classes))[1] & X.classes == levels(as.factor(X.classes))[2]), ], col = "red", pch = 19) # misclassified 1-> true 2
points(X.values[which(classification == levels(as.factor(X.classes))[2] & X.classes == levels(as.factor(X.classes))[1]), ], col = "blue", pch = 19) # misclassified 2 -> true 1
# Separation lines:
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = 200)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), length = 200)
xy <- expand.grid(Altezza = x, Peso = y)
z.q <- predict(lda.mod, xy)$post
z1.q <- z.q[, 1] - z.q[, 2]
z2.q <- z.q[, 2] - z.q[, 1]
contour(x, y, matrix(z1.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z2.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)

# c) QDA
library(MASS) # for lda/qda
qda.mod <- qda(as.matrix(X.values), X.classes) # remember to add prior if needed
qda.mod
classification <- predict(qda.mod)$class
# Plot:
x11()
plot(X.values, col = ifelse(classification == levels(as.factor(X.classes))[1], "red", "blue"), main = "QDA")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue"), bty = "n")
points(X.values[which(classification == levels(as.factor(X.classes))[1] & X.classes == levels(as.factor(X.classes))[2]), ], col = "red", pch = 19) # misclassified 1-> true 2
points(X.values[which(classification == levels(as.factor(X.classes))[2] & X.classes == levels(as.factor(X.classes))[1]), ], col = "blue", pch = 19) # misclassified 2 -> true 1
# Separation lines:
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = 200)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), length = 200)
xy <- expand.grid(Altezza = x, Peso = y)
names(xy) <- names(X.values)
z.q <- predict(qda.mod, xy)$post
z1.q <- z.q[, 1] - z.q[, 2]
z2.q <- z.q[, 2] - z.q[, 1]
contour(x, y, matrix(z1.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)
contour(x, y, matrix(z2.q, 200), levels = 0, drawlabels = F, add = T, lty = 2)

# d) APER
# APER <- 0
# for(g in 1:G)
#   APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]
misclass <- sum(classification != X.classes)
conf.mat <- table(true.class = X.classes, classifier = classification)
conf.mat
prior <- qda.mod$prior
APER <- conf.mat[1, 2] / (conf.mat[1, 2] + conf.mat[1, 1]) * prior[1] + conf.mat[2, 1] / (conf.mat[2, 1] + conf.mat[2, 2]) * prior[2]
APER


# e) AER through leave-one-out cross-validation
errors_CV <- 0
n <- dim(X)[1]
for (i in 1:n) {
  modCV.i <- qda(X.values[-i, ], X.classes[-i]) # remembere to use the right model, add the prior if needed
  errors_CV <- errors_CV + as.numeric(predict(modCV.i, X.values[i, ])$class != X.classes[i])
}
errors_CV
AERCV <- sum(errors_CV) / n
AERCV

# f)Probability of a new sample being classified as class i
i <- 1
Prob <- conf.mat[i, i] / (sum(conf.mat[i, ])) * prior[i] + t(conf.mat[-i, i] / (table(X.classes)[-i])) %*% prior[-i]
Prob

# g) SVM
library(e1071) # for svm
dat <- data.frame(X.values, y = as.factor(X.classes)) # to use svm you must build a data.frame where the class must be called y and must be factor,
# no other feature can be called y!!
svm.mod <- svm(y ~ ., data = dat, kernel = "radial", cost = 10, scale = FALSE) # y~. means use all the features, kernel= linear/radial...
summary(svm.mod)
classification <- predict(svm.mod, X.values)
# Plot 1:
x11()
par(mfrow = c(1, 2))
plot(svm.mod, dat, col = c("salmon", "light blue"), pch = 19) # notice that it switches the coordinates, the "x" are the support vectors
# Plot 2 & 3:
n.g <- 100
xgrid <- expand.grid(
  x.1 = seq(from = range(dat[[1]])[1], to = range(dat[[1]])[2], length = n.g),
  x.2 = seq(from = range(dat[[2]])[1], to = range(dat[[2]])[2], length = n.g)
)
names(xgrid) <- names(X.values)
ygrid <- predict(svm.mod, xgrid)
x11()
plot(xgrid, col = c("red", "blue")[as.numeric(ygrid)], pch = 20, cex = .2, main = "SVM")
points(X.values, col = c("red", "blue")[as.numeric(dat$y)], pch = 19)
points(X.values[svm.mod$index, ], pch = 5, cex = 2) # support vectors
x11()
plot(X.values, col = c("red", "blue")[as.numeric(dat$y)], pch = 19, main = "SVM")
contour(seq(from = range(dat[[1]])[1], to = range(dat[[1]])[2], length = n.g),
  seq(from = range(dat[[2]])[1], to = range(dat[[2]])[2], length = n.g),
  matrix(as.numeric(ygrid), n.g, n.g),
  level = 1.5, add = TRUE,
  drawlabels = F
)

# h) Tune the cost parameter in SVM (10-fold cross-validation)
set.seed(1)
tune.out <- tune(svm, y ~ .,
  data = dat, kernel = "linear",
  ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10, 100))
)
summary(tune.out)
bestmod <- tune.out$best.model # Extract the best model from the result of tune
summary(bestmod)
x11()
plot(bestmod, dat, col = c("salmon", "light blue"), pch = 19)

# i) Predict the class of a new unit
x0 <- data.frame(x = 10.8, y = 39.4)
names(x0) <- names(X.values)
predict(svm.mod, x0)

# j) k-nearest-neighbors
library(class) # for knn
k <- 4
knn.mod <- knn(train = X.values, test = X.values, cl = X.classes, k = k)
# Plot:
x11()
plot(X.values, col = ifelse(X.classes == levels(as.factor(X.classes))[1], "red", "blue"), main = "True classification")
legend("topright", levels(as.factor(X.classes)), fill = c("red", "blue"), bty = "n")
grid.pts <- 100
x <- seq(min(X.values[, 1]), max(X.values[, 1]), length = grid.pts)
y <- seq(min(X.values[, 2]), max(X.values[, 2]), lengthgrid.pts)
xy <- expand.grid(x, y)
names(xy) <- names(X.values)
knn.mod.plot <- knn(train = X.values, test = xy, cl = X.classes, k = k)
z <- as.numeric(knn.mod.plot)
contour(x, y, matrix(z, grid.pts), levels = c(1.5, 2.5), drawlabels = F, add = T)

# k) Tune k in k-n-n through leave-one-out cross-validation
set.seed(19)
n <- dim(X)[1]
ks <- 10:30
AERCV <- NULL
library(class) # for knn
for (k in ks) {
  errors_CV <- 0
  for (i in 1:n) {
    modCV.i <- knn(train = X.values[-i, ], test = X.values[i, ], X.classes[-i], k = k)
    errors_CV <- errors_CV + as.numeric(modCV.i != X.classes[i])
  }
  AERCV <- c(AERCV, sum(errors_CV) / n)
}
x11()
plot(AERCV)
best.k <- which.min(AERCV)
min(AERCV)

# l) Predict with knn
x0 <- data.frame(x = 10.8, y = 39.4)
names(x0) <- names(X.values)
knn(train = X.values, test = x0, cl = X.classes, k = k)
