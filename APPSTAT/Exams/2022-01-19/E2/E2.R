setwd("~/HPC/APPSTAT/Exams/2022-01-19/E2")
load("mcshapiro.test.RData")
library(MASS)
library(heplots)


fishes <- read.table("fish.txt", header = TRUE)
dim(fishes)
head(fishes)

attach(fishes)

coords <- cbind(x, y)
head(coords)

labels <- factor(abundance)

where.low <- which(abundance == "L")
where.high <- which(abundance == "H")

dim(coords[where.low, ])
dim(coords[where.high, ])
levels(labels)

mean.H <- apply(coords[where.high, ], 2, mean)
mean.L <- apply(coords[where.low, ], 2, mean)
means <- rbind(mean.H, mean.L)
means
svg("clusters.svg", width = 6, height = 6)
plot(coords[where.low, ], col = "red", pch = 19, xlab = "x", ylab = "y")
points(coords[where.high, ], col = "blue", pch = 19)
points(means, col = c("blue", "red"), pch = 3, cex = 2, lwd = 5)
legend("topright", legend = levels(labels), col = c("blue", "red"), pch = 19)
dev.off()

# question a+b)

# verify assumptions: normality and homoscedasticity
# normality (univariate within groups)
mvn(coords[where.low, ])$multivariateNormality$`p value`
mvn(coords[where.high, ])$multivariateNormality$`p value`
# I can assume the normality only for the high abundance group

# homoscedasticity (multivariate between groups)
summary(boxM(coords, labels))

# Check with MANOVA to see if I can use (x,y) to discriminate between
# the two groups
fit <- manova(coords ~ labels)
summary.manova(fit, test = "Wilks")
# very small p-value, so I can use (x,y) to discriminate between the two
# groups

# I can still try to use LDA and QDA and see what happens
fish.lda <- lda(coords, labels, data = fishes)
fish.lda

# plot the classification regions
plot(coords[where.low, ], col = "red", pch = 19, xlab = "x", ylab = "y")
points(coords[where.high, ], col = "blue", pch = 19)
points(means, col = c("blue", "red"), pch = 3, cex = 2, lwd = 5)
legend("topright", legend = levels(labels), col = c("blue", "red"), pch = 19)

x <- seq(min(coords[, 1]), max(coords[, 1]), length = 200)
y <- seq(min(coords[, 2]), max(coords[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)
z <- predict(fish.lda, newdata = xy)$post
z1 <- z[, 1] - z[, 2]
z2 <- z[, 2] - z[, 1]

contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)

# plot the LDA
plot(fish.lda, dimen = 1, type = "both", xlab = "x", ylab = "y")

lda.fish <- predict(fish.lda, newdata = coords)

table(class.true = labels, class.lda = lda.fish$class)
errors <- (lda.fish$class != labels)

APER <- sum(errors) / length(labels)
APER

(17 + 18) / 250

# Use cross validation:
fish.lda.cv <- lda(coords, labels, data = fishes, CV = TRUE)

table(class.true = labels, class.assignedCV = fish.lda.cv$class)
errors.cv <- (fish.lda.cv$class != labels)

APER.cv <- sum(errors.cv) / length(labels)
APER.cv

(17 + 20) / 250
# not that different from the APER without cross validation

# QDA
fish.qda <- qda(coords, labels, data = fishes)
fish.qda

plot(coords[where.low, ], col = "red", pch = 19, xlab = "x", ylab = "y")
points(coords[where.high, ], col = "blue", pch = 19)
points(means, col = c("blue", "red"), pch = 3, cex = 2, lwd = 5)
legend("topright", legend = levels(labels), col = c("blue", "red"), pch = 19)

x <- seq(min(coords[, 1]), max(coords[, 1]), length = 200)
y <- seq(min(coords[, 2]), max(coords[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)
z <- predict(fish.qda, newdata = xy)$post
z1 <- z[, 1] - z[, 2]
z2 <- z[, 2] - z[, 1]

contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)

qda.fish <- predict(fish.qda, newdata = coords)
table(class.true = labels, class.qda = qda.fish$class)

errors.qda <- (qda.fish$class != labels)
APER.qda <- sum(errors.qda) / length(labels)
APER.qda

(12 + 12) / 250

# Use cross validation:
fish.qda.cv <- qda(coords, labels, data = fishes, CV = TRUE)
table(class.true = labels, class.assignedCV = fish.qda.cv$class)

errors.qda.cv <- (fish.qda.cv$class != labels)
APER.qda.cv <- sum(errors.qda.cv) / length(labels)
APER.qda.cv

(12 + 13) / 250

# question c)
library(class)

set.seed(19)
n <- dim(coords)[1]
ks <- 10:30
AERCV <- NULL
for (k in ks) {
    errors_CV <- 0
    for (i in 1:n) {
        modCV.i <- knn(train = coords[-i, ], test = coords[i, ], labels[-i], k = k)
        errors_CV <- errors_CV + as.numeric(modCV.i != labels[i])
    }
    AERCV <- c(AERCV, sum(errors_CV) / n)
}

plot(ks, AERCV)
ks[which(min(AERCV) == AERCV)]

plot(coords[where.low, ], col = "red", pch = 19, xlab = "x", ylab = "y")
points(coords[where.high, ], col = "blue", pch = 19)
points(means, col = c("blue", "red"), pch = 3, cex = 2, lwd = 5)
legend("topright", legend = levels(labels), col = c("blue", "red"), pch = 19)

x <- seq(min(coords[, 1]), max(coords[, 1]), length = 200)
y <- seq(min(coords[, 2]), max(coords[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)
mod.plot <- knn(train = coords, test = xy, labels, k = 13)
z <- as.numeric(mod.plot)

contour(x, y, matrix(z, 200), levels = c(1.5, 2.5), drawlabels = FALSE, lty = 2, add = TRUE)


# question d)
new.location <- data.frame(x = 10.8, y = 39.4)

plot(coords[where.low, ], col = "red", pch = 19, xlab = "x", ylab = "y")
points(coords[where.high, ], col = "blue", pch = 19)
points(means, col = c("blue", "red"), pch = 3, cex = 2, lwd = 5)
points(new.location, col = "darkgreen", pch = 19)
legend("topright", legend = c(levels(labels), "new location"), col = c("blue", "red", "green"), pch = 19)

x <- seq(min(coords[, 1]), max(coords[, 1]), length = 200)
y <- seq(min(coords[, 2]), max(coords[, 2]), length = 200)
xy <- expand.grid(x = x, y = y)

# lda contour plot
z.lda <- predict(fish.lda, newdata = xy)$post
z1.lda <- z.lda[, 1] - z.lda[, 2]
z2.lda <- z.lda[, 2] - z.lda[, 1]

contour(x, y, matrix(z1.lda, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)
contour(x, y, matrix(z2.lda, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)

# qda contour plot
z.qda <- predict(fish.qda, newdata = xy)$post
z1.qda <- z.qda[, 1] - z.qda[, 2]
z2.qda <- z.qda[, 2] - z.qda[, 1]

contour(x, y, matrix(z1.qda, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)
contour(x, y, matrix(z2.qda, 200), levels = 0, drawlabels = FALSE, lty = 2, add = TRUE)

# knn controur plot
mod.plot <- knn(train = coords, test = xy, labels, k = 13)
z <- as.numeric(mod.plot)

contour(x, y, matrix(z, 200), levels = c(1.5, 2.5), drawlabels = FALSE, lty = 2, add = TRUE)


predict(fish.lda, newdata = new.location)
# result is H

predict(fish.qda, newdata = new.location)
# result is L



knn(train = coords, test = new.location, labels, k = 13)
# result is H


detach(fishes)
