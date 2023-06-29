### --------------------###
### LAB 7 (14/04/2023) ###
### --------------------###

### TOPICS:
### Linear and Quadratic Discriminant Analysis
### Support Vector Machines

options(rgl.printRglwidget = TRUE)

# _______________________________________________________________________________
###
### Example 1 (LDA with 2 classes, univariate)
### ----------------------------------
setwd("~/HPC/APPSTAT/Lab7")
cyto <- read.table("cytokines.txt", header = TRUE)
cyto

attach(cyto)

# Group A: favorable clinical outcome (the treatment has effect)
A <- which(group == "A")
# Group B: unfavorable clinical outcome (the treatment has no effect)
B <- which(group == "B")

plot(cyto[, 1], cyto[, 2],
   pch = 19, col = c(rep("blue", 8), rep("red", 5)),
   xlab = "Inf-g", ylab = "IL-5"
)

dev.off()

### Idea: we aim to find a "rule" to classify patients as Group A or Group B
### given the measurements of Inf-g (Interferon gamma) and IL-5 (Interleukin 5).
### For this example, we discard the IL-5 variable (which seems to have poor
### "discriminative power") and we consider the univariate case with Inf-g only

### LDA (univariate)
### ------------------
# Assumptions:
# 1) if L = i, X.i ~ N(mu.i, sigma.i^2), i = A,B
# 2) sigma.A = sigma.B
# 3) c(A|B) = c(B|A) (equal misclassification costs)

# verify assumptions 1) and 2):
# 1) normality (univariate) within the groups
shapiro.test(cyto[A, 1])
shapiro.test(cyto[B, 1])
# -> OK

# 2) equal variance (univariate)
# var.test is used instead of bartlett.test when there are only 2 groups
var.test(cyto[A, 1], cyto[B, 1])
# -> OK

nA <- length(A)
nB <- length(B)
n <- nA + nB

# Prior probabilities (estimated from the data, no prior knowledge)
PA <- nA / n
PB <- nB / n

PA
PB

### Recall:
# the classification region is obtained by comparing pA*f.A and pB*f.B

MA <- mean(Infg[A])
MB <- mean(Infg[B])

SA <- var(Infg[A])
SB <- var(Infg[B])
S <- ((nA - 1) * SA + (nB - 1) * SB) / (nA + nB - 2) # pooled estimate

x <- seq(-10, 35, 0.5) # include the range of Infg

par(mfrow = c(2, 1))
plot(x, PA * dnorm(x, MA, sqrt(S)),
   type = "l", col = "blue",
   ylab = expression(paste("estimated ", p[i] * f[i], "(x)")), main = "LDA"
)

points(x, PB * dnorm(x, MB, sqrt(S)), type = "l", col = "red")
points(Infg[A], rep(0, length(A)), pch = 16, col = "blue")
points(Infg[B], rep(0, length(B)), pch = 16, col = "red")

legend(-10, 0.03,
   legend = c(
      expression(paste("P(A)", f[A], "(x)")),
      expression(paste("P(B)", f[B], "(x)"))
   ),
   col = c("blue", "red"), lty = 1, cex = 0.7
)

plot(x,
   PA * dnorm(x, MA, sqrt(S)) / (PA * dnorm(x, MA, sqrt(S)) + PB * dnorm(x, MB, sqrt(S))),
   type = "l", col = "blue", ylab = "estimated posterior"
)

points(x,
   PB * dnorm(x, MB, sqrt(S)) / (PA * dnorm(x, MA, sqrt(S)) + PB * dnorm(x, MB, sqrt(S))),
   type = "l", col = "red"
)

points(Infg[A], rep(0, length(A)), pch = 16, col = "blue")
points(Infg[B], rep(0, length(B)), pch = 16, col = "red")

legend(-10, 0.9,
   legend = c("P(A|X=x)", "P(B|X=x)"), col = c("blue", "red"), lty = 1,
   cex = 0.7
)

dev.off()
### end Recall

### LDA with R
library(MASS)
help(lda)

cyto.lda <- lda(group ~ Infg)
cyto.lda
# Note: if we don't specify the prior probabilities, they are estimated
# from the sample

# posterior probability and classification for x=0
x <- data.frame(Infg = 0)
x
# The command predict() returns a list containing (see the help of predict.lda):
# - the class associated with the highest posterior probability
predict(cyto.lda, x)$class
# - the posterior probabilities for the classes
predict(cyto.lda, x)$posterior
# - in lda: the coordinates of the canonical analysis of Fisher
#           (Fisher's discriminant scores)
predict(cyto.lda, x)$x


plot(Infg[A], rep(0, length(A)),
   pch = 16, col = "blue", ylim = c(0, 1),
   xlab = "x", ylab = "estimated posterior", main = "LDA", xlim = range(Infg)
)
points(Infg[B], rep(0, length(B)), pch = 16, col = "red")
abline(v = 0, col = "grey")
points(c(0, 0), c(predict(cyto.lda, data.frame(Infg = 0))$posterior),
   col = c("blue", "red"), pch = "*", cex = 2.5
)
abline(v = predict(cyto.lda, data.frame(Infg = 0))$x, col = "blue")


# posterior probability for a grid of x's
x <- data.frame(Infg = seq(-10, 35, 0.5))

head(predict(cyto.lda, x)$posterior)
# posterior probability for class A
cyto.LDA.A <- predict(cyto.lda, x)$posterior[, 1]
# posterior probability for class B
cyto.LDA.B <- predict(cyto.lda, x)$posterior[, 2]

predict(cyto.lda, x)$class

lines(
   x[, 1],
   cyto.LDA.A,
   type = "l",
   col = "blue",
   xlab = "x",
   ylab = "estimated posterior",
   main = "LDA"
)
lines(x[, 1], cyto.LDA.B, type = "l", col = "red")
abline(h = 0.5)
legend(
   -10, 0.9,
   legend = c("P(A|X=x)", "P(B|X=x)"),
   fill = c("blue", "red"),
   cex = 0.7
)

dev.off()


# set prior probabilities
# Pay attention to the order in which the levels appear in factor(group)!
factor(group)
cyto.lda.1 <- lda(group ~ Infg, prior = c(0.95, 0.05))
cyto.lda.1

x <- data.frame(Infg = seq(-10, 35, 0.5))

# posterior probability for class A
cyto.LDA.A.1 <- predict(cyto.lda.1, x)$posterior[, 1]
# posterior probability for class B
cyto.LDA.B.1 <- predict(cyto.lda.1, x)$posterior[, 2]


plot(
   x[, 1],
   cyto.LDA.A.1,
   type = "l",
   col = "blue",
   xlab = "x",
   ylab = "estimated posterior",
   main = "LDA",
   ylim = c(0, 1)
)
points(x[, 1], cyto.LDA.B.1, type = "l", col = "red")
abline(h = 0.5)
legend(
   -10, 0.9,
   legend = c("P(A|X=x)", "P(B|X=x)"),
   fill = c("blue", "red"),
   cex = 0.7
)
points(Infg[A], rep(0, length(A)), pch = 16, col = "blue")
points(Infg[B], rep(0, length(B)), pch = 16, col = "red")

points(x[, 1], cyto.LDA.A, type = "l", col = "grey")
points(x[, 1], cyto.LDA.B, type = "l", col = "grey")



### k-nearest neighbor classifier
### -------------------------------

library(class)
help(knn)

cyto.knn <- knn(train = Infg, test = x, cl = group, k = 3, prob = TRUE)
summary(cyto.knn)
cyto.knn.class <- cyto.knn == "B"
cyto.knn.class

cyto.knn.B <- ifelse(cyto.knn.class == TRUE,
   attributes(cyto.knn)$prob,
   1 - attributes(cyto.knn)$prob
)


plot(
   x[, 1],
   cyto.LDA.B,
   type = "l",
   col = "red",
   lty = 2,
   xlab = "x",
   ylab = "estimated posterior"
)
points(x[, 1], cyto.knn.B, type = "l", col = "black", lty = 1)
abline(h = 0.5)
legend(
   -10, 0.75,
   legend = c("LDA", "knn"),
   lty = c(2, 1),
   col = c("red", "black")
)

# let's change k

par(mfrow = c(3, 4))
for (k in 1:12)
{
   cyto.knn <- knn(train = Infg, test = x, cl = group, k = k, prob = TRUE)
   cyto.knn.class <- (cyto.knn == "B") + 0
   cyto.knn.B <-
      ifelse(
         cyto.knn.class == 1,
         attributes(cyto.knn)$prob,
         1 - attributes(cyto.knn)$prob
      )
   plot(
      x[, 1],
      cyto.LDA.B,
      type = "l",
      col = "red",
      lty = 2,
      xlab = "x",
      ylab = "estimated posterior",
      main = k
   )
   points(x[, 1], cyto.knn.B, type = "l", col = "black", lty = 1, lwd = 2)
   abline(h = 0.5)
}

detach(cyto)
graphics.off()

# _____________________________________________________________________________
### Example 2 (3 classes, bivariate)
### -------------------------------------------

library(MASS)
library(MVN) # Can need separate installation of "gsl" software library
library(car)
library(mvtnorm)
library(mvnormtest)
library(heplots)

help(iris)
rm(list = ls())
### We consider only the first two variables, Sepal.Length and Sepal.Width
### (p=2, g=3); we aim to build a classifier based on the characteristic
### of the sepal that identifies the iris species.

attach(iris)
species.name <- factor(Species)

g <- length(levels(species.name))
g

i1 <- which(species.name == "setosa")
i2 <- which(species.name == "versicolor")
i3 <- which(species.name == "virginica")

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n <- n1 + n2 + n3

detach(iris)

iris2 <- iris[, 1:2]
iris2

# Jittering
set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd = 0.025))
iris2.jitt

# plot the data
par(mfrow = c(1, 2))
plot(iris2, main = "Iris Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 19)
points(iris2[i1, ], col = "red", pch = 19)
points(iris2[i2, ], col = "green", pch = 19)
points(iris2[i3, ], col = "blue", pch = 19)
legend("topright", legend = levels(species.name), fill = c("red", "green", "blue"))

plot(iris2.jitt, main = "Iris Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 19)
points(iris2.jitt[i1, ], col = "red", pch = 19)
points(iris2.jitt[i2, ], col = "green", pch = 19)
points(iris2.jitt[i3, ], col = "blue", pch = 19)
legend("topright", legend = levels(species.name), fill = c("red", "green", "blue"))

dev.off()

iris2 <- iris2.jitt

m <- colMeans(iris2)
m1 <- colMeans(iris2[i1, ])
m2 <- colMeans(iris2[i2, ])
m3 <- colMeans(iris2[i3, ])

mvn(iris2[i1, ])$multivariateNormality$`p value`
mvn(iris2[i2, ])$multivariateNormality$`p value`
mvn(iris2[i3, ])$multivariateNormality$`p value`

S1 <- cov(iris2[i1, ])
S2 <- cov(iris2[i2, ])
S3 <- cov(iris2[i3, ])
Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2 + (n3 - 1) * S3) / (n - g)

par(mfrow = c(1, 3))
image(
   S1,
   col = heat.colors(100),
   main = "Cov. S1",
   asp = 1,
   axes = FALSE,
   breaks = quantile(rbind(S1, S2, S3), (0:100) / 100, na.rm = TRUE)
)
image(
   S2,
   col = heat.colors(100),
   main = "Cov. S2",
   asp = 1,
   axes = FALSE,
   breaks = quantile(rbind(S1, S2, S3), (0:100) / 100, na.rm = TRUE)
)
image(
   S3,
   col = heat.colors(100),
   main = "Cov. S3",
   asp = 1,
   axes = FALSE,
   breaks = quantile(rbind(S1, S2, S3), (0:100) / 100, na.rm = TRUE)
)

summary(boxM(iris2, species.name))
# small pvalue but this is too restrictive for MANOVA anyway

# One-way MANOVA (See LAB 6)
# I use MANOVA to check if I can use the first two variables to discriminate
# the three species.
fit <- manova(as.matrix(iris2) ~ species.name)
summary.manova(fit, test = "Wilks")

# Linear Discriminant Analysis (LDA)
lda.iris <- lda(iris2, species.name)

# "coefficients of linear discriminants" and "proportion of trace":
#                     Fisher discriminant analysis.
# In particular:
# - coefficients of linear discriminants: versors of the canonical directions
#   [to be read column-wise]
# - proportion of trace: proportion of variance (between wrt to within)
#   explained by the corresponding canonical direction
# As we can see from the proportion of trace, the first canonical direction
# explains 96.28% of the variance, while the second one only 3.72%.

Lda.iris <- predict(lda.iris, iris2)
Lda.iris

# Estimate of AER (actual error rate):
# 1) APER (apparent error rate)
# 2) estimate of AER by cross-validation

# 1) Compute the APER
Lda.iris$class # assigned classes
species.name # true labels
table(class.true = species.name, class.assigned = Lda.iris$class)

errors <- (Lda.iris$class != species.name)
errors
sum(errors)
length(species.name)

APER <- sum(errors) / length(species.name)
APER

(1 + 14 + 15) / 150

# Remark: this is correct only if we estimate the prior with the empirical
#         frequencies! Otherwise:
# prior <- c(1/3,1/3,1/3)
# G <- 3
# misc <- table(class.true=species.name, class.assigned=Lda.iris$class)
# APER <- 0
# for(g in 1:G)
#   APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior[g]

# 2) Compute the estimate of the AER by leave-one-out cross-validation

### Recall: LOO CV

errors_CV <- 0
for (i in 1:150) {
   LdaCV.i <- lda(iris2[-i, ], species.name[-i], prior = c(50, 50, 50) / 150)
   errors_CV <- errors_CV + as.numeric(predict(LdaCV.i, iris2[i, ])$class != species.name[i])
}
errors_CV

AERCV <- sum(errors_CV) / length(species.name)
AERCV

####

# with R:
LdaCV.iris <- lda(iris2, species.name, CV = TRUE) # specify the argument CV

LdaCV.iris$class
species.name
table(class.true = species.name, class.assignedCV = LdaCV.iris$class)

errorsCV <- (LdaCV.iris$class != species.name)
errorsCV
sum(errorsCV)

AERCV <- sum(errorsCV) / length(species.name)
AERCV
# Remark: correct only if we estimate the priors through the sample frequencies!


# Plot the partition induced by LDA
plot(iris2, main = "Iris Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20)
points(iris2[i1, ], col = "red", pch = 20)
points(iris2[i2, ], col = "green", pch = 20)
points(iris2[i3, ], col = "blue", pch = 20)
legend("topright", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

points(lda.iris$means, pch = 4, col = c("red", "green", "blue"), lwd = 2, cex = 1.5)

x <- seq(min(iris[, 1]), max(iris[, 1]), length = 200)
y <- seq(min(iris[, 2]), max(iris[, 2]), length = 200)
xy <- expand.grid(Sepal.Length = x, Sepal.Width = y)

z <- predict(lda.iris, xy)$post # these are P_i*f_i(x,y)
z1 <- z[, 1] - pmax(z[, 2], z[, 3]) # P_1*f_1(x,y)-max{P_j*f_j(x,y)}
z2 <- z[, 2] - pmax(z[, 1], z[, 3]) # P_2*f_2(x,y)-max{P_j*f_j(x,y)}
z3 <- z[, 3] - pmax(z[, 1], z[, 2]) # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

# Plot the contour line of level (levels=0) of z1, z2, z3:
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z3, 200), levels = 0, drawlabels = FALSE, add = TRUE)

library(rgl)
library(mvtnorm)
open3d()
points3d(iris2[i1, 1], iris2[i1, 2], 0, col = "red", pch = 15)
points3d(iris2[i2, 1], iris2[i3, 2], 0, col = "green", pch = 15)
points3d(iris2[i3, 1], iris2[i2, 2], 0, col = "blue", pch = 15)
surface3d(x, y, dmvnorm(xy, m1, Sp) / 3, alpha = 0.4, color = "red")
surface3d(x, y, dmvnorm(xy, m2, Sp) / 3, alpha = 0.4, color = "green", add = TRUE)
surface3d(x, y, dmvnorm(xy, m3, Sp) / 3, alpha = 0.4, color = "blue", add = TRUE)
box3d()


### Quadratic Discriminant Analysis (QDA)
### ---------------------------------------

help(qda)

qda.iris <- qda(iris2, species.name)
qda.iris

Qda.iris <- predict(qda.iris, iris2)

# compute the APER
Qda.iris$class
species.name
table(class.true = species.name, class.assigned = Qda.iris$class)

errorsq <- (Qda.iris$class != species.name)
errorsq

APERq <- sum(errorsq) / length(species.name)
APERq

(15 + 13 + 1) / 150

# Remark: correct only if we estimate the priors through the sample frequencies!

# Compute the estimate of the AER by leave-one-out cross-validation
QdaCV.iris <- qda(iris2, species.name, CV = TRUE)
QdaCV.iris$class
species.name
table(class.true = species.name, class.assignedCV = QdaCV.iris$class)

errorsqCV <- (QdaCV.iris$class != species.name)
errorsqCV

AERqCV <- sum(errorsqCV) / length(species.name)
AERqCV
# Remark: correct only if we estimate the priors through the sample frequencies!

# Plot the partition induced by QDA

plot(iris2, main = "Iris Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20)
points(iris2[i1, ], col = "red", pch = 20)
points(iris2[i2, ], col = "green", pch = 20)
points(iris2[i3, ], col = "blue", pch = 20)
legend("topright", legend = levels(species.name), fill = c("red", "green", "blue"))

points(qda.iris$means, col = c("red", "green", "blue"), pch = 4, lwd = 2, cex = 1.5)

x <- seq(min(iris[, 1]), max(iris[, 1]), length = 200)
y <- seq(min(iris[, 2]), max(iris[, 2]), length = 200)
xy <- expand.grid(Sepal.Length = x, Sepal.Width = y)

z <- predict(qda.iris, xy)$post
z1 <- z[, 1] - pmax(z[, 2], z[, 3])
z2 <- z[, 2] - pmax(z[, 1], z[, 3])
z3 <- z[, 3] - pmax(z[, 1], z[, 2])

contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z3, 200), levels = 0, drawlabels = FALSE, add = TRUE)

open3d()
points3d(iris2[i1, 1], iris2[i1, 2], 0, col = "red", pch = 15)
points3d(iris2[i2, 1], iris2[i3, 2], 0, col = "green", pch = 15)
points3d(iris2[i3, 1], iris2[i2, 2], 0, col = "blue", pch = 15)
surface3d(x, y, dmvnorm(xy, m1, S1) / 3, alpha = 0.4, color = "red")
surface3d(x, y, dmvnorm(xy, m2, S2) / 3, alpha = 0.4, color = "green", add = TRUE)
surface3d(x, y, dmvnorm(xy, m3, S3) / 3, alpha = 0.4, color = "blue", add = TRUE)
box3d()


### knn-classifier
### ----------------
# Plot the partition induced by knn
library(class)
k <- 7

plot(iris2, main = "Iris.Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20)
points(iris2[i1, ], col = 2, pch = 20)
points(iris2[i3, ], col = 4, pch = 20)
points(iris2[i2, ], col = 3, pch = 20)
legend("topright", legend = levels(species.name), fill = c(2, 3, 4))

x <- seq(min(iris[, 1]), max(iris[, 1]), length = 200)
y <- seq(min(iris[, 2]), max(iris[, 2]), length = 200)
xy <- expand.grid(Sepal.Length = x, Sepal.Width = y)

iris.knn <- knn(train = iris2, test = xy, cl = iris$Species, k = k)

z <- as.numeric(iris.knn)

contour(x, y, matrix(z, 200), levels = c(1.5, 2.5), drawlabels = FALSE, add = TRUE)

graphics.off()

# _______________________________________________________________________________
##### Problem 2 of 9/09/2009
##### --------------------------
# The vending machines of the Exxon fuel contain optical detectors able
# to measure the size of the banknotes inserted. Knowing that 0.1% of the
# 10$ banknotes in circulation are counterfeit, Exxon would like to implement a
# software to identify false 10$ banknotes, as to minimize the economic losses.
# Assuming that:
#  - both the populations of real and false banknotes follow a normal
#    distribution (with different mean and covariance matrices);
#  - accepting a false banknote leads to an economic loss of 10$;
#  - rejecting a true banknote brings a economic loss quantifiable in 5 cents;
# satisfy the following requests of the Exxon:
# a) build an appropriate classifier, estimating the unknown parameters
#    starting from the two datasets moneytrue.txt and moneyfalse.txt, containing
#    data about 100 true banknotes and 100 counterfeit banknotes (in mm).
#    Qualitatively show the two classification regions in a graph;
# b) calculate the APER of the classifier and, based on the APER, estimate the
#    expected economic damage of the classifier;
# c) what is the estimated probability that the first 10$ banknote inserted in the
#    machine is rejected?
rm(list = ls())
library(MVN)
load("mcshapiro.test.RData")

true <- read.table("moneytrue.txt", header = TRUE)
false <- read.table("moneyfalse.txt", header = TRUE)

banknotes <- rbind(true, false)
vf <- factor(rep(c("true", "false"), each = 100), levels = c("true", "false"))

plot(banknotes[, 1:2], main = "Banknotes", xlab = "V1", ylab = "V2", pch = 20)
points(false, col = "red", pch = 20)
points(true, col = "blue", pch = 20)
legend("bottomleft", legend = levels(vf), fill = c("blue", "red"), cex = .7)

# question a)
mcshapiro.test(true)
mcshapiro.test(false)
mvn(true)$multivariateNormality
mvn(false)$multivariateNormality

# misclassification costs
# label as true a false banknote
c.vf <- 10
# label as false a true banknote
c.fv <- 0.05

# prior probabilities
pf <- 0.001
pt <- 1 - 0.001

# Prior modified to account for the misclassification costs
prior.c <-
   c(
      pt * c.fv / (c.vf * pf + c.fv * pt),
      pf * c.vf / (c.vf * pf + c.fv * pt)
   )
prior.c

# QDA
qda.m <- qda(banknotes, vf, prior = prior.c)
qda.m


plot(banknotes[, 1:2], main = "Banknotes", xlab = "V1", ylab = "V2", pch = 20)
points(false, col = "red", pch = 20)
points(true, col = "blue", pch = 20)
legend("bottomleft", legend = levels(vf), fill = c("blue", "red"), cex = .7)

points(qda.m$means, pch = 4, col = c("red", "blue"), lwd = 2, cex = 1.5)

x <- seq(min(banknotes[, 1]), max(banknotes[, 1]), length = 200)
y <- seq(min(banknotes[, 2]), max(banknotes[, 2]), length = 200)
xy <- expand.grid(V1 = x, V2 = y)

z <- predict(qda.m, xy)$post
z1 <- z[, 1] - z[, 2]
z2 <- z[, 2] - z[, 1]

contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, add = TRUE)

m1 <- colMeans(true)
m2 <- colMeans(false)

S1 <- cov(true)
S2 <- cov(false)

open3d()
points3d(true[, 1], true[, 2], 0, col = "blue", pch = 15)
points3d(false[, 1], false[, 2], 0, col = "red", pch = 15)
surface3d(x, y, matrix(dmvnorm(xy, m1, S1) * prior.c[1] / 100, 200), alpha = 0.6, color = "blue")
surface3d(x, y, matrix(dmvnorm(xy, m2, S2) * prior.c[2] / 100, 200), alpha = 0.6, color = "red", add = TRUE)
box3d()

# LDA
lda.m <- lda(banknotes, vf, prior = prior.c)
lda.m

Lda.m <- predict(lda.m)
table(class.true = vf, class.assigned = Lda.m$class)

plot(banknotes[, 1:2], main = "Banknotes", xlab = "V1", ylab = "V2", pch = 20)
points(false, col = "red", pch = 20)
points(true, col = "blue", pch = 20)
legend("bottomleft", legend = levels(vf), fill = c("blue", "red"), cex = .7)

points(lda.m$means, pch = 4, col = c("red", "blue"), lwd = 2, cex = 1.5)

z <- predict(lda.m, xy)$post
z1 <- z[, 1] - z[, 2]
z2 <- z[, 2] - z[, 1]

contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, add = TRUE)

# LDA can't discriminate between the two classes!

# question b)

# APER
Qda.m <- predict(qda.m)
table(class.true = vf, class.assigned = Qda.m$class)

APER <- (2 * pt + 80 * pf) / 100
APER

# Expected economic loss:
2 / 100 * pt * c.fv + 80 / 100 * pf * c.vf
# here 2 / 100 is P[rejected | true]
#      80 / 100 is P[accepted | false]

# question c)
# P[rejected] = P[rejected | true]P[true] + P[rejected | false]P[false]
2 / 100 * pt + 20 / 100 * pf

# LOO cross-validation
qda.m.loo <- qda(banknotes, vf, prior = prior.c, CV = TRUE)
qda.m.loo
table(class.true = vf, class.assigned = qda.m.loo$class)

APER.loo <- (2 * pt + 82 * pf) / 100
APER.loo

2 / 100 * pt * c.fv + 82 / 100 * pf * c.vf

2 / 100 * pt + 18 / 100 * pf

### -------------------------------------------------------------------------
### FISHER DISCRIMINANT ANALYSIS
### -----------------------------------
### Let's change viewpoint: we look for the directions that highlight
### the discrimination among groups
### -> we look for the canonical directions

# Remark. Assumptions: homogeneity of the covariance structure
# [we relax the normality assumption]

# Let's consider again the iris dataset
rm(list = ls())
attach(iris)

species.name <- factor(Species, labels = c("setosa", "versicolor", "virginica"))

g <- 3

i1 <- which(species.name == "setosa")
i2 <- which(species.name == "versicolor")
i3 <- which(species.name == "virginica")

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n <- n1 + n2 + n3

detach(iris)

iris2 <- iris[, 1:2]
head(iris2)

set.seed(1)
iris2 <- iris2 + cbind(rnorm(150, sd = 0.025)) # jittering

m <- colMeans(iris2)
m1 <- colMeans(iris2[i1, ])
m2 <- colMeans(iris2[i2, ])
m3 <- colMeans(iris2[i3, ])

S1 <- cov(iris2[i1, ])
S2 <- cov(iris2[i2, ])
S3 <- cov(iris2[i3, ])
Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2 + (n3 - 1) * S3) / (n - g)

# covariance between groups (estimate)
B <- 1 / g * (cbind(m1 - m) %*% rbind(m1 - m) +
   cbind(m2 - m) %*% rbind(m2 - m) +
   cbind(m3 - m) %*% rbind(m3 - m))
B

# covariance within groups (estimate)
Sp

# how many coordinates?
g <- 3
p <- 2
s <- min(g - 1, p)
s

# Matrix Sp^(-1/2) by eigenvalue decomposition
val.Sp <- eigen(Sp)$val
vec.Sp <- eigen(Sp)$vec
invSp.2 <-
   1 / sqrt(val.Sp[1]) * vec.Sp[, 1] %*% t(vec.Sp[, 1]) +
   1 / sqrt(val.Sp[2]) * vec.Sp[, 2] %*% t(vec.Sp[, 2])
invSp.2

# spectral decomposition of Sp^(-1/2) B Sp^(-1/2)
spec.dec <- eigen(invSp.2 %*% B %*% invSp.2)

# first canonical coordinate
a1 <- invSp.2 %*% spec.dec$vec[, 1]
a1

# second canonical coordinate
a2 <- invSp.2 %*% spec.dec$vec[, 2]
a2

# compare with the output of lda():
lda.iris <- lda(iris2, species.name)
lda.iris
a1
a2
spec.dec$val / sum(spec.dec$val)

### How are the data classified?
# Compute the canonical coordinates of the data
cc1.iris <- as.matrix(iris2) %*% a1
cc2.iris <- as.matrix(iris2) %*% a2

coord.cc <- cbind(cc1.iris, cc2.iris)

# Compute the coordinates of the mean within groups along the canonical directions
cc.m1 <- c(m1 %*% a1, m1 %*% a2)
cc.m2 <- c(m2 %*% a1, m2 %*% a2)
cc.m3 <- c(m3 %*% a1, m3 %*% a2)

# Assign data to groups
f.class <- rep(0, n)
for (i in 1:n) { # for each datum
   # Compute the Euclidean distance of the i-th datum from mean within the groups
   dist.m <- c(
      d1 = sqrt(sum((coord.cc[i, ] - cc.m1)^2)),
      d2 = sqrt(sum((coord.cc[i, ] - cc.m2)^2)),
      d3 = sqrt(sum((coord.cc[i, ] - cc.m3)^2))
   )
   # Assign the datum to the group whose mean is the nearest
   f.class[i] <- which.min(dist.m)
}
f.class
table(class.true = species.name, class.assigned = f.class)

errors <- 150 - sum(diag(table(class.true = species.name, class.assigned = f.class)))
errors
length(species.name)

APERf <- errors / length(species.name)
APERf

### How do I classify a new observation?
x.new <- c(5.85, 2.90)
# compute the canonical coordinates
cc.new <- c(x.new %*% a1, x.new %*% a2)
# compute the distance from the means
dist.m <- c(
   d1 = sqrt(sum((cc.new - cc.m1)^2)),
   d2 = sqrt(sum((cc.new - cc.m2)^2)),
   d3 = sqrt(sum((cc.new - cc.m3)^2))
)
# assign to the nearest mean
which.min(dist.m)

# visually
color.species <- rep(c("red", "green", "blue"), each = 50)
par(mfrow = c(1, 2))
plot(iris2[, 1], iris2[, 2],
   main = "Plane of original coordinates",
   xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20, col = as.character(color.species)
)
legend("topleft", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)
points(x.new[1], x.new[2], col = "gold", pch = 19)
points(m1[1], m1[2], pch = 4, col = "red", lwd = 2, cex = 1.5)
points(m2[1], m2[2], pch = 4, col = "green", lwd = 2, cex = 1.5)
points(m3[1], m3[2], pch = 4, col = "blue", lwd = 2, cex = 1.5)

plot(cc1.iris, cc2.iris, main = "Plane of canonical coordinates", xlab = "first canonical coordinate", ylab = "second canonical coordinate", pch = 20, col = as.character(color.species))
legend("topleft", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

points(cc.m1[1], cc.m1[2], pch = 4, col = "red", lwd = 2, cex = 1.5)
points(cc.m2[1], cc.m2[2], pch = 4, col = "green", lwd = 2, cex = 1.5)
points(cc.m3[1], cc.m3[2], pch = 4, col = "blue", lwd = 2, cex = 1.5)

points(cc.new[1], cc.new[2], col = "gold", pch = 19)

segments(cc.m1[1], cc.m1[2], cc.new[1], cc.new[2])
segments(cc.m2[1], cc.m2[2], cc.new[1], cc.new[2])
segments(cc.m3[1], cc.m3[2], cc.new[1], cc.new[2])

dev.off()

### We plot the partition generated by the canonical coordinates
color.species <- species.name
levels(color.species) <- c("red", "green", "blue")


plot(cc1.iris, cc2.iris, main = "Fisher discriminant analysis", xlab = "first canonical coordinate", ylab = "second canonical coordinate", pch = 20, col = as.character(color.species))
legend("topleft", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

points(cc.m1[1], cc.m1[2], pch = 4, col = "red", lwd = 2, cex = 1.5)
points(cc.m2[1], cc.m2[2], pch = 4, col = "green", lwd = 2, cex = 1.5)
points(cc.m3[1], cc.m3[2], pch = 4, col = "blue", lwd = 2, cex = 1.5)

x.cc <- seq(min(cc1.iris), max(cc1.iris), len = 200)
y.cc <- seq(min(cc2.iris), max(cc2.iris), len = 200)
xy.cc <- expand.grid(cc1 = x.cc, cc2 = y.cc)

z <- cbind(sqrt(rowSums(scale(xy.cc, cc.m1, scale = FALSE)^2)), sqrt(rowSums(scale(xy.cc, cc.m2, scale = FALSE)^2)), sqrt(rowSums(scale(xy.cc, cc.m3, scale = FALSE)^2)))
z1.cc <- z[, 1] - pmin(z[, 2], z[, 3])
z2.cc <- z[, 2] - pmin(z[, 1], z[, 3])
z3.cc <- z[, 3] - pmin(z[, 1], z[, 2])

contour(x.cc, y.cc, matrix(z1.cc, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x.cc, y.cc, matrix(z2.cc, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x.cc, y.cc, matrix(z3.cc, 200), levels = 0, drawlabels = FALSE, add = TRUE)

dev.off()

# Plot LDA

plot(iris2, main = "Iris Sepal", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20)
points(iris2[i1, ], col = "red", pch = 20)
points(iris2[i2, ], col = "green", pch = 20)
points(iris2[i3, ], col = "blue", pch = 20)
legend("topright", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

points(rbind(m1, m2, m3), pch = 4, col = c("red", "green", "blue"), lwd = 2, cex = 1.5)

x <- seq(min(iris[, 1]), max(iris[, 1]), length = 200)
y <- seq(min(iris[, 2]), max(iris[, 2]), length = 200)
xy <- expand.grid(Sepal.Length = x, Sepal.Width = y)

z <- predict(lda.iris, xy)$post # these are P_i*f_i(x,y)
z1 <- z[, 1] - pmax(z[, 2], z[, 3]) # P_1*f_1(x,y)-max{P_j*f_j(x,y)}
z2 <- z[, 2] - pmax(z[, 1], z[, 3]) # P_2*f_2(x,y)-max{P_j*f_j(x,y)}
z3 <- z[, 3] - pmax(z[, 1], z[, 2]) # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

# Plot the contour line of level (levels=0) of z1, z2, z3:
# P_i*f_i(x,y)-max{P_j*f_j(x,y)}=0 i.e., boundary between R.i and R.j
# where j realizes the max.
contour(x, y, matrix(z1, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z2, 200), levels = 0, drawlabels = FALSE, add = TRUE)
contour(x, y, matrix(z3, 200), levels = 0, drawlabels = FALSE, add = TRUE)

dev.off()

###########################################################################
# Plot of the projections on the canonical directions (not orthogonal!)
plot(iris2, main = "Projections on the canonical directions", xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20, xlim = c(-3, 8), ylim = c(-3, 7))
points(iris2[i1, ], col = "red", pch = 20)
points(iris2[i2, ], col = "green", pch = 20)
points(iris2[i3, ], col = "blue", pch = 20)
legend("topleft", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

points(rbind(m1, m2, m3), pch = 4, col = c("red", "green", "blue"), lwd = 2, cex = 1.5)

abline(h = 0, v = 0, col = "grey35")

arrows(x0 = 0, y0 = 0, x1 = a1[1], y1 = a1[2], length = .1)
arrows(x0 = 0, y0 = 0, x1 = a2[1], y1 = a2[2], length = .1)

text(a1[1], a1[2], "a1", pos = 3)
text(a2[1], a2[2], "a2", pos = 2)

abline(coef = c(0, (a1[2] / a1[1])), col = "grey55", lty = 2)
abline(coef = c(0, (a2[2] / a2[1])), col = "grey55", lty = 2)

points(cc1.iris * a1[1] / (sum(a1^2)), cc1.iris * a1[2] / (sum(a1^2)), col = as.character(color.species))
points(cc2.iris * a2[1] / (sum(a2^2)), cc2.iris * a2[2] / (sum(a2^2)), col = as.character(color.species))

dev.off()

plot(cc1.iris, cc2.iris, main = "Coordinate system of the canonical coordinates", xlab = "first canonical coordinate", ylab = "second canonical coordinate", pch = 20, col = as.character(color.species))
legend("topleft", legend = levels(species.name), fill = c("red", "green", "blue"), cex = .7)

cc.m1 <- c(m1 %*% a1, m1 %*% a2)
cc.m2 <- c(m2 %*% a1, m2 %*% a2)
cc.m3 <- c(m3 %*% a1, m3 %*% a2)

points(cc.m1[1], cc.m1[2], pch = 4, col = "red", lwd = 2, cex = 1.5)
points(cc.m2[1], cc.m2[2], pch = 4, col = "green", lwd = 2, cex = 1.5)
points(cc.m3[1], cc.m3[2], pch = 4, col = "blue", lwd = 2, cex = 1.5)

graphics.off()



################################################################################################
### Support Vector Machines

rm(list = ls())

library(e1071)
help(svm)

###########
### Linear case

# Generate the data
set.seed(123)
x <- matrix(rnorm(20 * 2), ncol = 2)
y <- c(rep(-1, 10), rep(1, 10))
x[y == 1, ] <- x[y == 1, ] + 1

# The classes are not separable
plot(x,
   col = ifelse(y == 1, "blue", "red"),
   pch = 19, xlab = "x1", ylab = "x2", asp = 1
)

# Fit the Support Vector Classifier (kernel = "linear")
# given a cost C
dat <- data.frame(x = x, y = as.factor(y))
svmfit <- svm(y ~ ., data = dat, kernel = "linear", cost = 10, scale = FALSE)
summary(svmfit)


par(mfrow = c(1, 2))
plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19, asp = 1)

# support vectors are indicated with crosses
# they are:
svmfit$index

###
n.g <- 100

xgrid <- expand.grid(
   x.1 = seq(from = range(dat$x.1)[1], to = range(dat$x.1)[2], length = n.g),
   x.2 = seq(from = range(dat$x.2)[1], to = range(dat$x.2)[2], length = n.g)
)
ygrid <- predict(svmfit, xgrid)

plot(xgrid, col = c("red", "blue")[as.numeric(ygrid)], pch = 20, cex = .2)
points(x, col = c("red", "blue")[as.numeric(dat$y)], pch = 19)
points(x[svmfit$index, ], pch = 5, cex = 2)


plot(x, col = c("red", "blue")[as.numeric(dat$y)], pch = 19)
contour(seq(from = range(dat$x.1)[1], to = range(dat$x.1)[2], length = n.g),
   seq(from = range(dat$x.2)[1], to = range(dat$x.2)[2], length = n.g),
   matrix(as.numeric(ygrid), n.g, n.g),
   level = 1.5, add = TRUE,
   drawlabels = FALSE
)

###

# If we try to change the cost parameter we get more support points
# (higher bias, lower variance)
svmfit <- svm(y ~ ., data = dat, kernel = "linear", cost = 0.1, scale = FALSE)

plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19, asp = 1)

# To set the parameter C we can use the function tune(),
# which is based on cross-validation (10-fold)
set.seed(1)
tune.out <- tune(svm, y ~ .,
   data = dat, kernel = "linear",
   ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100))
)
summary(tune.out)

# Extract the best model from the result of tune
bestmod <- tune.out$best.model
summary(bestmod)

plot(bestmod, dat, col = c("salmon", "light blue"), pch = 19, asp = 1)

# Prediction for a new observation (command predict())
xtest <- matrix(rnorm(20 * 2), ncol = 2)
ytest <- sample(c(-1, 1), 20, rep = TRUE)
xtest[ytest == 1, ] <- xtest[ytest == 1, ] + 1
testdat <- data.frame(x = xtest, y = as.factor(ytest))

plot(xtest,
   col = ifelse(ytest == 1, "light blue", "salmon"),
   pch = 19, xlab = "x1", ylab = "x2", asp = 1
)

ypred <- predict(bestmod, testdat)
table(true.label = testdat$y, assigned.label = ypred)


# If the classes are separable, setting a high value for the cost function
# leads to the maximal margin classifier (i.e., it returns the classification
# provided by the best separating hyperplane)

species.name <- factor(iris$Species, labels = c("setosa", "versicolor", "virginica"))
set.seed(1)
iris2 <- iris[, 1:2] + cbind(rnorm(150, sd = 0.025)) # jittering

# setosa VS versicolor+virginica
y <- rep(0, 150)
y[which(species.name == "setosa")] <- 1


plot(iris2[, 1], iris2[, 2], xlab = "Sepal.Length", ylab = "Sepal.Width", pch = 20, col = as.character(y + 1))

dat <- data.frame(x = iris2[, c(2, 1)], y = as.factor(y))
svmfit <- svm(y ~ ., data = dat, kernel = "linear", cost = 1000, scale = FALSE)
summary(svmfit)


par(mfrow = c(1, 2))
plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19)


dat <- data.frame(x = iris2[, c(2, 1)], y = as.factor(y))
svmfit <- svm(y ~ ., data = dat, kernel = "linear", cost = 1, scale = FALSE)
summary(svmfit)


par(mfrow = c(1, 2))
plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19)


###########
### Non-linear case

# Generate the data
set.seed(1)
x <- matrix(rnorm(200 * 2), ncol = 2)
x[1:100, ] <- x[1:100, ] + 2
x[101:150, ] <- x[101:150, ] - 2
y <- c(rep(1, 150), rep(2, 50))
dat <- data.frame(x = x, y = as.factor(y))

# The classes are not separable
plot(x, col = ifelse(y == 1, "salmon", "light blue"), pch = 19, xlab = "x1", ylab = "x2", asp = 1)

# Randomly split in train and test
train <- sample(200, 100)

# Fit a Support Vector Machine (kernel = "radial") given a cost C
svmfit <- svm(y ~ ., data = dat[train, ], kernel = "radial", gamma = 1, cost = 0.5)
summary(svmfit)

# Plot the SVM
plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19, asp = 1)

# Misclassification error on the training set
table(true = dat[train, "y"], pred = predict(svmfit,
   newdata = dat[train, ]
))
(3 + 4) / 100

# Misclassification error on the test set
table(true = dat[-train, "y"], pred = predict(svmfit,
   newdata = dat[-train, ]
))
(3 + 10) / 100

# Increasing the cost decreases the errors on the training set,
# at the expense of a more irregular boundary
svmfit <- svm(y ~ ., data = dat[train, ], kernel = "radial", gamma = 1, cost = 1e5)


plot(svmfit, dat, col = c("salmon", "light blue"), pch = 19, asp = 1)

# Misclassification error on the training set
table(true = dat[train, "y"], pred = predict(svmfit,
   newdata = dat[train, ]
))
0

# Misclassification error on the test set
table(true = dat[-train, "y"], pred = predict(svmfit,
   newdata = dat[-train, ]
))
(5 + 20) / 100

# Set parameters via CV:
set.seed(1)
tune.out <- tune(svm, y ~ .,
   data = dat[train, ], kernel = "radial",
   ranges = list(
      cost = c(0.1, 1, 10, 100, 1000),
      gamma = c(0.5, 1, 2, 3, 4)
   )
)
summary(tune.out)

# Misclassification error with best model on the training set
table(true = dat[train, "y"], pred = predict(tune.out$best.model,
   newdata = dat[train, ]
))
(2 + 4) / 100

# Misclassification error with best model on the test set
table(true = dat[-train, "y"], pred = predict(tune.out$best.model,
   newdata = dat[-train, ]
))
(2 + 10) / 100



###########
### Application to Gene Expression Data (multiclass classification)
library(ISLR)

help(Khan)

is(Khan)
names(Khan)

dim(Khan$xtrain)
dim(Khan$xtest)
length(Khan$ytrain)
length(Khan$ytest)

table(Khan$ytrain)
table(Khan$ytest)

dat <- data.frame(x = Khan$xtrain, y = as.factor(Khan$ytrain))
out <- svm(y ~ ., data = dat, kernel = "linear", cost = 10)
summary(out)

table(out$fitted, dat$y)

dat.te <- data.frame(x = Khan$xtest, y = as.factor(Khan$ytest))
pred.te <- predict(out, newdata = dat.te)
table(pred.te, dat.te$y)
