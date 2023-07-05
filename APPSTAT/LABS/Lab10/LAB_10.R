### ---------------------###
### LAB 10 (16/05/2023)  ###
### ---------------------###

### TOPICS:
### Linear models

library(MASS)
library(car)
library(rgl)

library(glmnet)

options(rgl.printRglwidget = TRUE)

# _______________________________________________________________________________
##### Problem of colinearity
##### --------------------------

### Example 1: Multiple linear regression
### ----------------------------------------
### Dataset cars: distance taken to stop [ft] as a function of velocity [mph]
### for some cars in the 1920s

cars

plot(cars,
  xlab = "Speed",
  ylab = "Stopping distance",
  las = 1
)

n <- dim(cars)[[1]]
distance <- cars$dist
speed1 <- cars$speed
speed2 <- cars$speed^2

### Model:
### distance = beta_0 + beta_1 * speed + beta_2 * speed^2 + Eps
### (linear in the parameters!)

fm <- lm(distance ~ speed1 + speed2)
summary(fm)
# Note: the overall model is significant (F-test), but the individual
#       coefficients are not significant (t-test) -> colinearity
# Note: colinearity

# Variance inflation factor
help(vif)

vif(fm) # Rule of thumb -> problem when VIF exceeds 10 (or 5 sometimes)
# in this case VIF is ~ 24.5 for both regressors -> problem

# recall vif formula:
1 / (1 - summary(lm(speed1 ~ speed2))$r.squared)
1 / (1 - summary(lm(speed2 ~ speed1))$r.squared)

### A possible solution to colinearity: PCA
speed.pc <- princomp(cbind(speed1, speed2), scores = TRUE)
summary(speed.pc)
# Note: the first PC explains 99.99% of the variance
speed.pc$load
# Note: the loadings of the first principal component comes all from speed^2
#       and we're explaining 99.99% of the variance with it.

speed.pc$scores

sp1.pc <- speed.pc$scores[, 1] # first PC (IT'S BASICALLY SPEED^2!)
sp2.pc <- speed.pc$scores[, 2] # second PC (IT'S BASICALLY SPEED!)

# Now we estimate the model by inserting the PCs instead of the
# original regressors
# Model: y = b0 + b1*PC1+ b2*PC2 + eps, eps~N(0,sigma^2)
fm.pc <- lm(distance ~ sp1.pc + sp2.pc)
summary(fm.pc)
# The first component is the only one significant to the model -> speed^2

# Note: we could have performed dimensionality reduction before
# estimating the model and then considered only the first PC.

# We can re-write the model as:
# Model:
# y= b0 + b1*      PC1                 + b2*      PC2                 + eps =
#  = b0 + b1*(e11*(Alluminium-m1)+e21*(Silicate-m2)) + b2*(e12*(Alluminium-m1)+e22*(Silicate-m2)) + eps =
#  = b0 - b1*e11*m1 - b2*e12*m1 - b1*e21*m2 - b2*e22*m2 +
#                           + (b1*e11+b2*e12)*Alluminium + (b1*e21+b2*e22)*Silicate + eps
# where e.ij are the loadings, i=1,2, j=1,2.
# => We can compute the coefficients of the model which used the original
#    regressors
m1 <- mean(speed1)
m2 <- mean(speed2)

beta0 <- coefficients(fm.pc)[1] -
  coefficients(fm.pc)[2] * speed.pc$load[1, 1] * m1 -
  coefficients(fm.pc)[3] * speed.pc$load[1, 2] * m1 -
  coefficients(fm.pc)[2] * speed.pc$load[2, 1] * m2 -
  coefficients(fm.pc)[3] * speed.pc$load[2, 2] * m2
beta1 <- coefficients(fm.pc)[2] * speed.pc$load[1, 1] +
  coefficients(fm.pc)[3] * speed.pc$load[1, 2]
beta2 <- coefficients(fm.pc)[2] * speed.pc$load[2, 1] +
  coefficients(fm.pc)[3] * speed.pc$load[2, 2]

c(
  beta0 = as.numeric(beta0),
  beta1 = as.numeric(beta1),
  beta2 = as.numeric(beta2)
)
fm$coefficients
# and in fact as we can see the coefficients are not the same: the first
# coefficient, related to speed, is smaller in the model with the PCs, while the
# second coefficient, related to speed^2, is bigger in the model with the PCs.

x <- seq(0, 25, len = 100)
plot(cars,
  xlab = "Speed",
  ylab = "Stopping distance",
  las = 1
)
lines(x, beta0 + beta1 * x + beta2 * x^2)

# Reduce the model:
fm.pc <- lm(distance ~ sp1.pc)
summary(fm.pc)

# We can re-write the model as:
# Model: y = b0 + b1*      PC1                 + eps =
#          = b0 + b1*(e11*(Alluminium-m1)1+e21*(Silicate-m2)) + eps =
#          = b0 - b1*e11*m1 - b2*e21*m2 + b1*e11*Alluminium + b1*e21*Silicate + eps
beta0 <- coefficients(fm.pc)[1] -
  coefficients(fm.pc)[2] * speed.pc$load[1, 1] * m1 -
  coefficients(fm.pc)[2] * speed.pc$load[2, 1] * m2
beta1 <- coefficients(fm.pc)[2] * speed.pc$load[1, 1]
beta2 <- coefficients(fm.pc)[2] * speed.pc$load[2, 1]

c(
  beta0 = as.numeric(beta0),
  beta1 = as.numeric(beta1),
  beta2 = as.numeric(beta2)
)
fm$coefficients

plot(
  sp1.pc,
  distance,
  xlab = "PC1",
  ylab = "Stopping distance",
  las = 1,
  xlim = c(-250, 361),
  ylim = c(-5, 130)
)
x <- seq(-250, 361, by = 1)
b <- coef(fm.pc)
lines(x, b[1] + b[2] * x)

plot(
  speed1,
  distance,
  ylab = "Stopping distance",
  las = 1,
  ylim = c(-5, 130)
)
x <- seq(0, 25, by = 1)
lines(x, beta0 + beta1 * x + beta2 * x^2)

# diagnostics of the residuals
par(mfrow = c(2, 2))
plot(fm.pc)

shapiro.test(residuals(fm.pc))

dev.off()

# The first PC is basically speed2 (centered with respect to the mean),
# hence we could just consider the regressor speed2

fm.2 <- lm(distance ~ speed2)
summary(fm.2)
# Note: when using as regressors polynomials of the same variable (i.e., in
# the framework of polynomial regression), some authors recommend instead to
# keep all the orders lower than the maximum that one wants to keep (e.g., if
# keeping the maximum order 3, they recommend keeping also all the terms
# of order 2 and 1, for a "hierarchy" principle)

plot(
  speed2,
  distance,
  xlab = "speed2",
  ylab = "Stopping distance",
  las = 1,
  xlim = c(15, 626),
  ylim = c(-5, 130)
)
x <- seq(0, 650, by = 1)
b <- coef(fm.2)
lines(x, b[1] + b[2] * x)

plot(
  cars,
  xlab = "Speed",
  ylab = "Stopping distance",
  las = 1,
  xlim = c(0, 30),
  ylim = c(-5, 130)
)
x <- seq(0, 30, by = 0.1)
b <- coef(fm.2)
lines(x, b[1] + b[2] * x^2)

dev.off()

### Otherwise, another possible solution to colinearity: ridge regression
library(MASS)
help(lm.ridge)

# Fix lambda
lambda <- .5
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda)
# Note: to fit the model, R automatically centers X and y
# with respect to their mean.

coef.ridge <- coef(fit.ridge)
yhat.lm <-
  cbind(rep(1, n), speed1, speed2) %*% coef(fm) # LM fitted values
yhat.r <-
  cbind(rep(1, n), speed1, speed2) %*% coef.ridge # ridge fitted values

plot(
  speed1,
  yhat.lm,
  type = "l",
  lty = 4,
  lwd = 2,
  ylab = "Distance",
  xlab = "Speed"
)
points(speed1, distance, pch = 1, cex = .8)
matlines(
  speed1,
  yhat.r,
  type = "l",
  lty = 1,
  col = grey.colors(length(lambda)),
  lwd = 2
)
legend(
  "topleft",
  c("lm", "ridge"),
  lty = c(4, 1),
  col = c("black", grey.colors(length(lambda))),
  lwd = 2
)


# Repeat for a grid of lambda's
lambda.c <- seq(0, 10, 0.01)
fit.ridge <- lm.ridge(distance ~ speed1 + speed2, lambda = lambda.c)


par(mfrow = c(1, 3))
plot(
  lambda.c,
  coef(fit.ridge)[, 1],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[0])
)
abline(h = coef(fm)[1], lty = 2)
plot(
  lambda.c,
  coef(fit.ridge)[, 2],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[1])
)
abline(h = coef(fm)[2], lty = 2)
plot(
  lambda.c,
  coef(fit.ridge)[, 3],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[2])
)
abline(h = coef(fm)[3], lty = 2)

dev.off()

yhat.lm <- cbind(rep(1, n), speed1, speed2) %*% coef(fm)

plot(
  speed1,
  yhat.lm,
  type = "l",
  lty = 1,
  lwd = 2,
  ylab = "Distance",
  xlab = "Speed"
)
points(speed1, distance, pch = 1, cex = .8)
yhat.r <- NULL
for (i in 1:length(lambda.c)) {
  yhat.r <-
    cbind(yhat.r, cbind(rep(1, n), speed1, speed2) %*% coef(fit.ridge)[i, ])
}
matlines(speed1,
  yhat.r,
  type = "l",
  lty = 1,
  col = grey.colors(length(lambda.c))
)
lines(
  speed1,
  yhat.lm,
  type = "l",
  lty = 4,
  lwd = 2,
  ylab = "Distance",
  xlab = "Speed"
)

# Choice of the optimal lambda, e.g., via cross-validation
# GCV -> Generalized Cross Validation (a sort of LOO CV)
select(fit.ridge)

# or
lambda.opt <- lambda.c[which.min(fit.ridge$GCV)]
lambda.opt

par(mfrow = c(1, 3))
plot(
  lambda.c,
  coef(fit.ridge)[, 1],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[0])
)
abline(
  h = coef(fm)[1],
  lty = 1,
  col = "grey"
)
abline(v = lambda.opt, col = 2, lty = 2)
plot(
  lambda.c,
  coef(fit.ridge)[, 2],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[1])
)
abline(
  h = coef(fm)[2],
  lty = 1,
  col = "grey"
)
abline(v = lambda.opt, col = 2, lty = 2)
plot(
  lambda.c,
  coef(fit.ridge)[, 3],
  type = "l",
  xlab = expression(lambda),
  ylab = expression(beta[2])
)
abline(
  h = coef(fm)[3],
  lty = 1,
  col = "grey"
)
abline(v = lambda.opt, col = 2, lty = 2)

dev.off()

plot(
  speed1,
  distance,
  pch = 1,
  cex = .8,
  ylab = "Distance",
  xlab = "Speed"
)
matlines(speed1,
  yhat.r,
  type = "l",
  lty = 1,
  col = grey.colors(length(lambda.c))
)
lines(
  speed1,
  yhat.lm,
  type = "l",
  lty = 4,
  lwd = 2,
  ylab = "Distance",
  xlab = "Speed"
)
lines(
  speed1,
  yhat.r[, which.min(fit.ridge$GCV)],
  type = "l",
  lty = 1,
  lwd = 2,
  col = 2,
  ylab = "Distance",
  xlab = "Speed"
)
legend(
  "topleft",
  c("LM", "Ridge opt."),
  lty = c(4, 1),
  col = c(1, 2),
  lwd = 2
)

coef.ridge <- coef(fit.ridge)[which.min(fit.ridge$GCV), ]
coef.ridge

### Otherwise, another possible solution to colinearity, which also performs
# variables selection: lasso regression
library(glmnet)
help(glmnet)

# Build the matrix of predictors
x <- model.matrix(distance ~ speed1 + speed2)[, -1]
# Build the vector of response
y <- distance

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- 10^seq(5, -3, length = 100)
fit.lasso <-
  glmnet(x, y, lambda = lambda.grid) # default: alpha=1 -> lasso
# [note: if alpha=0 -> ridge regression]

plot(fit.lasso,
  xvar = "lambda",
  label = TRUE,
  col = rainbow(dim(x)[2])
)
legend(
  "topright",
  dimnames(x)[[2]],
  col = rainbow(dim(x)[2]),
  lty = 1,
  cex = 1
)

# Let's set lambda via cross validation
cv.lasso <-
  cv.glmnet(x, y, lambda = lambda.grid) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

plot(cv.lasso)
abline(v = log(bestlam.lasso), lty = 1)

# Get the coefficients for the optimal lambda
coef.lasso <-
  predict(fit.lasso, s = bestlam.lasso, type = "coefficients")
coef.lasso


# ______________________________________________________________________________________________________________________________________________________________
##### Multiple linear regression
##### --------------------------------------------------------
setwd("~/HPC/APPSTAT/Lab10")
data <- read.table("concrete.txt", header = TRUE)

head(data)
dim(data)
names(data)

pairs(data)

attach(data)

y <- Hardness_concrete

# detach(data)

par(mfrow = c(2, 2))
plot(
  Alluminium,
  y,
  main = "Hardness vs Alluminium",
  lwd = 2,
  xlab = "Alluminium",
  ylab = "Hardness concrete"
)
plot(
  Silicate,
  y,
  main = "Hardness vs Silicate",
  lwd = 2,
  xlab = "Silicate",
  ylab = "Hardness concrete"
)
plot(
  Alluminium_ferrite,
  y,
  main = "Hardness vs Alluminium ferrite",
  lwd = 2,
  xlab = "Alluminium ferrite",
  ylab = "Hardness concrete"
)
plot(
  Silicate_bicalcium,
  y,
  main = "Hardness vs Silicate bicalcium",
  lwd = 2,
  xlab = "Silicate bicalcium",
  ylab = "Hardness concrete"
)

## Multiple linear regression

result <-
  lm(y ~ Alluminium + Silicate + Alluminium_ferrite + Silicate_bicalcium)
summary(result)

library(car)
vif(result)

help(princomp)

#### PCA regression
result.pc <-
  princomp(cbind(Alluminium, Silicate, Alluminium_ferrite, Silicate_bicalcium),
    scores = TRUE
  )
summary(result.pc)
result.pc$load
# the first principal component show a contrast between Silicate and
# Silicate_bicalcium, while the second principal component show a contrast
# between Alluminium and Alluminium_ferrite

# Explained variance
layout(matrix(c(2, 3, 1, 3), 2, byrow = TRUE))
barplot(result.pc$sdev^2,
  las = 2,
  main = "Principal Components",
  ylab = "Variances"
)
barplot(
  c(
    sd(Alluminium),
    sd(Silicate),
    sd(Alluminium_ferrite),
    sd(Silicate_bicalcium)
  )^2,
  las = 2,
  main = "Original variables",
  ylab = "Variances"
)
plot(
  cumsum(result.pc$sdev^2) / sum(result.pc$sde^2),
  type = "b",
  axes = FALSE,
  xlab = "number of components",
  ylab = "contribution to the total variance",
  ylim = c(0, 1)
)
abline(h = 1, col = "blue")
abline(h = 0.9, lty = 2, col = "blue")
box()
axis(2, at = 0:10 / 10, labels = 0:10 / 10)
axis(1,
  at = 1:4,
  labels = 1:4,
  las = 2
)

materials <- cbind(
  Alluminium,
  Silicate,
  Alluminium_ferrite,
  Silicate_bicalcium
)
materials.scaled <- scale(materials, center = TRUE, scale = TRUE)

result.pc.scaled <-
  princomp(materials.scaled, scores = TRUE)
summary(result.pc.scaled)
result.pc.scaled$load

# Explained variance
layout(matrix(c(2, 3, 1, 3), 2, byrow = TRUE))
barplot(result.pc.scaled$sdev^2,
  las = 2,
  main = "Principal Components",
  ylab = "Variances"
)
barplot(
  c(
    sd(materials.scaled[, 1]),
    sd(materials.scaled[, 2]),
    sd(materials.scaled[, 3]),
    sd(materials.scaled[, 4])
  )^2,
  las = 2,
  main = "Original (scaled) variables",
  ylab = "Variances"
)
plot(
  cumsum(result.pc.scaled$sdev^2) / sum(result.pc.scaled$sdev^2),
  type = "b",
  axes = FALSE,
  xlab = "number of components",
  ylab = "contribution to the total variance",
  ylim = c(0, 1)
)
abline(h = 1, col = "blue")
abline(h = 0.9, lty = 2, col = "blue")
box()
axis(2, at = 0:10 / 10, labels = 0:10 / 10)
axis(1,
  at = 1:4,
  labels = 1:4,
  las = 2
)

# would still choose the first two components

# Loadings
par(mar = c(1, 4, 0, 2), mfcol = c(4, 2))
for (i in 1:4) {
  barplot(result.pc$load[, i], ylim = c(-1, 1))
}

# par(mar = c(1, 4, 0, 2), mfrow = c(4, 1))
for (i in 1:4) {
  barplot(result.pc.scaled$load[, i], ylim = c(-1, 1))
}

dev.off()

# Dimensionality reduction: select first two PCs:
pc1 <- result.pc$scores[, 1]
pc2 <- result.pc$scores[, 2]

# Now we estimate the model using the first two PCs as regressors.
# Model: y = b0 + b1*PC1+ b2*PC2 + eps, eps~N(0,sigma^2)
fm.pc <- lm(y ~ pc1 + pc2)

summary(fm.pc)

# We can re-write the model as:
# Model:
# y= b0 + b1*PC1 + b2*PC2 + eps =
#  = b0 + b1*(e11*(Alluminium-m1)+e21*(Silicate-m2)+e31*(Alluminium_ferrite-m3)+e41*(Silicate_bicalcium-m4)) +
#       + b2*(e12*(Alluminium-m1)+e22*(Silicate-m2)+e32*(Alluminium_ferrite-m3)+e42*(Silicate_bicalcium-m4)) + eps =
#  = b0 - b1*e11*m1 - b2*e12*m1 - b1*e21*m2 - b2*e22*m2 +
#       - b1*e31*m3 - b2*e32*m3 - b1*e41*m4 - b2*e42*m4 +
#       + (b1*e11+b2*e12)*Alluminium + (b1*e21+b2*e22)*Silicate +
#       + (b1*e31+b2*e32)*Alluminium_ferrite + (b1*e41+b2*e42)*Silicate_bicalcium + eps
# where e.ij are the loadings, i=1,2,3,4, j=1,2.
# => We can compute the coefficients of the model which used the original
#    regressors
m1 <- mean(Alluminium)
m2 <- mean(Silicate)
m3 <- mean(Alluminium_ferrite)
m4 <- mean(Silicate_bicalcium)
beta0 <- coefficients(fm.pc)[1] -
  coefficients(fm.pc)[2] * result.pc$load[1, 1] * m1 -
  coefficients(fm.pc)[3] * result.pc$load[1, 2] * m1 -
  coefficients(fm.pc)[2] * result.pc$load[2, 1] * m2 -
  coefficients(fm.pc)[3] * result.pc$load[2, 2] * m2 -
  coefficients(fm.pc)[2] * result.pc$load[3, 1] * m3 -
  coefficients(fm.pc)[3] * result.pc$load[3, 2] * m3 -
  coefficients(fm.pc)[2] * result.pc$load[4, 1] * m4 -
  coefficients(fm.pc)[3] * result.pc$load[4, 2] * m4
beta1 <- coefficients(fm.pc)[2] * result.pc$load[1, 1] +
  coefficients(fm.pc)[3] * result.pc$load[1, 2]
beta2 <- coefficients(fm.pc)[2] * result.pc$load[2, 1] +
  coefficients(fm.pc)[3] * result.pc$load[2, 2]
beta3 <- coefficients(fm.pc)[2] * result.pc$load[3, 1] +
  coefficients(fm.pc)[3] * result.pc$load[3, 2]
beta4 <- coefficients(fm.pc)[2] * result.pc$load[4, 1] +
  coefficients(fm.pc)[3] * result.pc$load[4, 2]

c(
  beta0 = as.numeric(beta0),
  beta1 = as.numeric(beta1),
  beta2 = as.numeric(beta2),
  beta3 = as.numeric(beta3),
  beta4 = as.numeric(beta4)
)
result$coefficients

# diagnostics of the residuals
par(mfrow = c(2, 2))
plot(fm.pc)

shapiro.test(residuals(fm.pc))
# residuals are normal

dev.off()

### Ridge and Lasso regression with glmnet

x <-
  model.matrix(y ~ Alluminium + Silicate + Alluminium_ferrite + Silicate_bicalcium)[, -1] # matrix of predictors
y <- y # vector of response
lambda.grid <- 10^seq(5, -3, length = 50)

# Ridge regression
fit.ridge <-
  glmnet(x, y, lambda = lambda.grid, alpha = 0) # alpha=0 -> ridge
# Also possible to let glmnet compute its own lambda sequence (length 100,
# between lambda.max*10^-4 and lambda.max where lambda.max is derived from
# the data)

plot(fit.ridge,
  xvar = "lambda",
  label = TRUE,
  col = rainbow(dim(x)[2])
)
legend(
  "topright",
  dimnames(x)[[2]],
  col = rainbow(dim(x)[2]),
  lty = 1,
  cex = 1
)

norm_l2 <- NULL
for (i in 1:50) {
  norm_l2 <- c(norm_l2, sqrt(sum((fit.ridge$beta[, i])^2)))
}

plot(log(lambda.grid), norm_l2)

# Let's set lambda via CV
set.seed(1)
cv.ridge <-
  cv.glmnet(x,
    y,
    alpha = 0,
    nfolds = 3,
    lambda = lambda.grid
  )

bestlam.ridge <- cv.ridge$lambda.min
bestlam.ridge

plot(cv.ridge)
abline(v = log(bestlam.ridge), lty = 1)

# Get the coefficients for the optimal lambda
coef.ridge <-
  predict(fit.ridge, s = bestlam.ridge, type = "coefficients")[1:5, ]
coef.ridge

plot(fit.ridge,
  xvar = "lambda",
  label = TRUE,
  col = rainbow(dim(x)[2])
)
abline(v = log(bestlam.ridge))


### Lasso regression

fit.lasso <-
  glmnet(x, y, lambda = lambda.grid, alpha = 1) # alpha=1 -> lasso

plot(fit.lasso,
  xvar = "lambda",
  label = TRUE,
  col = rainbow(dim(x)[2])
)
legend(
  "topright",
  dimnames(x)[[2]],
  col = rainbow(dim(x)[2]),
  lty = 1,
  cex = 1
)

norm_l1 <- NULL
for (i in 1:50) {
  norm_l1 <- c(norm_l1, sum(abs(fit.ridge$beta[, i])))
}

plot(log(lambda.grid), norm_l1)

# Let's set lambda via CV
set.seed(1)
cv.lasso <-
  cv.glmnet(x,
    y,
    alpha = 1,
    nfolds = 3,
    lambda = lambda.grid
  )

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

plot(cv.lasso)
abline(v = log(bestlam.lasso), lty = 1)

# Get the coefficients for the optimal lambda
coef.lasso <-
  predict(fit.lasso, s = bestlam.lasso, type = "coefficients")[1:5, ]
coef.lasso

plot(fit.lasso,
  xvar = "lambda",
  label = TRUE,
  col = rainbow(dim(x)[2])
)
abline(v = log(bestlam.lasso))


# Compare coefficients estimates for LS, Ridge and Lasso
plot(
  rep(0, dim(x)[2]),
  coef(lm(y ~ x))[-1],
  col = rainbow(dim(x)[2]),
  pch = 20,
  xlim = c(-1, 3),
  ylim = c(-1, 2),
  xlab = "",
  ylab = expression(beta),
  axes = FALSE
)
points(rep(1, dim(x)[2]),
  coef.ridge[-1],
  col = rainbow(dim(x)[2]),
  pch = 20
)
points(rep(2, dim(x)[2]),
  coef.lasso[-1],
  col = rainbow(dim(x)[2]),
  pch = 20
)
abline(
  h = 0,
  col = "grey41",
  lty = 1
)
box()
axis(2)
axis(1,
  at = c(0, 1, 2),
  labels = c("LS", "Ridge", "Lasso")
)
legend(
  "topright",
  dimnames(x)[[2]],
  col = rainbow(dim(x)[2]),
  pch = 20,
  cex = 1
)

# l2 norm
sqrt(sum((coef(lm(y ~ x))[-1])^2)) # LS
sqrt(sum((coef.ridge[-1])^2)) # ridge

# l1 norm
sum(abs(coef(lm(y ~ x))[-1])) # LS
sum(abs(coef.lasso[-1])) # lasso


### Variable selection
result <-
  lm(y ~ Alluminium + Silicate + Alluminium_ferrite + Silicate_bicalcium)
summary(result)

result1 <-
  lm(y ~ Alluminium + Alluminium_ferrite + Silicate_bicalcium)
summary(result1)

result2 <- lm(y ~ Alluminium + Silicate_bicalcium)
summary(result2)

# By doing this I'm not reducing too much the R^2, so I'm not losing too much
# information, but I'm reducing the number of variables, so I'm reducing the
# complexity of the model, and hence it is preferable.

# alternatively I can test the significance of the coefficients of Silicate and
# Alluminium_ferrite: is it okay to remove them? can I consider the reduced
# model (they're both zero) or do I have to keep them (at least one of them
# is not zero)?
linearHypothesis(result, rbind(c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0)), c(0, 0))
# pvalue is 0.19, so I cannot reject the null hypothesis that the coefficients
# are 0 -> I can consider them as zero and remove them from the model



# diagnostics of the residuals
par(mfrow = c(2, 2))
plot(result2)

shapiro.test(residuals(result2))


# _____________________________________________________________________________
##### Model / Variable Selection
##### -------------

library(ISLR)

# Hitters dataset
help(Hitters)
names(Hitters)
dim(Hitters)

# remove NA's
sum(is.na(Hitters$Salary))
Hitters <- na.omit(Hitters)
dim(Hitters)
sum(is.na(Hitters))

###  Subset Selection Methods
library(leaps)

help(regsubsets)

# Best Subset Selection: Exhaustive Search
regfit.full <- regsubsets(Salary ~ ., data = Hitters)
summary(regfit.full)
# for each row i, it tells us which variables are included in the model
# if I'm considering a model with only i variables. For example, the first
# row tells us that the best model with only one variable is the one with
# only CRBI, the second row tells us that the best model with only two
# variables is the one with CRBI and Hits, and so on.

regfit.full <- regsubsets(Salary ~ ., data = Hitters, nvmax = 19)
summary(regfit.full)

reg.summary <- summary(regfit.full)
names(reg.summary)

reg.summary$which

reg.summary$rsq # r-squared
reg.summary$adjr2 # adjusted r-squared
reg.summary$rss # residual sum of squares


par(mfrow = c(1, 3))
plot(reg.summary$rsq,
  xlab = "Number of Variables",
  ylab = "R-squared",
  type = "b"
)
plot(reg.summary$adjr2,
  xlab = "Number of Variables",
  ylab = "Adjusted RSq",
  type = "b"
)
plot(reg.summary$rss,
  xlab = "Number of Variables",
  ylab = "RSS",
  type = "b"
)

# extract coefficient estimates associated with the models
max.adjr2.idx <- which.max(reg.summary$adjr2)
coef(regfit.full, max.adjr2.idx)

coef(regfit.full, 6)

# graphical table of best subsets
help(plot.regsubsets)

plot(regfit.full, scale = "r2", main = "Exhaustive search")

plot(regfit.full, scale = "adjr2", main = "Exhaustive search")


# Forward and Backward Stepwise Selection
regfit.fwd <-
  regsubsets(Salary ~ .,
    data = Hitters,
    nvmax = 19,
    method = "forward"
  )
summary(regfit.fwd)


par(mfrow = c(1, 3))
plot(
  summary(regfit.fwd)$rsq,
  xlab = "Number of Variables",
  ylab = "R-squared",
  type = "b"
)
plot(
  summary(regfit.fwd)$adjr2,
  xlab = "Number of Variables",
  ylab = "Adjusted RSq",
  type = "b"
)
plot(
  summary(regfit.fwd)$rss,
  xlab = "Number of Variables",
  ylab = "RSS",
  type = "b"
)

plot(regfit.fwd, scale = "r2", main = "Forward Stepwise Selection")

plot(regfit.fwd, scale = "adjr2", main = "Forward Stepwise Selection")

# Backward Stepwise Selection
regfit.bwd <-
  regsubsets(Salary ~ .,
    data = Hitters,
    nvmax = 19,
    method = "backward"
  )
summary(regfit.bwd)


par(mfrow = c(1, 3))
plot(
  summary(regfit.bwd)$rsq,
  xlab = "Number of Variables",
  ylab = "R-squared",
  type = "b"
)
plot(
  summary(regfit.bwd)$adjr2,
  xlab = "Number of Variables",
  ylab = "Adjusted RSq",
  type = "b"
)
plot(
  summary(regfit.bwd)$rss,
  xlab = "Number of Variables",
  ylab = "RSS",
  type = "b"
)

plot(regfit.bwd, scale = "r2", main = "Backward Stepwise Selection")

plot(regfit.bwd, scale = "adjr2", main = "Backward Stepwise Selection")


coef(regfit.full, 7) # Exhaustive search
coef(regfit.fwd, 7) # Forward Stepwise Selection
coef(regfit.bwd, 7) # Backward Stepwise Selection

graphics.off()

### Choosing among models using the k-fold cross-validation approach
### (exhaustive search)
k <- 10

set.seed(1)
folds <- sample(1:k, nrow(Hitters), replace = TRUE)
folds
table(folds)

# function that performs the prediction for regsubsets
predict.regsubsets <- function(object, newdata, id) {
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id = id)
  xvars <- names(coefi)
  mat[, xvars] %*% coefi
}

cv.errors <- matrix(NA, k, 19, dimnames = list(NULL, paste(1:19)))

for (j in 1:k) {
  # First remove the data, then perform selection procedure and so on
  best.fit <-
    regsubsets(Salary ~ ., data = Hitters[folds != j, ], nvmax = 19)
  for (i in 1:19) {
    pred <- predict.regsubsets(best.fit, Hitters[folds == j, ], id = i)
    cv.errors[j, i] <- mean((Hitters$Salary[folds == j] - pred)^2)
  }
}

cv.errors
root.mean.cv.errors <-
  sqrt(apply(cv.errors, 2, mean)) # average over the columns
root.mean.cv.errors

plot(root.mean.cv.errors, type = "b")

which.min(root.mean.cv.errors)
points(
  which.min(root.mean.cv.errors),
  root.mean.cv.errors[which.min(root.mean.cv.errors)],
  col = "red",
  pch = 19
)

# estimation on the full dataset
reg.best <- regsubsets(Salary ~ ., data = Hitters, nvmax = 19)
coef(reg.best, 10)


### k-fold cross validation after model selection (WRONG WAY!)
best.fit <- regsubsets(Salary ~ ., data = Hitters, nvmax = 19)
summary(best.fit)

cv.errors_wrong <-
  matrix(NA, k, 19, dimnames = list(NULL, paste(1:19)))

for (j in 1:k) {
  # first perform selection procedure, then do cross validation on those
  for (i in 1:19) {
    covariate <- which(summary(best.fit)$which[i, -c(1, 19)])
    mod <-
      lm(Salary ~ ., data = Hitters[folds != j, c(covariate, 19)])
    pred <- predict(mod, Hitters)[folds == j]
    cv.errors_wrong[j, i] <-
      mean((Hitters$Salary[folds == j] - pred)^2)
  }
}

cv.errors_wrong
root.mean.cv.errors_wrong <-
  sqrt(apply(cv.errors_wrong, 2, mean)) # average over the columns
root.mean.cv.errors_wrong

plot(root.mean.cv.errors_wrong, type = "b")

which.min(root.mean.cv.errors_wrong)
points(
  which.min(root.mean.cv.errors_wrong),
  root.mean.cv.errors_wrong[which.min(root.mean.cv.errors_wrong)],
  col = "red",
  pch = 19
)

points(root.mean.cv.errors, type = "b", col = "blue")


### Ridge and Lasso regression with glmnet
x <- model.matrix(Salary ~ ., Hitters)[, -1] # predictor matrix
y <- Hitters$Salary # response
grid <- 10^seq(10, -2, length = 100) # grid of lambda

# Ridge regression

ridge.mod <- glmnet(x, y, alpha = 0, lambda = grid)

plot(ridge.mod, xvar = "lambda", label = TRUE)

# choosing the parameter lambda
set.seed(123)
cv.out <- cv.glmnet(x,
  y,
  alpha = 0,
  nfold = 3,
  lambda = grid
)

plot(cv.out)

bestlam.ridge <- cv.out$lambda.min
bestlam.ridge
log(bestlam.ridge)

abline(v = log(bestlam.ridge))

plot(ridge.mod, xvar = "lambda", label = TRUE)
abline(v = log(bestlam.ridge))

# Lasso regression

lasso.mod <- glmnet(x, y, alpha = 1, lambda = grid)

plot(lasso.mod, xvar = "lambda", label = TRUE)

# choosing the parameter lambda
set.seed(123)
cv.out <- cv.glmnet(x,
  y,
  alpha = 1,
  nfold = 3,
  lambda = grid
)

plot(cv.out)

bestlam.lasso <- cv.out$lambda.min
bestlam.lasso
log(bestlam.lasso)

abline(v = log(bestlam.lasso))

plot(lasso.mod, xvar = "lambda", label = TRUE)
abline(v = log(bestlam.lasso))


# Compare coefficients estimates for LS, Ridge and Lasso
coef.ridge <-
  predict(ridge.mod, s = bestlam.ridge, type = "coefficients")[1:20, ]
coef.lasso <-
  predict(lasso.mod, s = bestlam.lasso, type = "coefficients")[1:20, ]
coef.ridge
coef.lasso

plot(
  rep(0, dim(x)[2]),
  coef(lm(y ~ x))[-1],
  col = rainbow(dim(x)[2]),
  pch = 20,
  xlim = c(-1, 3),
  ylim = c(-1, 2),
  xlab = "",
  ylab = expression(beta),
  axes = FALSE
)
points(rep(1, dim(x)[2]),
  coef.ridge[-1],
  col = rainbow(dim(x)[2]),
  pch = 20
)
points(rep(2, dim(x)[2]),
  coef.lasso[-1],
  col = rainbow(dim(x)[2]),
  pch = 20
)
abline(
  h = 0,
  col = "grey41",
  lty = 1
)
box()
axis(2)
axis(1,
  at = c(0, 1, 2),
  labels = c("LS", "Ridge", "Lasso")
)

# l2 norm
sqrt(sum((coef(lm(
  y ~ x
))[-1])^2)) # LS
sqrt(sum((coef.ridge[-1])^2)) # ridge

# l1 norm
sum(abs(coef(lm(y ~ x))[-1])) # LS
sum(abs(coef.lasso[-1])) # lasso


##################################################################################
## Logistic Regression ##
help(glm)

##  Pb2 of 23/07/2009  ##
# A store of an appliance chain sells two types of televisions: with 4:3
# screen and with 16:9 screen. The TV.txt file shows the number of
# televisions sold annually for both types of televisions from 1999
# to 2008. By introducing an appropriate logistic model and estimating
# the parameters with the maximum likelihood method:
# a) comment the residual deviance by comparing it with that of the model
#    without regressor.
# Assuming that the store is representative of the Italian situation:
# b) provide a pointwise estimate for the proportion of 16:9 televisions
#    sold in Italy in 2009;
# c) provide a pointwise estimate for the year in which sales of 16:9
#    televisions exceeded those of 4:3;
# d) provide a pointwise estimate for the year in which 16:9 televisions
#    will cover (or have covered) 99% of the Italian market.

TV <- read.table("TV.txt", header = TRUE)
head(TV)
dim(TV)

# a)
fit <- glm(factor(Tipo) ~ Anno, data = TV, family = "binomial")
summary(fit)

plot(TV$Anno, as.numeric(factor(TV$Tipo)) - 1)
lines(
  seq(1997, 2010, by = 0.1),
  predict(fit, data.frame(Anno = seq(1997, 2010, by = 0.1)), type = "response")
)

# null model
Freq.tot <-
  table(TV$Tipo)[2] / (table(TV$Tipo)[1] + table(TV$Tipo)[2])
abline(h = Freq.tot, col = "blue", lty = 2)

# b)
predict(fit, data.frame(Anno = 2009), type = "response")

# c)
(log(0.5 / 0.5) - coefficients(fit)[1]) / coefficients(fit)[2]
abline(h = 0.5, col = "red")

# d)
(log(0.99 / 0.01) - coefficients(fit)[1]) / coefficients(fit)[2]
abline(h = 0.99, col = "red")




# _____________________________________________________________________________
##### Classification and regression trees
##### --------------------------

### Classification Trees: Carseats dataset

help(Carseats)
dim(Carseats)
names(Carseats)

attach(Carseats)

hist(Sales)
abline(v = 8, lwd = 3, col = "red")

High <- ifelse(Sales <= 8, "No", "yes")

detach(Carseats)

Carseats <- data.frame(Carseats, High)

library(tree)
help(tree)

tree.carseats <- tree(factor(High) ~ . - Sales, Carseats)
summary(tree.carseats)

plot(tree.carseats)
text(tree.carseats, pretty = 0)

tree.carseats

# use cross validation to prune the tree optimally
help(cv.tree)

set.seed(1)
cv.carseats <- cv.tree(tree.carseats, FUN = prune.misclass)

names(cv.carseats)
cv.carseats

plot(
  cv.carseats$size,
  cv.carseats$dev,
  type = "b",
  xlab = "size",
  ylab = "misclass"
)

help(prune.misclass)
prune.carseats <- prune.misclass(tree.carseats, best = 12)

plot(prune.carseats)
text(prune.carseats, pretty = 0)


### Regression Trees: Boston housing dataset

help(Boston)
dim(Boston)
names(Boston)

tree.boston <- tree(medv ~ ., Boston)
summary(tree.boston)

plot(tree.boston)
text(tree.boston, pretty = 0)

cv.boston <- cv.tree(tree.boston)

plot(
  cv.boston$size,
  cv.boston$dev,
  type = "b",
  xlab = "size",
  ylab = "deviance"
)

prune.boston <- prune.tree(tree.boston, best = 4)

plot(prune.boston)
text(prune.boston, pretty = 0)




#### Pb3 of 09/02/2022  ####
# The file wine.txt reports the data on the alcohol content in 179 bottles of wine.
# For the alcohol content consider a linear model, accounting for the sugar content
# of grapes, and for type of wine ('Red', 'Rose', 'White'):
#   alcoholg = b0,g + b1,g*sugar + eps,
# with eps ??? N(0, sigma^2) and g the grouping structure induced by the type of wine.
# a) Estimate the parameters of the model ({b0,g, b1,g, sigma}). Verify the model assumptions,
#    reporting any plot you consider important.
# b) Perform two statistical tests - each at level 1% - to verify if
#     - there is a significant dependence of the mean alcohol content on the type of wine;
#     - there is a significant dependence of the mean alcohol content on the sugar content.
# c) Based on tests (b) or any other test deemed relevant, reduce the model and report
#    the updated model parameters.
# d) Build a prediction interval at 99% for a new bottle of red wine made with grapes
#    with 20 g of sugar.

data <- read.table("wine.txt")
head(data)
data$type <- as.factor(data$type)

# a)
D1 <- ifelse(data$type == "Red", 1, 0)
D2 <- ifelse(data$type == "White", 1, 0)
fit <-
  lm(alcohol ~ D1 + D2 + sugar + D1:sugar + D2:sugar, data = data)
summary(fit)

fit2 <- lm(alcohol ~ type + sugar + type:sugar, data = data)
summary(fit2)

beta0 <-
  c(coef(fit)[1], coef(fit)[1] + coef(fit)[2], coef(fit)[1] + coef(fit)[3])
beta1 <-
  c(coef(fit)[4], coef(fit)[4] + coef(fit)[5], coef(fit)[4] + coef(fit)[6])
beta0
beta1


shapiro.test(residuals(fit))

par(mfrow = c(2, 2))
plot(fit)

# b)
library(car)
linearHypothesis(
  fit,
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 1)
  ),
  c(0, 0, 0, 0)
)

linearHypothesis(
  fit,
  rbind(
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 1)
  ),
  c(0, 0, 0)
)

# c)
linearHypothesis(
  fit,
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0)
  ),
  c(0, 0)
)


fit2 <- lm(alcohol ~ sugar + D1:sugar + D2:sugar, data = data)
summary(fit2)

c(coef(fit2)[1])
c(
  coef(fit2)[2],
  coef(fit2)[2] + coef(fit2)[3],
  coef(fit2)[2] + coef(fit2)[4]
)

# d)
predict(fit2,
  data.frame(sugar = 20, D1 = 1, D2 = 0),
  interval = "prediction",
  level = 0.99
)
