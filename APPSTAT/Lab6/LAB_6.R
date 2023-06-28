### --------------------###
### LAB 6 (31/03/2023) ###
### --------------------###

### TOPICS:
### One-way ANOVA and MANOVA
### Two-ways ANOVA and MANOVA

library(MVN)
library(car)
library(heplots)

# _____________________________________________________________________________
##### One-way ANOVA
##### (p=1, g=6) -> one varible, six groups
##### ---------------
# Dataset Chicken Weights:
# Data from an experiment to measure and compare the effectiveness
# of various feed supplements on the growth rate of chickens.
# (71 observations)
rm(list = ls())

help(chickwts)

head(chickwts)
dim(chickwts)
summary(chickwts)

attach(chickwts)

# Group-wise boxplot
plot(
  feed,
  weight,
  xlab = "treat",
  ylab = "weight",
  col = "grey85",
  main = "Dataset Chicken Weights"
)

feed.levels <- levels(feed)
feed.levels
for (i in 1:length(feed.levels)) {
  m <- mean(weight[feed == feed.levels[i]])
  abline(h = m, col = rainbow(6)[i], lwd = 2)
}

dev.off()

### Model is: weigth.ij = mu + tau.i + eps.ij; eps.ij ~ N(0,sigma^2)
### where: mu = mean of the population
###        tau.i = effect of the i-th feed supplement
###        eps.ij = residual term specific to the i-th feed supplement and
###                      the j-th chicken
### RMK: eps.ij are independent and identically normally distributed (iid)
###   all with the same variance sigma^2 (homoscedasticity)
### Test:
### hyp. H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0
###         i.e., all the true means are equal
### hyp. H1: (H0)^c
### i.e.,
### H0: The feed supplements don't have effect
###     (= "chickens belong to a single population")
### H1: At least one feed supplement has effect
###     (= "chickens belong to 2, 3, 4, 5 or 6 populations")

# Plot showing the estimator for the mean in H0 and a particular case of H1
# where we assume that there are as much populations as groups

par(mfrow = c(1, 2))
barplot(
  rep(mean(weight), 6),
  names.arg = levels(feed),
  ylim = c(0, max(weight)),
  las = 2,
  col = "grey85",
  main = "Model under H0"
)
barplot(
  tapply(weight, feed, mean),
  names.arg = levels(feed),
  ylim = c(0, max(weight)),
  las = 2,
  col = rainbow(6),
  main = "Model under a special case of
        H1 with 6 populations"
)

dev.off()

# This is a case of one-way ANOVA: one variable (weight) observed
# over g=6 levels (feed)
n <- length(feed) # total number of obs.
ng <- table(feed) # number of obs. in each group
treat <- levels(feed) # levels of the treatment
g <- length(treat) # number of levels (i.e., of groups)

### verify the assumptions:
# 1) normality (univariate) in each group (6 tests)
Ps <- c(
  shapiro.test(weight[feed == treat[1]])$p,
  shapiro.test(weight[feed == treat[2]])$p,
  shapiro.test(weight[feed == treat[3]])$p,
  shapiro.test(weight[feed == treat[4]])$p,
  shapiro.test(weight[feed == treat[5]])$p,
  shapiro.test(weight[feed == treat[6]])$p
)
Ps
# We can assume normality for each group (p-values are above 10%)


# 2) same covariance structure (= same sigma^2)
Var <- c(
  var(weight[feed == treat[1]]),
  var(weight[feed == treat[2]]),
  var(weight[feed == treat[3]]),
  var(weight[feed == treat[4]]),
  var(weight[feed == treat[5]]),
  var(weight[feed == treat[6]])
)
Var

# test of homogeneity of variances for normal samples (Bartlett's test)
# H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6
# H1: there exist i,j s.t. sigma.i!=sigma.j
# WARNING: Test extremely sensitive to departures from normality (low robustness)
# Normality assumption mandatory
bartlett.test(weight, feed)

### One-way ANOVA
### -----------------
help(aov)

fit <- aov(weight ~ feed)

summary(fit)

### How to read the summary:
#              Df   Sum Sq      Mean Sq      F value     Pr(>F)
#  treat      (g-1) SStreat  SStreat/(g-1)  Fstatistic  p-value [H0: tau.i=0 for every i]
#  Residuals  (n-g) SSres     SSres/(n-g)
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
###
# Reject the test -> very small p-value -> there is evidence to state that
# at least one treatment (feed supplement) has an effect on the growth rate
# of chicken
# manual computation
alpha <- 0.001
SStreat <- sum(ng * (tapply(weight, feed, mean) - mean(weight))^2)
SStreat
# vector of the means of each group
m <- tapply(weight, feed, mean)
SSres <- sum((weight - m[feed])^2)
SSres

# F-statistic
Fstat <- (SStreat / (g - 1)) / (SSres / (n - g))
Fstat

# critical value for the F-test
cfr.fisher <- qf(1 - alpha, g - 1, n - g) # critical value for the F-test
Fstat < cfr.fisher
# -> we reject the null hypothesis

p.value <- 1 - pf(Fstat, g - 1, n - g) # p-value (equivalent to summary(fit)$"Pr(>F)"[1]

### We reject the test, i.e., we have evidence to state that the
### treatment (feed supplement) has an effect on the growth rate
### of chicken.
### Which supplement is responsible for this? To see this, we need to
### do g*(g-1)/2 comparisons.
### We use Bonferroni
k <- g * (g - 1) / 2
k
alpha <- 0.05

Mediag <- tapply(weight, feed, mean) # group-wise means
Mediag
SSres <- sum(residuals(fit)^2)
S <- SSres / (n - g)
S

# Example: CI for the difference "casein - horsebean"
paste(treat[1], "-", treat[2])
as.numeric(c(
  Mediag[1] - Mediag[2] - qt(1 - alpha / (2 * k), n - g) *
    sqrt(S * (1 / ng[1] + 1 / ng[2])),
  Mediag[1] - Mediag[2] + qt(1 - alpha / (2 * k), n - g) *
    sqrt(S * (1 / ng[1] + 1 / ng[2]))
))

# CI for all the differences
ICrange <- NULL
for (i in 1:(g - 1)) {
  for (j in (i + 1):g) {
    lower <- Mediag[i] - Mediag[j] - qt(1 - alpha / (2 * k), n - g) *
      sqrt(S * (1 / ng[i] + 1 / ng[j]))
    upper <- Mediag[i] - Mediag[j] + qt(1 - alpha / (2 * k), n - g) *
      sqrt(S * (1 / ng[i] + 1 / ng[j]))
    print(paste(treat[i], "-", treat[j]))
    print(as.numeric(c(lower, upper)))
    ICrange <- rbind(ICrange, as.numeric(c(lower, upper)))
  }
}

ICrange


par(mfrow = c(1, 2))
plot(
  feed,
  weight,
  xlab = "treat",
  ylab = "weight",
  col = rainbow(6),
  las = 2
)

h <- 1
plot(
  c(1, g * (g - 1) / 2),
  range(ICrange),
  pch = "",
  xlab = "pairs treat",
  ylab = "Conf. Int. tau weight"
)
for (i in 1:(g - 1)) {
  for (j in (i + 1):g) {
    ind <- (i - 1) * g - i * (i - 1) / 2 + (j - i)
    lines(c(h, h), c(ICrange[ind, 1], ICrange[ind, 2]), col = "grey55")

    points(h, Mediag[i] - Mediag[j], pch = 16, col = "grey55")

    points(h, ICrange[ind, 1], col = rainbow(6)[j], pch = 16)

    points(h, ICrange[ind, 2], col = rainbow(6)[i], pch = 16)

    h <- h + 1
  }
}
abline(h = 0)

# As we can see from the plot feeding chickes with casein and horsebean is
# significantly different (the CI does not contain 0)
# The same holds for all the pairs for which the CI does not contain 0,
# while if the CI contains 0 we can't say that the two treatments are
# significantly different

## RMK: even if the CI contains 0, we can't say that the two treatments
# are equal, but only that we can't reject the null hypothesis of
# equality. In fact looking for example at the pair horsebean - linseed,
# we can see that the CI contains 0, but chiens fed with horsebean are
# heavier than those fed with linseed: however don't have enough statistical
# evidence to reject the null hypothesis of equality.

### Be careful to apply transitivity in the statistical context!!
### {A not significantly (stat.) different from B} [A = B]   (red vs cyan)
###                       and
### {B not significantly (stat.) different from C} [B = C]   (cyan vs green)
###                  does NOT imply
### {A not significantly (stat.) different from C} [A = C]   (red vs green)
### Note. If we don't reject H0, we are not proving that A=B but
###       we are saying that we can't prove that A!=B

dev.off()

# Let's change the criterion to control the univariate rejection
# (multiple testing)

# We build k one-at-a-time confidence intervals, each of level alpha
# (not including the Bonferroni correction)
Auni <- matrix(0, 6, 6)
for (i in 1:6) {
  for (j in i:6) {
    Auni[i, j] <-
      Mediag[i] - Mediag[j] + qt(1 - alpha / 2, n - g) * sqrt(S * (1 / ng[i] + 1 /
        ng[j]))
  }
  for (j in 1:i) {
    Auni[i, j] <-
      Mediag[j] - Mediag[i] - qt(1 - alpha / 2, n - g) * sqrt(S * (1 / ng[i] + 1 /
        ng[j]))
  }
  Auni[i, i] <- 0
}


par(mfrow = c(1, 2))
h <- 1
plot(
  c(1, g * (g - 1) / 2),
  range(Auni),
  pch = "",
  xlab = "pairs treat",
  ylab = "CI delta weight",
  main = "Univariate Conf. Int.",
  col = "grey55"
)
for (i in 1:5) {
  for (j in (i + 1):6) {
    lines(c(h, h), c(Auni[i, j], Auni[j, i]))

    points(h, Mediag[i] - Mediag[j], pch = 16, col = "grey55")

    points(h, Auni[i, j], col = rainbow(6)[i], pch = 16)

    points(h, Auni[j, i], col = rainbow(6)[j], pch = 16)

    h <- h + 1
  }
}
abline(h = 0)

# Bonferroni is less restrictive than the one-at-a-time criterion

# We compute the p-values of the univariate tests

# Matrix of tests for the difference between all the pairs
P <- matrix(0, 6, 6)
for (i in 1:6) {
  for (j in i:6) {
    P[i, j] <-
      (1 - pt(abs((
        Mediag[i] - Mediag[j]
      ) / sqrt(
        S * (1 / ng[i] + 1 / ng[j])
      )), n - g)) * 2
  }
  for (j in 1:i) {
    P[i, j] <-
      (1 - pt(abs((
        Mediag[i] - Mediag[j]
      ) / sqrt(
        S * (1 / ng[i] + 1 / ng[j])
      )), n - g)) * 2
  }
  P[i, i] <- 0
}
P

# Vector of p-values
p <- c(P[1, 2:6], P[2, 3:6], P[3, 4:6], P[4, 5:6], P[5, 6])
p

plot(
  1:15,
  p,
  ylim = c(0, 1),
  type = "b",
  pch = 16,
  col = "grey55",
  xlab = "pairs treat",
  main = "P-values"
)
abline(h = alpha, lty = 2)

# Bonferroni correction
p.bonf <- p.adjust(p, "bonf")
lines(1:15,
  p.bonf,
  col = "blue",
  pch = 16,
  type = "b"
)

# Correction according to the false discovery rate (Benjamini-Hockberg)
p.fdr <- p.adjust(p, "fdr")
lines(1:15,
  p.fdr,
  col = "red",
  pch = 16,
  type = "b"
)

legend(
  "topleft",
  c("Not corr.", "Bonf.", "BH"),
  col = c("grey55", "blue", "red"),
  pch = 16
)

which(p.bonf < alpha)
which(p.fdr < alpha)

graphics.off()
detach(chickwts)

# _______________________________________________________________________________
##### One-way MANOVA
##### (p=4, g=3)
##### ---------------

library(MVN)
library(car)
library(heplots)

help(iris)
head(iris)
dim(iris)

### Variables: Sepal and Petal Length and Sepal and Petal Width of iris
### (p = 4)
### Groups: species (setosa, versicolor, virginica; g = 3)
### n1 = n2 = n3 = 50 (balanced design)

attach(iris)

species.name <-
  factor(Species, labels = c("setosa", "versicolor", "virginica"))
iris4 <- iris[, 1:4] # variables

### Data exploration
colore <- rep(rainbow(3), each = 50)

pairs(iris4, col = colore, pch = 16)

dev.off()

i1 <- which(species.name == "setosa")
i2 <- which(species.name == "versicolor")
i3 <- which(species.name == "virginica")

par(mfrow = c(1, 3))
boxplot(iris4[i1, ],
  main = "SETOSA",
  ylim = c(0, 8),
  col = rainbow(4)
)
boxplot(iris4[i2, ],
  main = "VERSICOLOR",
  ylim = c(0, 8),
  col = rainbow(4)
)
boxplot(iris4[i3, ],
  main = "VIRGINICA",
  ylim = c(0, 8),
  col = rainbow(4)
)


par(mfrow = c(1, 4))
boxplot(
  iris4[, 1] ~ species.name,
  main = "Sepal Length",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 2] ~ species.name,
  main = "Sepal Width",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 3] ~ species.name,
  main = "Petal Length",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 4] ~ species.name,
  main = "Petal Width",
  ylim = c(0, 8),
  col = rainbow(3)
)

graphics.off()

### Model: X.ij = mu + tau.i + eps.ij; eps.ij ~ N_p(0, Sigma), X.ij, mu, tau.i and eps.ij in R^4
### Test:
### H0: tau.1 = tau.2 = tau.3  = (0,0,0,0)'
### H1: (H0)^c
### that is
### H0: The membership to an iris species hasn't any significant effect on the mean
###     of X.ij (in any direction of R^4)
### H1: There exists at least one direction in R^4 along which at least two species
###     have some feature significantly different

# This is a case of one-way MANOVA: four variables (Sepal.Length, Sepal.Width,
# Petal.Length, Petal.Width) observed over g=3 levels (setosa, versicolor, virginica)

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n <- n1 + n2 + n3

g <- length(levels(species.name))
p <- 4


### Verify the assumptions:
# 1)  normality (multivariate) in each group (3 tests)
Ps <- NULL
for (i in 1:g) {
  mvn.test <- mvn(data = iris[get(paste("i", i, sep = "")), 1:4])
  Ps <- c(Ps, mvn.test$multivariateNormality$`p value`)
}
Ps

# I would reject at 5% the null hypothesis of multivariate normality for the
# first group (setosa) but not for the other two groups (versicolor and
# virginica)

# 2) same covariance structure (= same covariance matrix Sigma)
S <- cov(iris4)
S1 <- cov(iris4[i1, ])
S2 <- cov(iris4[i2, ])
S3 <- cov(iris4[i3, ])

# Qualitatively:
round(S1, digits = 1)
round(S2, digits = 1)
round(S3, digits = 1)

# Box's M test for homogeneity of covariance matrices
# Should be done only if n_i > ~20 for all i, p < 5 and g < 5
# WARNING: Very sensitive to departure from normality and too restrictive for ANOVA
# especially if we have a high number of samples
summary(boxM(iris4, species.name))
# p-value is very small -> reject the null hypothesis of homogeneity of covariance matrices
# very stron evidence that the covariance matrices are not the same
# -> anyway this won't be a problem for the MANOVA test
# Other way to compare: heatmap on the covariance matrices

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

dev.off()

# In any case, we decide to carry on the MANOVA

# Note: We can verify the assumptions a posteriori on the residuals of
#       the estimated model

### One-way MANOVA
### ----------------

help(manova)
help(summary.manova)

fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit, test = "Wilks")
# p-value is very small -> reject the null hypothesis of no effect of the
# membership to a species on the mean of the features
# Wilks lambda: ratio of the determinants of the within-group and the
#              total covariance matrices
# Exact tests for p<=2 or g<=3 already implemented in R

# Note: since g=3 the test is exact
#       (cfr. JW pag.300)
### Very small pvalue, hence we
### reject the test, i.e., we have statistical evidence to state that
### the factor "Species" has an effect on the mean features
### of the flowers.
### Who's the responsible for this?

### Via ANOVA: for each of the p=4 variables we perform an ANOVA test
###            to verify if the membership to a group has influence
###            on the mean of the variable (we explore separately the
###            4 axes directions in R^4)
summary.aov(fit)

# Each of the 4 variables is significantly influenced by the
# factor species.
# Note: this analysis does NOT say:
#       a) which group differ
#       b) with respect to which variables the groups in (a) differ
# => As for the ANOVA, we build confidence intervals (many more!)

### Via Bonferroni
alpha <- 0.05
k <- p * g * (g - 1) / 2
qT <- qt(1 - alpha / (2 * k), n - g)

summary.manova(fit, test = "Wilks")$SS

W <- summary.manova(fit)$SS$Residuals
m <- sapply(iris4, mean) # estimates mu
m1 <- sapply(iris4[i1, ], mean) # estimates mu.1 = mu + tau.1
m2 <- sapply(iris4[i2, ], mean) # estimates mu.2 = mu + tau.2
m3 <- sapply(iris4[i3, ], mean) # estimates mu.3 = mu + tau.3

inf12 <- m1 - m2 - qT * sqrt(diag(W) / (n - g) * (1 / n1 + 1 / n2))
sup12 <- m1 - m2 + qT * sqrt(diag(W) / (n - g) * (1 / n1 + 1 / n2))
inf13 <- m1 - m3 - qT * sqrt(diag(W) / (n - g) * (1 / n1 + 1 / n3))
sup13 <- m1 - m3 + qT * sqrt(diag(W) / (n - g) * (1 / n1 + 1 / n3))
inf23 <- m2 - m3 - qT * sqrt(diag(W) / (n - g) * (1 / n2 + 1 / n3))
sup23 <- m2 - m3 + qT * sqrt(diag(W) / (n - g) * (1 / n2 + 1 / n3))

CI <- list(
  setosa_versicolor    = cbind(inf12, sup12),
  setosa_virginica     = cbind(inf13, sup13),
  versicolor_virginica = cbind(inf23, sup23)
)
CI

# Now we have a complete frame (intervals for all the components of
# tau_1, tau_2 e tau_3): it is fault of all the groups and all the
# variables!


par(mfrow = c(2, 4))
boxplot(
  iris4[, 1] ~ species.name,
  main = "Sepal Length",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 2] ~ species.name,
  main = "Sepal Width",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 3] ~ species.name,
  main = "Petal Length",
  ylim = c(0, 8),
  col = rainbow(3)
)
boxplot(
  iris4[, 4] ~ species.name,
  main = "Petal Width",
  ylim = c(0, 8),
  col = rainbow(3)
)

mg <- rbind(m1, m2, m3)
sp.name <-
  c("Sepal Length", "Sepal Width", "Petal Length", "Petal Width")
for (k in 1:4) {
  plot(
    c(1, g * (g - 1) / 2),
    ylim = c(-4, 4),
    xlim = c(1, 3),
    pch = "",
    xlab = "pairs treat",
    ylab = paste("CI tau", k),
    main = paste("CI tau", sp.name[k])
  )
  lines(c(1, 1), c(CI[[1]][k, 1], CI[[1]][k, 2]))

  points(1, mg[1, k] - mg[2, k], pch = 16)

  points(1, CI[[1]][k, 1], col = rainbow(g)[2], pch = 16)

  points(1, CI[[1]][k, 2], col = rainbow(g)[1], pch = 16)

  lines(c(2, 2), c(CI[[2]][k, 1], CI[[2]][k, 2]))

  points(2, mg[1, k] - mg[3, k], pch = 16)

  points(2, CI[[2]][k, 1], col = rainbow(g)[3], pch = 16)

  points(2, CI[[2]][k, 2], col = rainbow(g)[1], pch = 16)

  lines(c(3, 3), c(CI[[3]][k, 1], CI[[3]][k, 2]))

  points(3, mg[2, k] - mg[3, k], pch = 16)

  points(3, CI[[3]][k, 1], col = rainbow(g)[3], pch = 16)

  points(3, CI[[3]][k, 2], col = rainbow(g)[2], pch = 16)

  abline(h = 0)
}

dev.off()

# For each variable, and for each pair of species, we have all confidence
# intervals that are not crossing the zero line.
# There is evidence that there is an effect for each of the specimen
# and for each of the variables.
# This means that based on these features, I should be able to
# distinguish the three species, i.e. using a classifier.

# _____________________________________________________________________________
##### Two-ways ANOVA
##### (p=1, g=2, b=2)
##### ----------------------

### Problem 4 of 14/09/06
### ------------------------
# In a small village in Switzerland there are two gas stations:
# one of Esso and one of Shell. Both sell either gasoline 95 octanes
# and 98 octanes.
# A young statistician wants to find out which is the best gas station
# and the best gasoline to refuel his car, in order to maximize the
# number of kilometers covered with a single refueling.
# After 8 refuelings, the measured performances are:
# km/l  : (18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
# distr.: ('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell')
# benz. : ('95','95','98','98','95','95','98','98')
# (a) Via a two-ways ANOVA identify which is the best station and the
#     best gasoline for the young statistician to refuel his car.
# (b) Is there an interaction between the gas station and the gasoline?

### Variables: distance covered [km/l]
### factor1:   Gas station   (0=Esso, 1=Shell)
### factor2:   Gasoline      (0=95,   1=98)
### Balanced design

rm(list = ls())

km <- c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
distr <-
  factor(c(
    "Esso",
    "Esso",
    "Esso",
    "Esso",
    "Shell",
    "Shell",
    "Shell",
    "Shell"
  ))
benz <- factor(c("95", "95", "98", "98", "95", "95", "98", "98"))
distr_benz <-
  factor(c(
    "Esso95",
    "Esso95",
    "Esso98",
    "Esso98",
    "Shell95",
    "Shell95",
    "Shell98",
    "Shell98"
  ))

# g = 2, number of factor 1 levels
g <- length(levels(distr))
# b = 2, number of factor 2 levels
b <- length(levels(benz))
# n = 2, group sizes (2 observations for each combination of the two factors)
n <- length(km) / (g * b)

# overall mean
M <- mean(km)
M
# Mean per factor 1 level (gas station)
Mdistr <- tapply(km, distr, mean)
Mdistr
# Mean per factor 2 level (type of gasoline)
Mbenz <- tapply(km, benz, mean)
Mbenz
# Mean per factor 1 and factor 2 level
Mdistr_benz <- tapply(km, distr_benz, mean)
Mdistr_benz

par(mfrow = c(2, 3), las = 2)
barplot(
  rep(M, 4),
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  main = "No factor"
)
barplot(
  rep(Mdistr, each = 2),
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  col = rep(c("blue", "red"), each = 2),
  main = "Only Fact. Stat."
)
barplot(
  rep(Mbenz, times = 2),
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  col = rep(c("darkgreen", "orange"), times = 2),
  main = "Only Fact. Gas"
)
barplot(
  c(
    Mdistr[1] + Mbenz[1] - M,
    Mdistr[1] + Mbenz[2] - M,
    Mdistr[2] + Mbenz[1] - M,
    Mdistr[2] + Mbenz[2] - M
  ),
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  col = rep(c("darkgreen", "orange"), times = 2),
  density = rep(10, 4),
  angle = 135,
  main = "Additive model Stat.+Gas"
)
barplot(
  c(
    Mdistr[1] + Mbenz[1] - M,
    Mdistr[1] + Mbenz[2] - M,
    Mdistr[2] + Mbenz[1] - M,
    Mdistr[2] + Mbenz[2] - M
  ),
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  col = rep(c("blue", "red"), each = 2),
  density = rep(10, 4),
  add = TRUE
)
barplot(
  Mdistr_benz,
  names.arg = levels(distr_benz),
  ylim = c(0, 24),
  col = rainbow(5)[2:5],
  main = "Model with Interact. Stat.+Gas."
)
plot(
  distr_benz,
  km,
  col = rainbow(5)[2:5],
  ylim = c(0, 24),
  xlab = ""
)

### RMK: all these are just estimates! we have to do some statistical tests
### to quantify the goodness of these estimates

### Two-ways ANOVA
### ----------------

### Model with interaction (complete model):
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk ---
###     with eps.ijk~N(0,sigma^2),
###     i=1,2 (effect station), j=1,2 (effect gasoline)

fit.aov2.int <- aov(km ~ distr + benz + distr:benz)
summary.aov(fit.aov2.int)

### Test:
### 1) H0: gamma.11 = gamma.12 = gamma.21 = gamma.22 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: There is no significant interaction between the factors station
###        and gasoline in terms of performances
###    H1: There exists a significant interaction between the factors station
###        and gasoline in terms of performances
###
### 2) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gas station" doesn't significantly influence performances
###    H1: The effect "gas station" significantly influences performances
###
### 3) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
###    i.e.,
###    H0: The effect "gasoline" doesn't significantly influence performances
###    H1: The effect "gasoline" significantly influences performances

# Test 1): Let's focus on the row of the summary distr:benz :
#             Df Sum Sq Mean Sq F value  Pr(>F)
# distr:benz   1  10.35   10.35   6.884 0.05857 .
# The P-value of test 1) is 0.05857. Reject at 10%, don't reject at 1%,5% -> ?

# Test 2): Let's focus on the row of the summary distr:
#             Df Sum Sq Mean Sq F value  Pr(>F)
# distr        1   1.53    1.53   1.018 0.37001
# The P-value of test 2) is 0.37001. Don't reject at 10%, 5%, 1% -> not significant

# Test 3): Let's focus on the row of the summary benz:
#             Df Sum Sq Mean Sq F value  Pr(>F)
# benz         1  66.70   66.70  44.357 0.00264 **
# The P-value of test 3) is 0.00264. Reject at 10%, 5%, 1% -> significant

# From test 1): We don't have strong evidence that the interaction has effect
# => remove the interaction term and estimate the model without interaction

### Additive model:
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2),
###     i=1,2 (effect station), j=1,2 (effect gasoline)
fit.aov2.ad <- aov(km ~ distr + benz)
summary.aov(fit.aov2.ad)
# Remark: by removing the interaction, the residual degrees of freedom increase!

# Test: 2bis) H0: tau.1 = tau.2 = 0    vs   H1: (H0)^c
# From the summary:
#             Df Sum Sq Mean Sq F value  Pr(>F)
# distr        1   1.53    1.53   0.468 0.52440
# The P-value of test 2bis) is 0.52440. Don't reject at 10%, 5%, 1% -> not significant

# Test: 3bis) H0: beta.1 = beta.2 = 0    vs   H1: (H0)^c
# From the summary:
#             Df Sum Sq Mean Sq F value  Pr(>F)
# benz         1  66.70   66.70  20.378 0.00632 **
# The P-value of test 2bis) is 0.00632. Reject at 10%, 5%, 1% -> significant

### Note: These aren't the only tests we can do!
### Example: global test for the significance of the two treatments
###          (model without interaction)
SSdistr <-
  sum(n * b * (Mdistr - M)^2) # or from the summary: 1.53
SSbenz <-
  sum(n * g * (Mbenz - M)^2) # or from the summary: 66.70
SSres <-
  sum((km - M)^2) - (SSdistr + SSbenz) # or from the summary: 16.37

Ftot <-
  ((SSdistr + SSbenz) / ((g - 1) + (b - 1))) / (SSres / (n * g * b - g -
    b + 1))
# (Variability between treatments / dofs) /
#                               (Variability within treatments / dofs)
Ptot <-
  1 - pf(Ftot, (g - 1) + (b - 1), n * g * b - g - b + 1) # attention to the dof!
Ptot

# Test 2bis): there is no evidence that the factor "gas station" has
#             effect on the performances (don't reject at any reasonable
#             level [high p-value!])
# => we remove the variable "station" and reduce to a one-way ANOVA

### Reduced additive model (ANOVA one-way, b=2):
### X.jk = mu + beta.j + eps.jk; eps.jk~N(0,sigma^2), only leave the gasoline
###     j=1,2 (effect gasoline)
fit.aov1 <- aov(km ~ benz)
summary.aov(fit.aov1)

SSres <- sum(residuals(fit.aov1)^2)

### Interval at 95% for the differences (reduced additive model)
### [b=2, thus one interval only]

t.stat <- qt(0.95, (n * g - 1) * b)
IC <-
  c(
    diff(Mbenz) - t.stat *
      sqrt(SSres / ((n * g - 1) * b) * (1 / (n * g) + 1 / (n * g))),
    diff(Mbenz) + t.stat *
      sqrt(SSres / ((n * g - 1) * b) * (1 / (n * g) + 1 / (n * g)))
  )
names(IC) <- c("Inf", "Sup")
IC # IC for mu(98)-mu(95)

### Note: we should have verified the hypotheses of normality and variance
###       homogeneity for the complete model, but with only 2 data for each
###       group we can't perform the tests.
### => we verify the assumptions on the reduced model (one-way ANOVA)
# 1) normality (univariate) in each group (2 tests)
Ps <- c(
  shapiro.test(km[benz == levels(benz)[1]])$p,
  shapiro.test(km[benz == levels(benz)[2]])$p
)
Ps

# 2) homogeneity of variances
bartlett.test(km, benz)

graphics.off()

# _____________________________________________________________________________
##### Two-ways MANOVA
##### (p=3, g=2, b=2)
##### ----------------------

### Example 6.13 JW (p.319)
### ------------------------
# The optimum conditions for extruding plastic films have been
# examined using a technique called Evolutionary Operation. In
# the course of the study, three responses were measured:
# X.1 = Tear Resistance
# X.2 = Gloss ("Shininess")
# X.3 = Opacity ("Transparency")
# at two levels of the factors:
# factor.1 = Rate of Extrusion     (0=Low level, 1=High level)
# factor.2 = Amount of an Additive (0=Low level, 1=High level)
# The measurements were repeated n=5 times at each combination
# of the factor levels.
rm(list = ls())

setwd("~/HPC/APPSTAT/Lab6")
library(MVN)
plastic <-
  read.table("T6-4.dat", col.names = c("Ex", "Ad", "Tr", "Gl", "Op"))
plastic

Ex <- factor(plastic$Ex, labels = c("L", "H")) # Treat.1
Ad <- factor(plastic$Ad, labels = c("L", "H")) # Treat.2
Ex
Ad

ExAd <- Ex
levels(ExAd) <- c("LL", "LH", "HL", "HH")
ExAd[Ex == "L" & Ad == "L"] <- "LL"
ExAd[Ex == "L" & Ad == "H"] <- "LH"
ExAd[Ex == "H" & Ad == "L"] <- "HL"
ExAd[Ex == "H" & Ad == "H"] <- "HH"

plastic3 <- plastic[, 3:5]

### Graphical exploration of the data
# effect of the treatments + their interaction on the first variable
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE))
boxplot(plastic3[, 1] ~ ExAd,
  main = "Model with Interac. Extrusion+Additive (Tear Resistance)",
  ylab = "Tr",
  col = "grey95"
)
boxplot(
  plastic3[, 1] ~ Ex,
  main = "Only Factor Extrusion",
  ylab = "Tr",
  col = c("red", "blue")
)
boxplot(
  plastic3[, 1] ~ Ad,
  main = "Only Factor Additive",
  ylab = "Tr",
  col = c("forestgreen", "gold")
)

# effect of the treatments + their interaction on the second variable
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE))
boxplot(plastic3[, 2] ~ ExAd,
  main = "Model with Interac. Extrusion+Additive (Gloss)",
  ylab = "Gl",
  col = "grey95"
)
boxplot(
  plastic3[, 2] ~ Ex,
  main = "Only Factor Extrusion",
  ylab = "Gl",
  col = c("red", "blue")
)
boxplot(
  plastic3[, 2] ~ Ad,
  main = "Only Factor Additive",
  ylab = "Gl",
  col = c("forestgreen", "gold")
)

# effect of the treatments + their interaction on the third variable
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE))
boxplot(plastic3[, 3] ~ ExAd,
  main = "Model with Interac. Extrusion+Additive (Opacity)",
  ylab = "Op",
  col = "grey95"
)
boxplot(
  plastic3[, 3] ~ Ex,
  main = "Only Factor Extrusion",
  ylab = "Op",
  col = c("red", "blue")
)
boxplot(
  plastic3[, 3] ~ Ad,
  main = "Only Factor Additive",
  ylab = "Op",
  col = c("forestgreen", "gold")
)

dev.off()

### Model with interaction (complete model):
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect Additive),
###     X.ijs, mu, tau.i, beta.j, gamma.ij in R^3

### Verify the assumptions (although we only have 5 data in each group!)
# 1) normality (multivariate) in each group (4 test)
Ps <-
  c(
    mvn(plastic3[ExAd == levels(ExAd)[1], ])$multivariateNormality$`p value`,
    mvn(plastic3[ExAd == levels(ExAd)[2], ])$multivariateNormality$`p value`,
    mvn(plastic3[ExAd == levels(ExAd)[3], ])$multivariateNormality$`p value`,
    mvn(plastic3[ExAd == levels(ExAd)[4], ])$multivariateNormality$`p value`
  )
Ps
# Can't do this, too few data

# 2) homogeneity of the covariance (qualitatively)
S1 <- cov(plastic3[ExAd == levels(ExAd)[1], ])
S2 <- cov(plastic3[ExAd == levels(ExAd)[2], ])
S3 <- cov(plastic3[ExAd == levels(ExAd)[3], ])
S4 <- cov(plastic3[ExAd == levels(ExAd)[4], ])



par(mfrow = c(1, 4))
image(
  S1,
  col = heat.colors(100),
  main = "Cov. S1",
  asp = 1,
  axes = FALSE,
  breaks = quantile(rbind(S1, S2, S3, S4), (0:100) / 100, na.rm = TRUE)
)
image(
  S2,
  col = heat.colors(100),
  main = "Cov. S2",
  asp = 1,
  axes = FALSE,
  breaks = quantile(rbind(S1, S2, S3, S4), (0:100) / 100, na.rm = TRUE)
)
image(
  S3,
  col = heat.colors(100),
  main = "Cov. S3",
  asp = 1,
  axes = FALSE,
  breaks = quantile(rbind(S1, S2, S3, S4), (0:100) / 100, na.rm = TRUE)
)
image(
  S4,
  col = heat.colors(100),
  main = "Cov. S4",
  asp = 1,
  axes = FALSE,
  breaks = quantile(rbind(S1, S2, S3, S4), (0:100) / 100, na.rm = TRUE)
)

dev.off()

### Two-ways MANOVA
### -----------------
### Model with interaction (complete model):
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect Additive),
###     X.ijs, mu, tau.i, beta.j, gamma.ij in R^3
fit <- manova(as.matrix(plastic3) ~ Ex + Ad + Ex:Ad)
summary.manova(fit, test = "Wilks")

### Model without interaction (additive model):
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N_p(0,Sigma), [p=3]
###     i=1,2 (effect Extrusion), j=1,2 (effect additive),
###     X.ijs, mu, tau.i, beta.j, in R^3
fit2 <- manova(as.matrix(plastic3) ~ Ex + Ad)
summary.manova(fit2, test = "Wilks")

# Both the treatments have a significant effect on the mean (but not
# their interaction, that we could remove)

# Let's verify if this is true for all the variables through appropriate
# conf. int. and tests on the components:

# ANOVA on the components (we look at the 3 axes-directions in R^3
#                          separately)
summary.aov(fit2)

# I reject the null hypothesis if the p-value is < alpha

# TEAR RESISTENCE
# The effect of extrusion has very significant evidence that the mean
# of the first variable (tear resistence) changes wrt the level of the
# treatment.
# Also the effect of additive has a significant evidence (we would reject
# the null hypothesis at a confidence level of 90% and 95%, but no at 99%)

# GLOSS
# The effect of extrusion has significant evidence that gloss changes wrt the
# level of the treatment (We would reject the null hypothesis at a confidence
# level of 90% and 95% (but no at 99%))
# The effect of additive has some evidence that gloss changes, but not
# that much (we would reject the null hypothesis at a confidence level of 90%,
# but not at 95% and 99%)

# OPACITY
# Both the effects of extrusion and additive do not have significant evidence
# that opacity changes wrt the level of the treatment (we would not reject
# the null hypothesis that the mean stays the same at any significant
# confidence level)


# Bonferroni
alpha <- 0.05
g <- 2
b <- 2
p <- 3
n <- 5
N <- n * g * b # 20

W <- summary.manova(fit2)$SS$Residuals

# how many comparisons?
k <- g * (g - 1) / 2 * p + b * (b - 1) / 2 * p
# because we have: g levels on the first treatment on p components
#                  b levels on the second treatment on p components
k

qT <- qt(1 - alpha / (2 * k), g * b * n - g - b + 1)
# the degrees of freedom of the residuals on the additive model are
# g*b*n-g-b+1

mExL <- sapply(plastic3[Ex == "L", ], mean)
mExH <- sapply(plastic3[Ex == "H", ], mean)
infEx <-
  mExH - mExL - qT *
    sqrt(diag(W) / (g * b * n - g - b + 1) * (2 / N + 2 / N))
supEx <-
  mExH - mExL + qT *
    sqrt(diag(W) / (g * b * n - g - b + 1) * (2 / N + 2 / N))

mAdL <- sapply(plastic3[Ad == "L", ], mean)
mAdH <- sapply(plastic3[Ad == "H", ], mean)
infAd <-
  mAdH - mAdL - qT *
    sqrt(diag(W) / (g * b * n - g - b + 1) * (2 / N + 2 / N))
supAd <-
  mAdH - mAdL + qT *
    sqrt(diag(W) / (g * b * n - g - b + 1) * (2 / N + 2 / N))

IC2 <-
  list(
    ExH_ExL = cbind(infEx, supEx),
    AdH_AdL = cbind(infAd, supAd)
  )
IC2

# With the bonferroni confidence intervals I would reject only the null
# hypothesis about the effect of extrusion on the tear resistence, while I
# would accept all the others.
# That's because computing the anova on each component separately is like
# computing the one at a time confidence intervals, while with Bonferroni
# we are computing the simultaneous confidence intervals.

par(mfrow = c(3, 4))
boxplot(
  plastic3[, 1] ~ Ex,
  main = "Fact.: Extrusion (Tear Resistance)",
  ylab = "Tr",
  col = rainbow(2 * 6)[c(1, 2)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[1]][1, ]),
  pch = "",
  main = "IC (tau.1-tau.2)[1]",
  xlab = "pairs treat",
  ylab = "IC (tau.1-tau.2)[1]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[1]][1, 1], IC2[[1]][1, 2]), col = "grey55")

points(1, (mExH - mExL)[1], pch = 16, col = "grey55")

points(1, IC2[[1]][1, 1], col = rainbow(2 * 6)[1], pch = 16)

points(1, IC2[[1]][1, 2], col = rainbow(2 * 6)[2], pch = 16)

abline(h = 0)

boxplot(
  plastic3[, 1] ~ Ad,
  main = "Fact.: Additive (Tear Resistance)",
  ylab = "Tr",
  col = rainbow(2 * 6)[c(7, 8)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[2]][1, ]),
  pch = "",
  main = "IC (beta.1-beta.2)[1]",
  xlab = "pairs treat",
  ylab = "IC (beta.1-beta.2)[1]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[2]][1, 1], IC2[[2]][1, 2]), col = "grey55")

points(1, (mAdH - mAdL)[1], pch = 16, col = "grey55")

points(1, IC2[[2]][1, 1], col = rainbow(2 * 6)[7], pch = 16)

points(1, IC2[[2]][1, 2], col = rainbow(2 * 6)[8], pch = 16)

abline(h = 0)

boxplot(
  plastic3[, 2] ~ Ex,
  main = "Fact.: Extrusion (Gloss)",
  ylab = "Gl",
  col = rainbow(2 * 6)[c(3, 4)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[1]][2, ]),
  pch = "",
  main = "IC (tau.1-tau.2)[2]",
  xlab = "pairs treat",
  ylab = "IC (tau.1-tau.2)[2]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[1]][2, 1], IC2[[1]][2, 2]), col = "grey55")

points(1, (mExH - mExL)[2], pch = 16, col = "grey55")

points(1, IC2[[1]][2, 1], col = rainbow(2 * 6)[3], pch = 16)

points(1, IC2[[1]][2, 2], col = rainbow(2 * 6)[4], pch = 16)

abline(h = 0)

boxplot(
  plastic3[, 2] ~ Ex,
  main = "Fact.: Additive (Gloss)",
  ylab = "Gl",
  col = rainbow(2 * 6)[c(9, 10)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[2]][2, ]),
  pch = "",
  main = "IC (beta.1-beta.2)[2]",
  xlab = "pairs treat",
  ylab = "IC (beta.1-beta.2)[2]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[2]][2, 1], IC2[[2]][2, 2]), col = "grey55")

points(1, (mAdH - mAdL)[2], pch = 16, col = "grey55")

points(1, IC2[[2]][2, 1], col = rainbow(2 * 6)[9], pch = 16)

points(1, IC2[[2]][2, 2], col = rainbow(2 * 6)[10], pch = 16)

abline(h = 0)

boxplot(
  plastic3[, 3] ~ Ex,
  main = "Fact.: Extrusion (Opacity)",
  ylab = "Op",
  col = rainbow(2 * 6)[c(5, 6)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[1]][3, ]),
  pch = "",
  main = "IC (tau.1-tau.2)[3]",
  xlab = "pairs treat",
  ylab = "IC (tau.1-tau.2)[3]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[1]][3, 1], IC2[[1]][3, 2]), col = "grey55")

points(1, (mExH - mExL)[3], pch = 16, col = "grey55")

points(1, IC2[[1]][3, 1], col = rainbow(2 * 6)[5], pch = 16)

points(1, IC2[[1]][3, 2], col = rainbow(2 * 6)[6], pch = 16)

abline(h = 0)

boxplot(
  plastic3[, 3] ~ Ex,
  main = "Only Factor Additive (Opacity)",
  ylab = "Op",
  col = rainbow(2 * 6)[c(11, 12)],
  ylim = c(-2, 10)
)
plot(
  c(1, g * (g - 1) / 2),
  range(IC2[[2]][3, ]),
  pch = "",
  main = "IC (beta.1-beta.2)[3]",
  xlab = "pairs treat",
  ylab = "IC beta.1[3]",
  ylim = c(-2, 10)
)
lines(c(1, 1), c(IC2[[2]][3, 1], IC2[[2]][3, 2]), col = "grey55")

points(1, (mAdH - mAdL)[3], pch = 16, col = "grey55")

points(1, IC2[[2]][3, 1], col = rainbow(2 * 6)[11], pch = 16)

points(1, IC2[[2]][3, 2], col = rainbow(2 * 6)[12], pch = 16)

abline(h = 0)

dev.off()

# _______________________________________________________________________________
##### Problem 3 of 18/02/09
##### ------------------------
# During the austral summer three species of penguins nest in the Antarctic
# Peninsula: Chinstrap, Adelie and Gentoo.
# Some biologists of the Artowski basis measured the weight [kg] of 90 adults:
# 15 males and 15 females for each of the three species (penguins.txt file).
# a) By using an ANOVA model with two factors, claim if gender and / or species
#    belonging significantly affect the weight.
# b) Using an appropriate model (possibly reduced), provide estimates (global
#    90% confidence) of means and variances of the groups identified at point
#    (a).
rm(list = ls())

penguins <- read.table("penguins.txt", header = TRUE)
head(penguins)

attach(penguins)

### question a)
# Model with interaction (complete model):
# X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2),
#     i=1,2 (effect gender), j=1,2,3 (effect species)
fit1 <- aov(weight ~ gender + species + species:gender)
summary(fit1)
# Very high p-value for interaction term, so we can remove it from the model

# Model without interaction (additive model):
# X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2),
#     i=1,2 (effect gender), j=1,2,3 (effect species)
fit2 <- aov(weight ~ gender + species, penguins)
summary(fit2)
# Very high p-value also for the gender term, so we can remove it from the
# model

# Reduced model (one-way):
# X.jk = mu + beta.j + eps.jk; eps.jk~N(0,sigma^2),
#     j=1,2,3 (effect species)
fit3 <- aov(weight ~ species, penguins)
summary(fit3)
# Very small p-value for species term, so we can reject the null hypothesis
# and conclude that species significantly affect the weight

# verify assumptions (on the last model):
# 1) normality
Ps <- c(
  shapiro.test(weight[species == levels(factor(species))[1]])$p,
  shapiro.test(weight[species == levels(factor(species))[2]])$p,
  shapiro.test(weight[species == levels(factor(species))[3]])$p
)
Ps
# I can assume normality on the three groups

# 2) homogeneity of variances
bartlett.test(weight, species)
# I can assume homogeneity of variances (even though I would reject at 1%,
# but not at 5%)


### question b)
n <- 15 # number of observations for each group
g <- 3 # number of species
b <- 2 # number of genders
DF <- fit3$df # degrees of freedom n*g*b - 1 - (g-1) = 15*3*2-1-2 = 90-3 = 87
Spooled <- sum(fit3$res^2) / DF
Spooled

means <- as.vector(tapply(weight, species, mean))
names(means) <- levels(factor(species))
means

alpha <- 0.10
k <-
  4 # g + 1 = 4 (g Conf Int for the means and 1 for the variance)
BF <-
  rbind(
    cbind(
      means - sqrt(Spooled / 30) * qt(1 - alpha / (2 * k), DF),
      means + sqrt(Spooled / 30) * qt(1 - alpha / (2 * k), DF)
    ),
    c(
      Spooled * DF / qchisq(1 - alpha / (2 * k), DF),
      Spooled * DF / qchisq(alpha / (2 * k), DF)
    )
  )
rownames(BF)[4] <- "Var."
BF

detach(penguins)
