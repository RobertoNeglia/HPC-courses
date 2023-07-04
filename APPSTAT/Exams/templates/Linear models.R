X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2017/2017-07-03/garden.txt")
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-06-16/Exercise 3/danceability.txt", header=TRUE)
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-01-19/tattoo.txt")
# X <- read.table("~/GitHub/Applied-Statistics-Exam//Exams of previous years/2022/2022-02-09/wine.txt")
n <- dim(X)[1]

# a) Linear model
# categ.factor <- ifelse(X$type=="Yes",1,0) # categorical factor
type.f <- factor(X$type)
lmod <- lm(alcohol ~ -1 + type.f + type.f:sugar, data = X)
summary(lmod)
Betas <- lmod$coefficients
Betas
sigma <- summary(lmod)$sigma
sigma
# Assumptions (gaussianity of residuals, homosckedasticity, independence from regressors (lack of patterns))
x11()
par(mfrow = c(2, 2))
plot(lmod)
shapiro.test(lmod$residuals)$p
# qqnorm(lmod$residuals) # already present in plot(lmod)
# qqline(lmod$residuals)
x11()
plot(X$sugar, scale(lmod$residuals, scale = T, center = F), main = "Residuals vs regressor") # repeat plot for other regressors

# b) Hypothesis testing on the coefficients
library(car)
linearHypothesis(lmod, rbind(c(1, 0, 0, 0, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0)), c(0, 0, 0))

# c) Model reduction
library(regclass)
vif(lmod) # variance inflation factor: if big (above 5) risk of collinearity
lmod1 <- lm(extension ~ carps + maple + stones, data = X) # remove one at the time starting from highest p-value
summary(lmod1)
anova(lmod, lmod1) # test if we lose in explainability, we do not want to reject (high p-value = good)
# remember to check assumptions again

# d) Prediction and confidence intervals
x.new <- data.frame(Va = 35, Vi = 25, wind.factor = 1)
k <- 1 # num intervals, Bonferroni correction
alpha <- 0.01 / k
predict(lmod1, x.new, interval = "prediction", level = 1 - alpha) # interval=prediction/confidence

# e) Confidence intervals for the maximum of the means
max.mean.pos <- which.max(lmod$fitted.values)
x0 <- data.frame(pos.dummy = 1, season.dummy = 0, t = 2) # fill with the values of the regressors in X[max.mean.pos,]
alpha <- 0.05
predict(lmod, x0, interval = "confidence", level = 1 - alpha)

# f) Random intercept and PVRE
library(nlmeU)
library(corrplot)
library(nlme)
library(lattice)
library(plot.matrix)
library(lme4)
library(insight)
lmoder <- lmer(danceability ~ loudness + tempo + (1 | genre), data = X)
summary(lmoder)
confint(lmoder, oldNames = TRUE)
sigma2_eps <- as.numeric(get_variance_residual(lmoder))
sigma2_b <- as.numeric(get_variance_random(lmoder))
PVRE <- sigma2_b / (sigma2_b + sigma2_eps)
PVRE
x11()
dotplot(ranef(fm16.1mer, condVar = T)) # visualization of the random intercepts with their 95% confidence intervals

# g) Ridge regression
library(car)
lambda.c <- 10^seq(5, -3, length = 100)
fit.ridge <- lm.ridge(alcohol ~ -1 + type.f + type.f:sugar, data = X, lambda = lambda.c)
select(fit.ridge)

# h) Lasso regression
library(glmnet)
x <- model.matrix(alcohol ~ -1 + type.f + type.f:sugar, data = X)[, -1] # matrix of predictors
y <- X$alcohol # vector of response
lambda.grid <- 10^seq(5, -3, length = 100) # grid of candidate lambda's for the estimate
fit.lasso <- glmnet(x, y, lambda = lambda.grid) # default: alpha=1, if alpha=0 -> ridge regression
plot(fit.lasso, xvar = "lambda", label = TRUE, col = rainbow(dim(x)[2]))
legend("topright", dimnames(x)[[2]], col = rainbow(dim(x)[2]), lty = 1, cex = 1)
cv.lasso <- cv.glmnet(x, y, lambda = lambda.grid) # set lambda via cross validation, default: 10-fold CV
bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso
plot(cv.lasso)
abline(v = log(bestlam.lasso), lty = 1)
coef.lasso <- as.matrix(predict(fit.lasso, s = bestlam.lasso, type = "coefficients")) # coefficients for the optimal lambda
coef.lasso
