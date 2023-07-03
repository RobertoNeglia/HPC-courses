##############################################################
#############  Applied Statistics 2022/2023  #################
##############  Linear Mixed-effects models  #################
##############################################################

install.packages(c("plot.matrix", "insight"))

library(nlmeU)
library(corrplot)
library(nlme)
library(lattice)
library(plot.matrix)
library(lme4)
library(insight)

rm(list = ls())
graphics.off()

# Topics:
#   LINEAR MIXED MODELS WITH HOMOSCEDASTIC RESIDUALS
#   1. Linear Models with random intercept (q=0)
#   2. Linear Models with random intercept + slope (q=1)
#      2.1 general structure of D
#      2.2 diagonal D
#   3. Interpretation of random effects and PVRE
#   4. Prediction
#   5. Diagnostic
#   6. Models comparison
#
#   SUPPLEMENTARY MATERIAL: LINEAR MIXED MODELS WITH HETEROSCEDASTIC RESIDUALS (VarPower())
#   1. Linear Models with random intercept (q=0)
#   2. Linear Models with random intercept + slope (q=1)
#      2.1 general structure of D
#      2.2 diagonal D



data(armd) # Age-Related Macular Degeneration


## Linear Mixed-Effects models (LMM)
## In the LMM approach, the hierarchical structure of the data is directly addressed,
## with random effects that describe the contribution of the variability at different levels
## of the hierarchy to the total variability of the observations.

## Two main R packages:
## 1. 'lme4' --> it does not handle heteroscedastic residuals but it has a lot of "accessories"
## 2. 'nlme' --> it handles heteroscedastic residuals but it has less "accessories"

## We will use lmer() function in lme4 package for LMM models with homogeneous residuals and
## lme() function in nlme package for LMM models with heteroscedastic residuals

##########################################################
#### LINEAR MIXED MODELS WITH HOMOSCEDASTIC RESIDUALS ####
############ lme4 package --> lmer() function ############
##########################################################

######################################################################################################
# Model 1. Random intercept, homoscedastic residuals

# We now treat time as a numeric variable

fm16.1mer <- lmer(visual ~ visual0 + time * treat.f + (1 | subject),
     data = armd
)

summary(fm16.1mer)
confint(fm16.1mer, oldNames = TRUE)

## Var-Cov matrix of fixed-effects
vcovb <- vcov(fm16.1mer)
vcovb
corb <- cov2cor(vcovb)
nms <- abbreviate(names(fixef(fm16.1mer)), 5)
rownames(corb) <- nms
corb

## Var-Cov matrix of random-effects and errors
print(vc <- VarCorr(fm16.1mer), comp = c("Variance", "Std.Dev."))

sigma2_eps <- as.numeric(get_variance_residual(fm16.1mer))
sigma2_eps
sigma2_b <- as.numeric(get_variance_random(fm16.1mer))
sigma2_b


## Let's compute the conditional and marginal var-cov matrix of Y
sgma <- summary(fm16.1mer)$sigma

A <-
     getME(fm16.1mer, "A") # A  --> N x n, A represents the D (not italic)
I.n <- Diagonal(ncol(A)) # IN  --> n x n

## the conditional variance-covariance matrix of Y (diagonal matrix)
SigmaErr <- sgma^2 * (I.n)
SigmaErr[3:6, 3:6] ## visualization of individual 2
# Conditioned to the random effects b_i, we observe the var-cov of the errors
# that are independent and homoscedastic

## we visualize the first 20 rows/columns of the matrix
plot(as.matrix(SigmaErr[1:20, 1:20]), main = "Conditional estimated Var-Cov matrix of Y")

## the marginal variance-covariance matrix of Y (block-diagonal matrix)
V <-
     sgma^2 * (I.n + crossprod(A)) # V = s^2*(I_N+A*A) --> s^2*(I_N) is the error part, s^2*(A*A) is the random effect part
V[3:6, 3:6] #-> V is a block-diagional matrix, the marginal var-cov matrix

# visualization of the first 20 rows/columns
plot(as.matrix(V[1:20, 1:20]), main = "Marginal estimated Var-Cov matrix of Y")


# Another way to interpret the variance output is to note percentage of the subject variance out
# of the total, i.e. the Percentage of Variance explained by the Random Effect (PVRE).
# This is also called the intraclass correlation (ICC), because it is also an estimate of the within
# cluster correlation.

PVRE <- sigma2_b / (sigma2_b + sigma2_eps)
PVRE # 51% is very high!

## visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i for i=1,...,234
dotplot(ranef(fm16.1mer, condVar = T))

# The dotplot shows the point and interval estimates for the random effects,
# ordering them and highlighting which are significantly different from the mean (0)


# Prediction
#-------------
# Let's now examine standard predictions vs. subject-specific predictions.
# As with most R models, we can use the predict function on the model object.

# Prediction from regression model
lm1 <-
     lm(visual ~ -1 + visual0 + time.f + treat.f:time.f, data = armd)
summary(lm1)

predict_lm <- predict(lm1)
head(predict_lm)

# Prediction from mixed model on the training set:
# 1) Without random effects ->  re.form=NA
predict_no_re <- predict(fm16.1mer, re.form = NA)
head(predict_no_re) # (almost) same predictions
# 2) With random effects
predict_re <- predict(fm16.1mer)
head(predict_re)

# Prediction from mixed model on a test observation from a subject present in the training set:
test.data <-
     data.frame(
          subject = "234",
          treat.f = "Active",
          visual0 = 63,
          time = 12
     )

# 1) Without random effects ->  re.form=NA
predict_no_re <-
     predict(fm16.1mer, newdata = test.data, re.form = NA)
predict_no_re # (9.28808 + 0.82644*test.data$visual0 -0.21222*test.data$time -2.42200  -0.04959*test.data$time )

# 2) With random effects
predict_re <- predict(fm16.1mer, newdata = test.data)
predict_re # (9.28808 + 0.82644*test.data$visual0 -0.21222*test.data$time -2.42200  -0.04959*test.data$time -3.33466872 )

# where -3.33466872 comes from the random intercept vector and corresponds to the subject 234
re <- ranef(fm16.1mer)[[1]]
re[row.names(re) == test.data$subject, ]

# Prediction from mixed model on a test observation from a subject not present in the training set:
test.data <-
     data.frame(
          subject = "400",
          treat.f = "Active",
          visual0 = 63,
          time = 12
     )

# 1) Without random effects ->  re.form=NA
predict_no_re <-
     predict(fm16.1mer, newdata = test.data, re.form = NA)
predict_no_re # the same as before

# 2) With random effects
predict_re <- predict(fm16.1mer, newdata = test.data)
# it does not recognize the subject --> allow.new.levels = T
predict_re <-
     predict(fm16.1mer, newdata = test.data, allow.new.levels = T)
predict_re # the same as before, it uses the average of the random intercept, i.e. 0



# Diagnostic plots
#--------------------
# 1) Assessing Assumption on the within-group errors
x11()
plot(fm16.1mer) ## Pearson and raw residuals are the same now

x11()
qqnorm(resid(fm16.1mer))
qqline(resid(fm16.1mer), col = "red", lwd = 2)

# 2) Assessing Assumption on the Random Effects
x11()
qqnorm(unlist(ranef(fm16.1mer)$subject), main = "Normal Q-Q Plot - Random Effects on Intercept")
qqline(unlist(ranef(fm16.1mer)$subject), col = "red", lwd = 2)





###################################################################################################
## Model 2: random intercept + slope and homoscedastic residuals
## Model 2.1: general D

fm16.2mer <-
     lmer(
          visual ~ visual0 + time * treat.f + (1 + time | subject),
          data = armd,
          control = lmerControl(
               optimizer = "bobyqa",
               optCtrl = list(maxfun = 2e5)
          )
     )

summary(fm16.2mer)
confint(fm16.2mer, oldNames = TRUE)

vcovb <- vcov(fm16.2mer)
vcovb
corb <- cov2cor(vcovb)
nms <- abbreviate(names(fixef(fm16.2mer)), 5)
rownames(corb) <- nms
corb

## Var-Cov matrix of random-effects and errors
print(vc <- VarCorr(fm16.2mer), comp = c("Variance", "Std.Dev."))


## Let's compute the conditional and marginal var-cov matrix of Y
sgma <- summary(fm16.2mer)$sigma

A <- getME(fm16.2mer, "A") # A : N*2 x n
I.n <- Diagonal(ncol(A)) # IN: n x n

## the conditional variance-covariance matrix of Y (diagonal matrix)
SigmaErr <- sgma^2 * (I.n)
SigmaErr[3:6, 3:6] ## visualization of individual 2
# Conditioned to the random effects b_i, we observe the var-cov of the errors
# that are independent and homoscedastic

plot(as.matrix(SigmaErr[1:20, 1:20]), main = "Conditional estimated Var-Cov matrix of Y")

## the marginal variance-covariance matrix of Y (block-diagonal matrix)
V <- sgma^2 * (I.n + crossprod(A)) # V = s^2*(I_N+A*A)
V[3:6, 3:6] #-> V is a block-diagional matrix, the marginal var-cov matrix

# visualization of the first 20 rows/columns
plot(as.matrix(V[1:20, 1:20]), main = "Marginal estimated Var-Cov matrix of Y")


# PVRE
#--------------------
# In this case the variance of random sigma2_R effects represents the mean random
# effect variance of the model and is given by
# sigma2_b = Var(b0,b1) = sigma2_b0 + 2Cov(b0,b1)*mean(w) + sigma2_b1*mean(w^2)
# See equation (10) in Johnson (2014), Methods in Ecology and Evolution, 5(9), 944-946.
sigma2_eps <- as.numeric(get_variance_residual(fm16.2mer))
sigma2_eps
sigma2_b <-
     as.numeric(get_variance_random(fm16.2mer)) # 49.933917 + 0.074552*mean(armd$time^2) +2*0.143*7.06639*0.27304*mean(armd$time)
sigma2_b

PVRE <- sigma2_b / (sigma2_b + sigma2_eps)
PVRE # 72% is very high!

## visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i, b_1i for i=1,...,234
dotplot(ranef(fm16.2mer, condVar = T))


# Comparing models
#------------------
# The anova function, when given two or more arguments representing fitted models,
# produces likelihood ratio tests comparing the models.
anova(fm16.1mer, fm16.2mer)

# The p-value for the test is essentially zero -> we prefer fm16.2mer


# Diagnostic plots
#--------------------
# 1) Assessing Assumption on the within-group errors
x11()
plot(fm16.2mer)

x11()
qqnorm(resid(fm16.2mer))
qqline(resid(fm16.2mer), col = "red", lwd = 2)

# 2) Assessing Assumption on the Random Effects
x11()
qqnorm(unlist(ranef(fm16.2mer)$subject[, 1]), main = "Normal Q-Q Plot - Random Effects on Intercept")
qqline(unlist(ranef(fm16.2mer)$subject[, 1]), col = "red", lwd = 2)

x11()
qqnorm(unlist(ranef(fm16.2mer)$subject[, 2]), main = "Normal Q-Q Plot - Random Effects on Slope")
qqline(unlist(ranef(fm16.2mer)$subject[, 2]), col = "red", lwd = 2)


## We observe that the correlation between d_11 and d_22 id very low,
## we fit a new model with a diagonal D matrix

###################################################################################################
## Model 2.2: diagonal D

fm16.2dmer <-
     lmer(
          visual ~ visual0 + time * treat.f + (1 |
               subject) + (0 + time | subject),
          data = armd,
          control = lmerControl(
               optimizer = "bobyqa",
               optCtrl = list(maxfun = 2e5)
          )
     )

summary(fm16.2dmer)
confint(fm16.2dmer, oldNames = TRUE)

vcovb <- vcov(fm16.2dmer)
vcovb
corb <- cov2cor(vcovb)
nms <- abbreviate(names(fixef(fm16.2dmer)), 5)
rownames(corb) <- nms
corb

## Var-Cov matrix of random-effects and errors
print(vc <- VarCorr(fm16.2dmer), comp = c("Variance", "Std.Dev."))


## Let's compute the conditional and marginal var-cov matrix of Y
sgma <- summary(fm16.2dmer)$sigma

A <- getME(fm16.2dmer, "A") # A
I.n <- Diagonal(ncol(A)) # IN

## the conditional variance-covariance matrix of Y (diagonal matrix)
SigmaErr <- sgma^2 * (I.n)
SigmaErr[3:6, 3:6] ## visualization of individual 2
# Conditioned to the random effects b_i, we observe the var-cov of the errors
# that are independent and homoscedastic

plot(as.matrix(SigmaErr[1:20, 1:20]), main = "Conditional estimated Var-Cov matrix of Y")

## the marginal variance-covariance matrix of Y (block-diagonal matrix)
V <- sgma^2 * (I.n + crossprod(A)) # V = s^2*(I_N+A*A)
V[3:6, 3:6] #-> V is a block-diagional matrix, the marginal var-cov matrix

# visualization of the first 20 rows/columns
plot(as.matrix(V[1:20, 1:20]), main = "Marginal estimated Var-Cov matrix of Y")


# PVRE
#--------------------
# In this case the variance of random sigma2_R effects represents the mean random
# effect variance of the model and is given by
# sigma2_b = Var(b0,b1) = sigma2_b0 + 0 + sigma2_b1*mean(z^2)
# See equation (10) in Johnson (2014), Methods in Ecology and Evolution, 5(9), 944-946.
sigma2_eps <- as.numeric(get_variance_residual(fm16.2dmer))
sigma2_eps
sigma2_b <-
     as.numeric(get_variance_random(fm16.2dmer)) + mean(armd$time^2) * as.numeric(get_variance_slope(fm16.2dmer)) # 54.07117 + 0.07935904*mean(armd$time^2)
sigma2_b

PVRE <- sigma2_b / (sigma2_b + sigma2_eps)
PVRE # 72% is very high!

## visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i, b_1i for i=1,...,234
dotplot(ranef(fm16.2mer, condVar = T))


# Comparing models
#------------------
# The anova function, when given two or more arguments representing fitted models,
# produces likelihood ratio tests comparing the models.
anova(fm16.2mer, fm16.2dmer)

# The p-value for the test is essentially zero -> we prefer fm16.2dmer






##########################################################
### LINEAR MIXED MODELS WITH HETEROSCEDASTIC RESIDUALS ###
############ nlme package --> lme() function #############
##########################################################

######################################################################################################
## Model 1. Random intercept, heteroscedastic residuals (varPower of time)

## fixed-effects formula
lm2.form <-
     formula(visual ~ visual0 + time + treat.f + treat.f:time)

# LMM with homoscedastic residuals
fm16.1 <- lme(lm2.form, random = ~ 1 | subject, data = armd)

# update fm16.1 including heteroscedastic residuals
fm16.2 <- update(fm16.1,
     weights = varPower(form = ~time),
     data = armd
)

summary(fm16.2)

VarCorr(fm16.2)

## var-cov matrix of the errors (i.e. of Y, conditional to the random effects), that are independent but heteroscedastic
fm16.2ccov <-
     getVarCov(fm16.2, type = "conditional", individual = "2")
fm16.2ccov

plot(as.matrix(fm16.2ccov[[1]]), main = expression(paste(
     "Conditional estimated Var-Cov matrix of ", Y[2]
)))

## var-cov matrix of Y_i
fm16.2cov <- getVarCov(fm16.2, type = "marginal", individual = "2")
fm16.2cov # (90.479 = 31.103 + 59.37555; 121.440 = 62.062 + 59.37555; ...)
plot(as.matrix(fm16.2cov[[1]]), main = expression(paste("Marginal estimated Var-Cov matrix of ", Y[2])))

# var-cov matrix of y_i is the same for each subject i,
# except for the number of observations, ranging from 1 to 4

## correlation matrix of Y_i
cov2cor(fm16.2cov[[1]])

## ANALYSIS OF RESIDUALS
# Default residual plot of conditional Pearson residuals
plot(fm16.2)

# Plots (and boxplots) of Pearson residuals per time and treatment
plot(fm16.2, resid(., type = "pearson") ~ time | treat.f,
     id = 0.05
)
bwplot(resid(fm16.2, type = "p") ~ time.f | treat.f,
     panel = panel.bwplot,
     # User-defined panel (not shown)
     data = armd
)
# Despite standardization, the variability of the residuals seems to vary a bit.

# Normal Q-Q plots of Pearson residuals
qqnorm(fm16.2, ~ resid(.) | time.f)

## ANALYSIS OF RANDOM EFFECTS
# Normal Q-Q plots of predicted random effects
qqnorm(fm16.2, ~ ranef(.))

## Computing predictions comparing population average predictions with patient-specific predictions

aug.Pred <- augPred(
     fm16.2,
     primary = ~time,
     # Primary covariate
     level = 0:1,
     # fixed/marginal (0) and subj.-spec.(1)
     length.out = 2
) # evaluated in two time instants (4 e 52 wks)

plot(aug.Pred, layout = c(4, 4, 1))


######################################################################################################
## Model 2.1. random intercept + slope (correlated), heteroscedastic residuals (varPower of time)

fm16.3 <- update(fm16.2,
     random = ~ 1 + time | subject,
     data = armd
)
summary(fm16.3)

getVarCov(fm16.3, individual = "2") # D_i italic (i=2)

intervals(fm16.3, which = "var-cov") # Estimate of theta_D, delta e sigma

######################################################################################################
## Model 2.2. random intercept + slope independent, heteroscedastic residuals (varPower of time)
fm16.4 <- update(fm16.3,
     random = list(subject = pdDiag(~time)), # Diagonal D
     data = armd
)
summary(fm16.4) ## results suggest to remove the Treat and Time interaction

intervals(fm16.4)

anova(fm16.4, fm16.3) # We test if d_12 = 0 --> d_12 is not statistically different from 0, we can simplify the D structure in diagonal

qqnorm(fm16.4, ~ ranef(.)) # to be interpreted with caution since it might not reflect the real unknown distribution

plot(fm16.4, resid(., type = "pearson") ~ time | treat.f,
     id = 0.05
)
bwplot(resid(fm16.4, type = "p") ~ time.f | treat.f,
     panel = panel.bwplot,
     # User-defined panel (not shown)
     data = armd
)

## We make predictions comparing population average predictions with patient specific predictions
aug.Pred <- augPred(
     fm16.4,
     primary = ~time,
     # Primary covariate
     level = 0:1,
     # Marginal(0) and subj.-spec.(1)
     length.out = 2
) # Evaluated in two time instants (4 e 52 wks)

plot(aug.Pred, layout = c(4, 4, 1), columns = 2)


## let's compare the 4 fitted models
AIC(fm16.1, fm16.2, fm16.3, fm16.4)
anova(fm16.1, fm16.2, fm16.3, fm16.4)
