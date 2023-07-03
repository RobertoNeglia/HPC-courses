##############################################################
#############  Applied Statistics 2022/2023  #################
##############  Linear Mixed-effects models  #################
##############################################################

# Topics:
#   1. Linear Models with homoscedastic and independent errors
#   2. Linear Models with heteroscedastic and independent errors
#      2.1 VarIdent()
#      2.2 VarPower()
#   3. Linear Models with heteroscedastic and dependent errors
#      3.1 CorCompSym()
#      3.2 AR(1)
#      3.3 general

rm(list = ls())
graphics.off()


library(nlmeU) ## --> for the dataset
library(nlme) ## --> for models implementation

library(corrplot)
library(lattice)
library(plot.matrix)

data(armd) # Age-Related Macular Degeneration: dataset of interest
data(armd0) # Age-Related Macular Degeneration: dataset for visualization
help(armd)

one <- armd[armd$subject == 2, ]
head(one)

# The ARMD data arise from a randomized multi-center clinical trial comparing an experimental treatment (interferon-alpha)
# versus placebo for patients diagnosed with ARMD. Patients with ARMD progressively lose vision.
# The dataset contains information about 234 subjects, for which the visual level is measured up to 4 times.
# We are mainly interested in the effect of treatment on the visual acuity measurements.

## The ARMD0 dataset contains the same information, but with an extra row for each patient relative
## to the measurement at time 0. ARMD0 contains 240 subjects.

## Visual-acuity profiles for selected patients --> we visualize same of the trends
armd0.subset <-
   subset(armd0, as.numeric(subject) %in% seq(1, 240, 5)) # one each 5 patients

xy1 <- xyplot(
   visual ~ time | treat.f,
   groups = subject,
   data = armd0.subset,
   type = "l",
   lty = 1
)
update(xy1,
   xlab = "Time (in weeks)",
   ylab = "Visual acuity",
   grid = "h"
)
## We observe a decreasing trend in time, on average, but patients have very different trends.
## Also, Active patients have on average lower values of the response.


## sample means across time and treatment
flst <- list(armd$time.f, armd$treat.f)
tMn <- tapply(armd$visual, flst, FUN = mean)
tMn

## We confirm what we observe in the plot

## Box-plots for visual acuity by treatment and time
bw1 <- bwplot(visual ~ time.f | treat.f,
   data = armd0
)
xlims <- c("Base", "4\nwks", "12\nwks", "24\nwks", "52\nwks")
update(bw1, xlim = xlims, pch = "|")


##############################################################################################
# 1. Linear Models with homogeneous and independent errors: we start by considering all
#    observations as independent, with homogeneous variance

# LM: VISUAL_it = b_0t + b1 � VISUAL_0i + b_2t � TREAT_i + e_it

# b_0t, b_1, and b_2t denote the timepoint-specific intercept,
# baseline visual acuity effect, and timepoint-specific treatment effect.
# Thus, the model assumes a time-dependent treatment effect,
# with the time variable being treated as a factor.

# To obtain timepoint-specific intercepts at 4,12,24 and 52 weeks,
# the overall intercept is removed from the model by specifying the -1 term.

lm1.form <-
   lm(visual ~ -1 + visual0 + time.f + treat.f:time.f, data = armd)
summary(lm1.form)

# variance-covariance matrix of Y  --> it is a diagonal matrix with a value of 12.38^2
par(mar = c(4, 4, 4, 4))
plot(diag(x = 12.38^2, nrow = 30, ncol = 30), main = "Variance-covariance matrix of Y")

## residual analysis
plot(lm1.form$residuals) # they seem quite homoscedastic
abline(h = 0)

qqnorm(lm1.form$residuals)
qqline(lm1.form$residuals)

shapiro.test(lm1.form$residuals)

## But we know that observations are not independent and that the variance of the visual measurements increases in time

## let's color the residuals relative to different patients
colori <- rainbow(length(unique(armd$subject)))
num_sub <- table(armd$subject)
colori2 <- rep(colori, num_sub)
plot(lm1.form$residuals, col = colori2)
abline(h = 0) ## --> not very informative

boxplot(
   lm1.form$residuals ~ armd$subject,
   col = colori,
   xlab = "Subjects",
   ylab = "Residuals",
   main = "Distribution of residuals across patients"
) ## --> informative!

## let's color the residuals relative to different time instants
set.seed(1)
colori <- rainbow(4)
colori2 <-
   colori[armd$tp] # associate to each one of the 4 time instants a color
plot(lm1.form$residuals, col = colori2, ylab = "residuals")
abline(h = 0)
legend(
   650,
   -25,
   legend = c("time 4wks", "time 12wks", "time 24wks", "time 52wks"),
   col = colori,
   lty = 1,
   cex = 0.8
)

## Note: we observe that red points are the closest to 0, purple ones are the farthest
## We expect the residuals to be heterogeneous across different time instants observations


boxplot(
   lm1.form$residuals ~ armd$time.f,
   col = colori,
   xlab = "Time.f",
   ylab = "Residuals"
) ## -> the variance of th observations increases in time


# The model does not take into account the correlation
# between the visual acuity observations obtained from the same subject.
# It also does not take into account the heterogeneous variability
# present at different time points. Thus, it should not be used as a basis for inference.

############################################################################################
## 2. Linear models with heteroscedastic and independent errors

## We know that variance increases in time --> we model the variance as a function of time
## We have different possibilities

## gls() function allows the inclusion of dependency and heteroscedasticity

## 2.1 Option 1: VarIdent()

fm9.1 <-
   gls(
      visual ~ -1 + visual0 + time.f + treat.f:time.f,
      # the same as before
      weights = varIdent(form = ~ 1 |
         time.f),
      # Var. function; <delta, stratum>-group
      data = armd
   )
summary(fm9.1)

plot(fm9.1$residuals)

fm9.1$modelStruct$varStruct
intervals(fm9.1, which = "var-cov") ## 95% CI

# Visualization of Variance-covariance matrix of Y (first 30 observations)
par(mar = c(4, 4, 4, 4))
plot(
   diag(
      x = c(
         1.000000^2 * 8.244094^2,
         1.397600^2 * 8.244094^2,
         1.664321^2 * 8.244094^2,
         1.880852^2 * 8.244094^2
      ),
      nrow = 30,
      ncol = 30
   ),
   main = "Variance-covariance matrix of Y - VarIdent()"
)

## To formally test the hypothesis that the variances are timepoint specific,
## we apply the anova() function. The LR test tests the null hypothesis of homoscedasticity.

## The anova() function will take the model objects as arguments, and return an ANOVA testing
## whether the more complex model is significantly better at capturing the data than the simpler model.

anova(fm9.1, lm1.form) ## lm1.form C fm9.1

## 2.2 Option 2: VarPower()

## Now that we know the variance is increasing in time, we try a more parsimonious model

fm9.2 <-
   update(fm9.1, weights = varPower(form = ~time)) # Var. function; <delta, v_it>-group
summary(fm9.2)

fm9.2$modelStruct$varStruct
intervals(fm9.2, which = "var-cov")


# Visualization of Variance-covariance matrix of Y (first 30 observations)
par(mar = c(4, 4, 4, 4))
plot(
   diag(
      x = c(
         4^(2 * 0.2519332) * 5.974906^2,
         12^(2 * 0.2519332) * 5.974906^2,
         24^(2 * 0.2519332) * 5.974906^2,
         52^(2 * 0.2519332) * 5.974906^2
      ),
      nrow = 30,
      ncol = 30
   ),
   main = "Variance-covariance matrix of Y - VarIdent()"
)



# Test of the variance structure: power of time vs. timepoint-specific variances
anova(fm9.2, fm9.1)

AIC(fm9.2, fm9.1) # --> fm9.2 is better in terms of AIC and parsimony!




## Residual analysis --we assess the fit of the model using residual plots.

## raw residuals
plot(fm9.2, resid(., type = "response") ~ fitted(.)) # Raw vs. fitted
# We observe an asymmetric pattern, with large positive (negative) residuals present mainly for small (large) fitted values.
# but it can be a consequence of the fact that raw residuals are intrinsically heteroscedastic and correlated.

plot(fm9.2, resid(., type = "response") ~ time) # Raw vs. time (not shown)
bwplot(resid(fm9.2) ~ time.f, pch = "|", data = armd)
# The boxand-whiskers plots clearly show an increasing variance of the residuals.

## Pearson residuals
## Pearson residuals are obtained from the raw residuals by dividing the latter by an
## estimate of the appropriate residual standard deviation, so they should be more homoscedastic
plot(fm9.2, resid(., type = "pearson") ~ fitted(.)) # Pearson vs. fitted
plot(fm9.2, resid(., type = "pearson") ~ time)
bwplot(resid(fm9.2, type = "pearson") ~ time.f,
   # Pearson vs. time.f
   pch = "|",
   data = armd
)
## this plot illustrate the effect of scaling: the variance of the residuals is virtually constant.


#########################################################################################################
## 3. Linear models with heteroscedastic and dependent errors

## We now modify the model, so that the visual acuity measurements,
## obtained for the same individual, are allowed to be correlated.

## We can estimate the semivariogram to calculate correlation coefficients between Pearson
## residuals for every pair of timepoints, separately.

# The semivariogram function can be defined as the complement of the correlation function.

## Variogram per time difference
Vg1 <- Variogram(fm9.2, form = ~ time | subject)
Vg1
plot(Vg1,
   smooth = FALSE,
   xlab = "Time difference",
   ylim = c(0, 0.7)
)


## Variogram per time lag
Vg2 <- Variogram(fm9.2, form = ~ tp | subject)
Vg2
plot(Vg2,
   smooth = FALSE,
   xlab = "Time Lag",
   ylim = c(0, 0.7)
)

## From these two plots we see that correlation decreases with time lag/difference
## Therefore, a  correlation structure like, e.g., a compound symmetry, will most likely not fit the data well.
## A more appropriate structure might be, e.g., an autoregressive process of order 1 AR(1).

## Nevertheless, for illustrative purposes, we consider a model with a compound symmetry
## correlation structure.

## 3.1 Correlation 1: CorCompSym()
lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f)
fm12.1 <- gls(
   lm1.form,
   weights = varPower(form = ~time),
   correlation = corCompSymm(form = ~ 1 | subject),
   data = armd
)
summary(fm12.1)

intervals(fm12.1, which = "var-cov")
# With the estimates of rho, sigma and delta we can estimate the var-cov matrix

# The marginal variance-covariance structure
fm12.1vcov <-
   getVarCov(fm12.1, individual = "2") # estimate of R_i, e.g. i=2
nms <- c("4wks", "12wks", "24wks", "52wks")
dnms <- list(nms, nms) # Dimnames created
dimnames(fm12.1vcov) <- dnms # Dimnames assigned
print(fm12.1vcov)

matrix(as.numeric(fm12.1vcov), nrow = 4, ncol = 4, byrow = TRUE)

## on the diagonal we have (5.981515^2)*TIME^(2*0.2598167)
## out of the diagonal we have (5.981515^2)*TIME_1^(0.2598167)*TIME_2^(0.2598167)*rho


## Visualization of the marginal variance-covariance matrix of Y
# R_i <- rbind(
#    c(73.531, 56.077, 67.143, 82.081),
#    c(56.077, 130.140, 89.323, 109.200),
#    c(67.143, 89.323, 186.560, 130.740),
#    c(82.081, 109.200, 130.740, 278.810)
# )

R_i <- matrix(as.numeric(fm12.1vcov), nrow = 4, ncol = 4, byrow = TRUE)

R <- matrix(0, nrow = 28, ncol = 28)
for (i in 0:6) {
   R[(i * 4 + 1):(i * 4 + 4), (i * 4 + 1):(i * 4 + 4)] <- R_i
}
plot(R)

print(cov2cor(fm12.1vcov), corr = TRUE, stdevs = FALSE) ## Estimate of C_i (correlation matrix)


## Test of independence vs. compound-symmetry correlation structure
anova(fm9.2, fm12.1) # M9.2 C M12.1

# The result of the LR test is clearly statistically significant, indicating
# the importance of the adjustment for the correlation in modeling the data
# -> very small pvalue and also very much smaller AIC for fm12.1 wrt fm9.2

## 3.2 Correlation 2: AR(1)
fm12.2 <- update(fm9.2,
   correlation = corAR1(form = ~ tp | subject),
   data = armd
)
summary(fm12.2)
# Phi1 = 0.6573069 -> it's our estimate of rho
intervals(fm12.2, which = "var-cov")

# The marginal variance-covariance structure
fm12.2vcov <-
   getVarCov(fm12.2, individual = "2") # Estimate of R_i, e.g. i=2
dimnames(fm12.2vcov) <- dnms
fm12.2vcov

## on the diagonal we have (6.356295^2)*TIME^(2*0.2311874)
## out of the diagonal we have (6.35629^2)*4^(0.2311874)*12^(0.2311874)*0.6573069
##                             (6.35629^2)*4^(0.2311874)*24^(0.2311874)*0.6573069^2...


fm12.2cor <- cov2cor(fm12.2vcov) # Estimate of C_i
print(fm12.2cor,
   digits = 2,
   corr = TRUE,
   stdevs = FALSE
)

# Compound-symmetry vs. autoregressive correlation (nonnested models)
anova(fm12.1, fm12.2)
## we prefer AR(1)

## 3.3 Correlation 3: general correlation structure
fm12.3 <- update(fm12.2,
   correlation = corSymm(form = ~ tp |
      subject),
   ## the variance function is still VarPower()
   data = armd
)
summary(fm12.3)

intervals(fm12.3, # 95% CIs for rho, delta, sigma
   which = "var-cov"
)

fm12.3vcov <-
   getVarCov(fm12.3, individual = "2") ## Estimate of R_i (italic)
dimnames(fm12.3vcov) <- dnms
fm12.3vcov # fm12.3vcov[1,1] = sigma_2 * Lambda[1]^2 = 5.737927^2*(4^0.2712624)^2

fm12.3cor <- cov2cor(fm12.3vcov) ## Estimate of C_i
print(fm12.3cor, corr = TRUE, stdevs = FALSE)


## Autoregressive of order 1 vs. a general correlation structure
anova(fm12.2, fm12.3) # M12.2 C M12.3 --> we prefer M12.3

## Model-Fit Diagnostics

# (a) Plots (and boxplots) of raw residuals
panel.bwxplot0 <- function(x, y, subscripts, ...) {
   panel.grid(h = -1)
   panel.stripplot(x, y, col = "grey", ...)
   panel.bwplot(x, y, pch = "|", ...)
}
bwplot(
   resid(fm12.3) ~ time.f | treat.f,
   panel = panel.bwxplot0,
   ylab = "Residuals",
   data = armd
)
# The box-and-whiskers plots clearly show an increasing variance of the residuals with timepoint.
# This reflects the heteroscedasticity.


# (b) Plots of Pearson residuals vs. fitted values
# Pearson residuals are obtained from the raw residuals by dividing the latter by an
# estimate of the appropriate residual standard deviation, so they should be more homoscedastic

plot(fm12.3)
# Due to the correlation of the residuals corresponding to the measurements obtained
# for the same patient at different timepoints, the plot reveals a pattern, with a few
# large, positive residuals in the upper-left part and a few negative ones in the lower-right part.

## We therefore decide to visualize the residuals for each time instants
plot(
   fm12.3,
   resid(., type = "p") ~ fitted(.) | time.f
)
stdres.plot <-
   plot(
      fm12.3,
      resid(., type = "p") ~ jitter(time) | treat.f,
      id = 0.01,
      adj = c(-0.3, 0.5),
      grid = FALSE
   )
plot(update(
   stdres.plot,
   # Fig. 12.4
   xlim = c(-5, 59),
   ylim = c(-4.9, 4.9),
   grid = "h"
))
# The four scatterplots show a somewhat more balanced pattern.
