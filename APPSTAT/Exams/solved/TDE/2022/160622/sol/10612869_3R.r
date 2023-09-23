# We are interested in studing the danceability (Y ) of a song with respect to other features of the song. The file
# danceability.txt contains the values of danceability (the higher the value, the easier it is to dance to this song),
# loudness, energy, tempo and genre of 400 songs. Consider a linear model of the form:
#   Y = β0 + β1 · loudness + β2 · energy + β3 · tempo + ε
# with ε ∼ N (0, σ 2 ).
# a) Provide an estimate of the βi , i = 0, . . . , 3, and of σ.
# b) State and verify the model assumptions.
# c) Perform a test of level 5% to verify if loudness and energy can be both discarded from the model.
# d) Perform any other statistical tests that you consider useful to reduce the model, and update the estimates of
# its parameters.
# e) Let’s consider now the variable genre in the model as a random intercept. Fit a suitable model for accounting
# the hierarchy and compute and report the PVRE index.
# f) Report the dot plot of the estimated random intercepts. Net to the effect of fixed effect covariates, which is the
# genre associated to the highest danceability?

library(car)

danceability <- read.table('danceability.txt', header=T)

danceability$genre <- as.factor(danceability$genre)

# a) Provide an estimate of the βi , i = 0, . . . , 3, and of σ.

fit <- lm (danceability ~ . - genre, data=danceability)
summary(fit)
sqrt(sum(fit$residuals^2)/fit$df)

# b) State and verify the model assumptions.
par(mfrow=c(2,2))
plot(fit)
shapiro.test(fit$residuals)


# c) Perform a test of level 5% to verify if loudness and energy can be both discarded from the model.

#test H0: (β1, β2) == (0, 0) vs H1: (β1, β2) != (0, 0)
linearHypothesis(fit, rbind(c(0,1,0,0), c(0,0,1,0)), c(0,0)) 
#p-value is essentially zero, we reject the null hypothesis


# d) Perform any other statistical tests that you consider useful to reduce the model, and update the estimates of
# its parameters.
# We perform the test H0: β2 == 0 vs H1: β2 != 0
# p-value is 0.0731, we accept the null hypothesis: we erase from our model the dependency on energy

fit2 <- lm (danceability ~ . - genre - energy, data=danceability)
summary(fit2)
sqrt(sum(fit2$residuals^2)/fit2$df)


# e) Let’s consider now the variable genre in the model as a random intercept. Fit a suitable model for accounting
# the hierarchy and compute and report the PVRE index.

library(ggplot2)
library(insight)
library(lattice)
library(lme4)


# MODEL: danceability_ij = β0 + β1*loudness_ij + β2*energy_ij + β3*tempo_ij  + b_i + eps_ij
# eps_ij ~ N(0, sigma2_eps)
# b_i ~ N(0, sigma2_b)

lmm1 = lmer(danceability ~ loudness + energy + tempo + (1|genre), 
            data = danceability)
summary(lmm1)

sigma2_eps <- as.numeric(get_variance_residual(lmm1))
sigma2_eps
sigma2_b <- as.numeric(get_variance_random(lmm1))
sigma2_b

PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
PVRE

# f) Report the dot plot of the estimated random intercepts. Net to the effect of fixed effect covariates, which is the
# genre associated to the highest danceability?

dotplot(ranef(lmm1))
fixef(lmm1)





