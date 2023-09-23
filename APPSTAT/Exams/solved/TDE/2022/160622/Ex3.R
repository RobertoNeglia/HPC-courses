### Ex3


library(MASS)
library(car)
library(rgl)


dataset = read.table("danceability.txt", header = T)

attach(dataset)

n          <- dim(dataset)[[1]]

fm <- lm(danceability ~ loudness + energy + tempo)  # response ~ regressors

summary(fm)

b = fm$coefficients
b


# Assumption: Eps ~ N(0, sigma^2)

x11()
par(mfrow=c(2,2))
plot(fm)
# first plot: the red line is the estimated mean and we want it to be near 0
#             moreover we wand that sigma does not deoend on x
# Scale location: not imp
# Res vs lev: to see if there are leverages
# QQ plot

shapiro.test(residuals(fm))

C = rbind(c(0,1,0,0), c(0,0,1,0))
mu0 = c(0,0)
linearHypothesis(fm, C, mu0) 

fit <- lm(danceability ~ loudness  + tempo)
summary(fit)


library(nlmeU)
library(corrplot)
library(nlme)
library(lattice)
library(plot.matrix)
library(lme4)  # only handle independence an homosched residuals
# at the end there is an extension
library(insight)



# e) 

lmm.1 <- lmer(danceability ~ loudness + tempo  + (1|genre),    
              data = dataset)
summary(lmm.1)

confint(lmm.1,oldNames=TRUE)


sigma2_eps <- as.numeric(get_variance_residual(lmm.1))
sigma2_eps    # sigma^2 dei residui
sigma2_b <- as.numeric(get_variance_random(lmm.1))
sigma2_b 

PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
PVRE

fixef(lmm1)

ranef(lmm1)  # estrae tutte le random intercept
x11()
dotplot(ranef(lmm.1, condVar=T))

