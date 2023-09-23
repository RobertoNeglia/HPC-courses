### Ex3
library(MASS)
library(car)
library(rgl)


dataset = read.table("leaven.txt", header = T)

attach(dataset)

y = volume
z1 = time
z2 = time^2

n = dim(dataset)[1]

dummy_by = rep(0,n)
dummy_by[which(yeast == 'by')] = 1

fit <- lm(y ~ z1 + z2+ z1:dummy_by + z2:dummy_by)  # response ~ regressors

summary(fit)

b = fit$coefficients
b

Beta = rbind( by = c(b[1], b[2]+b[4], b[3]+b[5]), sd = c(b[1], b[2], b[3]))
Beta

sum(residuals(fit)^2)/fit$df

#####--------------------
x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))


C = rbind(c(0,0,0,1,0), c(0,0,0,0,1))
mu0 = c(0,0)
linearHypothesis(fit, C, mu0) 



fit2 <- lm(y ~ z1 + z1:dummy_by + z2:dummy_by)  # response ~ regressors

summary(fit2)

b = fit2$coefficients
b

Beta = rbind( by = c(b[1], b[2]+b[3], b[4]), sd = c(b[1], b[2], 0))
Beta

sum(residuals(fit2)^2)/fit2$df




# d)
Z0.by <- data.frame(z1=2, z2=4, dummy_by = 1)
Z0.sd <- data.frame(z1=2, z2=4, dummy_by = 0)

alpha = 0.05
Pred_by <- predict(fit2, Z0.by, interval='prediction', level=1-alpha)  
Pred_by

Pred_sd <- predict(fit2, Z0.sd, interval='prediction', level=1-alpha)  
Pred_sd

detach(dataset)


























