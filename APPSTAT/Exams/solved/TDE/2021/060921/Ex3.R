### Ex3

library(MASS)
library(car)
library(rgl)


dataset = read.table("boats.txt", header = T)

attach(dataset)
n = dim(dataset)[1]
dummy_w = rep(0,n)
dummy_w[which(material == 'wood')] = 1
y = price

fit = lm (y ~ dummy_w + length + power + draught + crew + year)
summary(fit)
b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( wood = c(b[1]+b[2], b[3], b[4], b[5], b[6], b[7]), fiberglass =c(b[1], b[3], b[4], b[5], b[6], b[7]) )
Beta

sum(residuals(fit)^2)/fit$df


x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))

# b)
A <- rbind(c(0,0,1,0,0,0,0), c(0,0,0,1,0,0,0), c(0,0,0,0,1,0,0))
# metto gli 1 in corrispondenza dei parametri della dummy e delle sue interazioni

mu0 <- c(0,0,0)   # lunghezza = regressori del modello che contengono la dummy
linearHypothesis(fit, A, mu0)


# b)
A <- rbind(c(0,1,0,0,0,0,0), c(0,0,0,0,0,1,0))
# metto gli 1 in corrispondenza dei parametri della dummy e delle sue interazioni

mu0 <- c(0,0)   # lunghezza = regressori del modello che contengono la dummy
linearHypothesis(fit, A, mu0)


fit = lm (y ~ dummy_w + length + power + crew + year)
summary(fit)
b = fit$coefficients

fit = lm (y ~ dummy_w + length + power + crew)
summary(fit)
b = fit$coefficients


Beta = rbind( wood = c(b[1]+b[2], b[3], b[4], b[5]), fiberglass =c(b[1], b[3], b[4], b[5]) )
Beta

sum(residuals(fit)^2)/fit$df


Z0.new <- data.frame(length=10, power=1070, crew=1, dummy_w=0)
alpha = 0.05
Pred <- predict(fit, Z0.new, interval='prediction', level=1-alpha)  
Pred











detach(dataset)










