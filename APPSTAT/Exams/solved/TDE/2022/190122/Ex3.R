### Ex3

library(MASS)
library(car)
library(rgl)

dataset = read.table("tattoo.txt", header = T)

attach(dataset)
n = dim(dataset)[1]



dummy_hand = rep(0,n)
dummy_hand[which(method == 'handmade')] = 1

fit <- lm(price ~ dimension + ncolors + dummy_hand + dimension:dummy_hand + ncolors:dummy_hand, data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( handmade = c(b[1]+b[4], b[2]+b[5], b[3]+b[6]), machine = c(b[1], b[2], b[3]))
Beta
sum(residuals(fit)^2)/fit$df


# b)
x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))

C = rbind(c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
mu0 = c(0,0,0)
linearHypothesis(fit, C, mu0)


C = rbind(c(0,0,1,0,0,0), c(0,0,0,0,0,1))
mu0 = c(0,0)
linearHypothesis(fit, C, mu0)


# d) 
fit <- lm(price ~ dimension + ncolors  + dimension:dummy_hand , data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( handmade = c(b[1], b[2]+b[4], b[3]), machine = c(b[1], b[2], b[3]))
Beta
sum(residuals(fit)^2)/fit$df


# d)
Z0.new <- data.frame(dimension=0, ncolors=0, dummy_hand = 1)

# Conf. int. for the mean
alpha = 0.05
Conf <- predict(fit, Z0.new, interval='confidence', level=1-alpha/2)  
Conf


Z0.new <- data.frame(dimension=6.5, ncolors=1, dummy_hand = 1)

# Conf. int. for the mean
alpha = 0.05
Conf <- predict(fit, Z0.new, interval='confidence', level=1-alpha/2)  
Conf














