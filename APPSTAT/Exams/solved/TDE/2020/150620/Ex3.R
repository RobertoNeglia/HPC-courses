### Ex3


library(MASS)
library(car)
library(rgl)


dataset = read.table("airfoil.txt", header = T)

attach(dataset)

n = dim(dataset)[1]
dummy = rep(0,n)
dummy[which(velocity == 'L')] = 1

fit <- lm(sound ~ frequency + dummy + frequency:dummy , data = dataset)
summary(fit)

b = fit$coefficients

Beta = rbind( L = c(b[1]+b[3], b[2]+b[4]), H = c(b[1], b[2]))
Beta

sum(residuals(fit)^2)/fit$df

# b) 
C = rbind(c(0,1,0,0), c(0,0,0,1))
mu0 = c(0,0)
linearHypothesis(fit, C, mu0) 

C = rbind(c(0,0,1,0), c(0,0,0,1))
mu0 = c(0,0)
linearHypothesis(fit, C, mu0) 


fit <- lm(sound ~ frequency + dummy  , data = dataset)
summary(fit)

Beta = rbind( L = c(b[1]+b[3], b[2]), H = c(b[1], b[2]))
Beta

sum(residuals(fit)^2)/fit$df


# d)
Z0.new <- data.frame(frequency=15000, dummy=0)

# Conf. int. for the mean
alpha = 0.05
Conf <- predict(fit, Z0.new, interval='confidence', level=1-alpha)  
Conf

detach(dataset)




















