### Ex3
library(MASS)
library(car)
library(rgl)


dataset = read.table("wine.txt", header = T)

attach(dataset)

n = dim(dataset)[1]

dummy_red = rep(0,n)
dummy_red[which(type == 'Red')] = 1

dummy_white = rep(0,n)
dummy_white[which(type == 'White')] = 1


fit <- lm(alcohol ~ sugar + dummy_red + dummy_white + sugar:dummy_red + sugar:dummy_white, data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( red = c(b[1]+b[3], b[2]+b[5]), white = c(b[1]+b[4], b[2]+b[6]), rose = c(b[1], b[2]))
Beta



x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))


C = rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
mu0 = c(0,0,0,0)
linearHypothesis(fit, C, mu0) 


C = rbind(c(0,1,0,0,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
mu0 = c(0,0,0)
linearHypothesis(fit, C, mu0) 



C = rbind( c(0,0,1,0,0,0), c(0,0,0,1,0,0))
mu0 = c(0,0)
linearHypothesis(fit, C, mu0) 

fit <- lm(alcohol ~ sugar  + sugar:dummy_red + sugar:dummy_white, data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( red = c(b[1], b[2]+b[3]), white = c(b[1], b[2]+b[4]), rose = c(b[1], b[2]))
Beta


# prediction
Z0.new <- data.frame(sugar =20, dummy_red=1 , dummy_white = 0)

# Conf. int. for the mean
alpha = 0.01

# Pred. int. for a new obs
Pred <- predict(fit, Z0.new, interval='prediction', level=1-alpha)  
Pred

detach(dataset)








