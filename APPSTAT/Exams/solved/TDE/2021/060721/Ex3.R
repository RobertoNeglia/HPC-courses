### EX3
library(MASS)
library(car)
library(rgl)

dataset = read.table("pc.txt", header = T)


# categorica con g livelli -> creo g-1 dummy
group = as.factor(dataset$OS)  # Tipo= nome delle colonna della var categorica
dummy_mac = ifelse(group == 'Mac', 1, 0)     # Red è uno dei gruppi    
dummy_linux = ifelse(group == 'Linux', 1, 0) 

attach(dataset)

fit <- lm(price ~ freq + cache_acc + dummy_mac:freq + dummy_linux:freq + dummy_mac:cache_acc + dummy_linux:cache_acc, data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( Mac = c(b[1], b[2]+b[4], b[3]+b[6]), Linux = c(b[1], b[2]+b[5], b[3]+b[7]), Windows= c(b[1], b[2], b[3]))
Beta

sum(residuals(fit)^2)/fit$df


#####--------------------
x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))


C <- rbind(c(0,0,0,1,0,0,0), c(0,0,0,0,1,0,0), c(0,0,0,0,0,1,0), c(0,0,0,0,0,0,1))

linearHypothesis(fit, C, c(0,0,0,0))

# remove cache
C <- rbind(c(0,0,1,0,0,0,0), c(0,0,0,0,0,1,0), c(0,0,0,0,0,0,1))

linearHypothesis(fit, C, c(0,0,0))


# ok

fit <- lm(price ~ freq + dummy_mac:freq + dummy_linux:freq , data = dataset)
summary(fit)


fit <- lm(price ~ freq + freq:dummy_mac , data = dataset)
summary(fit)

b = fit$coefficients

#                        cambia i numeri!
Beta = rbind( Mac = c(b[1], b[2]+b[3]), Linux = c(b[1], b[2]), Windows= c(b[1], b[2]))
Beta


sum(residuals(fit)^2)/fit$df

# e) 
Z0.new <- data.frame(freq=3.2, dummy_mac=0)

# Conf. int. for the mean
alpha = 0.01
Conf <- predict(fit, Z0.new, interval='confidence', level=1-alpha)  
Conf












