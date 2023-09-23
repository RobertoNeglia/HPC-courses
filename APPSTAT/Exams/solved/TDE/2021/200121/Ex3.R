## Ex3


library(MASS)
library(car)
library(rgl)


dataset = read.table("bikes.txt", header = T)

# categorica con 2 livelli (eg : sex )
group = as.factor(dataset$day)  # Tipo= nome delle colonna della var categorica
dummy = ifelse(group == 'Holiday', 1, 0)     # Red è uno dei gruppi    

# 0 -> no holiday

y = dataset$bike_count     # target
z1 = dataset$mean_temp     # regressor 1 
z2 = dataset$mean_wind     # regressor 2

fit <- lm(y ~ z1 + z2 + dummy + z1:dummy + z2:dummy, data = dataset)
summary(fit)

b = fit$coefficients

Beta = rbind( Holiday = c(b[1]+b[4], b[2]+b[5], b[3]+b[6]), No_holiday = c(b[1], b[2], b[3]))
Beta

##### Verify assumptions
#####--------------------
x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))

# test weather
A <- rbind(c(0,1,0,0,0,0), c(0,0,1,0,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
# metto gli 1 in corrispondenza dei parametri della dummy e delle sue interazioni

mu0 <- c(0,0,0,0)   # lunghezza = regressori del modello che contengono la dummy
linearHypothesis(fit, A, mu0)
# reject

# test weather
A <- rbind(c(0,0,0,1,0,0), c(0,0,0,0,1,0), c(0,0,0,0,0,1))
# metto gli 1 in corrispondenza dei parametri della dummy e delle sue interazioni

mu0 <- c(0,0,0)   # lunghezza = regressori del modello che contengono la dummy
linearHypothesis(fit, A, mu0)



# c)
vif(fit)

fit <- lm(y ~ z1 + z2 + dummy + z1:dummy , data = dataset)
summary(fit)

fit <- lm(y ~ z1 + z2 + dummy , data = dataset)
summary(fit)

fit <- lm(y ~ z1 + dummy , data = dataset)
summary(fit)

# ok tutti significati ma R^2 molto basso.. idee?

b = fit$coefficients

Beta = rbind( Holiday = c(b[1]+b[3], b[2]), No_holiday = c(b[1], b[2]))
Beta

fm = fit


# d)

Z0.new <- data.frame(z1=2, dummy=1)

Pred <- predict(fm, Z0.new, interval='prediction', level=1-0.05)  
Pred



















