setwd("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 10 - 10052021")

library(MASS)
library(car)

#_______________________________________________________________________________
##### Problem 2 of 28/2/2007
#####--------------------------------------
# Pb2.txt dataset shows average monthly temperature (°C) recorded in
# 2006 in three Canadian locations: Edmonton, Montreal and Resolute.
# It is common in meteorology to assume that the average monthly 
# temperatures fluctuate sinusoidally around an annual average value:
# Temp.g (t) =  beta0.g +beta1.g * sin (2pi / 12 * t) + 
#               beta2.g * cos (2pi / 12 * t) + eps
# with eps ~ N (0,2), t = 1, 2, 3,. . . , 12 (month) and g = Edmonton, 
# Resolute, Montreal (Location).
# (a) Using the least squares method is estimate the 10 parameters of the model
# (b) Verify the model assumptions.
# (c) Taking advantage of the known trigonometric relation
#     sin(alpha-beta) = sin(alpha) * cos(beta) - cos(alpha) * sin(beta)
#     and reinterpreting the model of the form:
#     Temp.g (t) = mu.g + A.g * sin (2pi / 12 * (t-phi.g) + eps
#     report the analytical relation between the new parameters
#     (mu.g, A.g, phi.g) and the old parameters (beta0.g, beta1.g, beta2.g).
# (d) Estimate the parameters of the new formulation, namely:
#     - The annual average values (mu.g).
#     - The oscillation amplitudes (A.g).
#     - The phases of the oscillations (phi.g).
# (e) Through the use of an appropriate statistical test (report the
#     corresponding p-value) justify the possibility of using a reduced model
#     that assume that the fluctuations have same amplitude and phase, but 
#     different annual means in the stations of in Edmonton and Montreal.

Temperature <- read.table('Pb2.txt')
Temperature

x11()
matplot(Temperature, type='l', lwd=2)
dev.off()

temp <- cbind(Tem = c(Temperature[,1],Temperature[,2],Temperature[,3]),
              Sin = sin(2*pi/12*c(1:12,1:12,1:12)),
              Cos = cos(2*pi/12*c(1:12,1:12,1:12)),
              Res = c(rep(0,12), rep(1,12), rep(0,12)), # dummy for Resolute
              Mon = c(rep(0,12), rep(0,12), rep(1,12))) # dummy for Montreal

temp <- data.frame(temp)
temp


### question (a)

fit <- lm(Tem ~ Res + Mon + Sin + Cos + Res*Sin + Res*Cos +Mon*Sin + Mon*Cos, data=temp)
summary(fit)

xplot <- rep(seq(1,12,len=100),3)
nuovi <- data.frame(Sin = sin(2*pi/12*xplot), Cos = cos(2*pi/12*xplot),
                    Res = c(rep(0,100), rep(1,100), rep(0,100)),
                    Mon = c(rep(0,200), rep(1,100)))

x11()
plot(xplot, predict(fit, newdata = nuovi) ,col=rep(1:3,each=100), pch=16)
points(rep(1:12,3), temp$Tem, col='blue', lwd=2)

dev.off()

coef(fit)
Beta <- rbind(
  Edmonton = coef(fit)[c(1,4,5)],
  Resolute = coef(fit)[c(1,4,5)] +  coef(fit)[c(2,6,7)],
  Montreal = coef(fit)[c(1,4,5)] +  coef(fit)[c(3,8,9)])
Beta

### question (b)
shapiro.test(residuals(fit))
shapiro.test(rstudent(fit))

x11()

plot(fitted(fit),residuals(fit))
plot(rep(1:12,3),residuals(fit))

dev.off()

### question (c)
nuovi <- cbind(media = Beta[,1], 
               ampiezza = sqrt(Beta[,2]^2+Beta[,3]^2),
               fase = -12/(2*pi)*atan(Beta[,3]/Beta[,2])+6)

### question (d)
nuovi

### question (e)
linearHypothesis(fit,
                 rbind(c(0,0,0,0,0,0,0,1,0),
                       c(0,0,0,0,0,0,0,0,1)),
                 c(0,0))

#_______________________________________________________________________________
##### Problem 4 of 1/7/2009
#####-------------------------
# The Index Librorum Prohibitorum (edition of 1948), lists about 10000 works 
# considered heretical by the Catholic Church lists. The file 'index.txt'
# shows, for the years ranging from 1300 to 1899, the number of works
# added annually to the Index. Most historians believe that the average
# number of works added each year decreased linearly in time during this
# period (model A). Recently, Prof. Langdon proposed a theory according to
# which the linear trend "momentarily" changed (in a discontinuous way) 
# during the French hegemony period (1768, the Treaty of Versailles, 1815, 
# Battle of Waterloo, included), during which a collapse of the works annually 
# added to the Index occurred (Model B). Defining as mu(t) the average number 
# of works added to the Index in year t, and formalizing the two models as 
# follows:
# Model A: mu(t) = alpha + beta * t;
# Model B: mu(t) = alpha1 + beta1 * t for 1768 <= t <= 1815
#          mu(t) = alpha2 + beta2 * t for t <= 1767 or t> = 1816;
# answer the following questions:
# a) estimate the parameters of both models using the method of least squares;
#    which assumptions needs to be introduced in order to get unbiased estimates?
# b) is there statistical evidence of a different linear trend in the period of
#    French hegemony?
# c) using the model (b) and Bonferroni's inequality, they provide 2 90% 
#    global confidence intervals for the mean and the variance of the number
#    of works included in the index in the year 1800;
# d) using the model (b), provide a 90% confidence interval for the difference 
#    between the average number of works added to the Index in 1816 and average 
#    number of works added in 1815.

index <- read.table('index.txt', header=TRUE)
index

rm(Anno); rm(Numero)
attach(index)
x11()
plot(Anno,Numero)

dev.off()

### question a)

# Model A
# Y = alpha + beta*t + eps, E[eps]=0, var(eps)=sigma^2
fitA <- lm(Numero ~ Anno)

summary(fitA)

# Model B
# Y = alpha.g + beta.g*t + eps =
#   = b0 + b1*D + b2*anno + b3*D*anno + eps, E[eps]=0, var(eps)=sigma^2

D <- ifelse(Anno>=1768 & Anno<=1815, 1, 0)

fitB <- lm(Numero ~ D + Anno + D*Anno )
summary(fitB)

alpha <- c(alpha1=coef(fitB)[1]+coef(fitB)[2],alpha2=coef(fitB)[1])
alpha
beta  <- c( beta1=coef(fitB)[3]+coef(fitB)[4], beta2=coef(fitB)[3])
beta

x11()
plot(Anno,Numero)
points(Anno, fitted(fitA), pch=19, col='blue')
points(Anno, fitted(fitB), pch=19)

dev.off()

### question b)

shapiro.test(residuals(fitB))

linearHypothesis(fitB, rbind(c(0,1,0,0),c(0,0,0,1)), c(0,0))

### question c)
k <- 2
alpha <- .1
n <- dim(index)[1]
r <- 3

Z0   <- data.frame(D=1, Anno=1800)
ICBmean <- predict(fitB, Z0, interval='confidence',level=1-alpha/k) 
ICBmean

e <- residuals(fitB)
ICBvar <- data.frame(L=t(e)%*%e/qchisq(1-alpha/(2*k),n-(r+1)),
                     U=t(e)%*%e/qchisq(alpha/(2*k),n-(r+1)))
ICBvar

### question d)
a <- c(0,-1,1,-1815)
Bf <- c('1816-1815_L'= t(a) %*% coefficients(fitB) - sqrt(t(a) %*% vcov(fitB) %*% a) * qt(1 - alpha/2, n-(r+1)),
        '1816-1815_U'= t(a) %*% coefficients(fitB) + sqrt(t(a) %*% vcov(fitB) %*% a) * qt(1 - alpha/2, n-(r+1)) )
Bf

detach(index)


#_______________________________________________________________________________
##### Problem 4 of 4/7/2007
#####-------------------------
# At the Tenaris steel mills, the relationship between length [m] and
# Temperature [°C] of some steel bars that will be sold to Pirelli
# is under study (the data are contained in tenaris.txt file). The relation
# is hypothesized of the kind:
#   L = L0 + C* T + D  * T ^ 2 + eps
# with L the length of the bar, T the temperature of the bar, L0 the length 
# of the bar at 0 °C, C the coefficient of of linear thermal expansion, D
# the coefficient of quadratic thermal expansion and eps a measurement error
# of zero mean.
# Answer the following questions using appropriate statistical arguments:
# a) Estimate the parameters L0, C, D and the variance of error eps.
# b) Based on the analysis of residuals, do you think that there are the
#    conditions to make inference on the coefficients based on a Gaussian
#    model? (In case of Yes proceed to step (c); in case of negative answer
#    identify the problem, remove it and return to point (a))
# c) Do you think that the model explains the possible dependence between 
#    the temperature T and the length L?
# d) do you deem plausible to consider that the length of the bars at 0 °C
#    is equal to 2?
# E) do you think that you can eliminate from the model the quadratic term?

ten <- read.table('tenaris.txt', header=TRUE)
ten

attach(ten)

### question a)

fit <- lm(L ~ T + I(T^2))

summary(fit)

e <- residuals(fit)
S2 <- t(e)%*%e / (df.residual(fit))
S2

### question b)

shapiro.test(e)

x11()
par(mfrow=c(2,2))
plot(fit)

x11()
plot(T,L)
points(ten[1,1],ten[1,2],pch=19)

graphics.off()

detach(ten)

# Remove the outlier
ten1 <- ten[-1,]

fit <- lm(L ~ T + I(T^2), data=ten1)
summary(fit)

e <- residuals(fit)
S2 <- t(e)%*%e / (df.residual(fit))
S2

shapiro.test(e)

x11()
par(mfrow=c(2,2))
plot(fit)

dev.off()

### question c)

attach(ten1)
x11()
plot(T,L)
points(T,fitted(fit),col='blue', pch=19)

dev.off()

summary(fit)$r.squared

### question d)

linearHypothesis(fit, c(1,0,0), 2)

### question e)
summary(fit)

# or
linearHypothesis(fit, c(0,0,1), 0)

detach(ten1)


#_______________________________________________________________________________
##### Problem 4 of 09/09/2009
#####--------------------------------------
# A zoologist is studying the temporal evolution of the height of a new breed
# of goat (dataset goat.txt). Using an exponential growth model for the i-th 
# individual of the following type: 
# h_i = A + B (1 - exp(-t_i)) + C eps_i, 
# with: 
# h_i = the height of the individual [cm],
# t_i = the age of the individual [years], 
# eps_i = a random term distributed according to a standard normal distribution, 
# A and B = parameters exclusively dependent on the gender of the individual, 
# C = a parameter  equal for the whole population.
# a) Estimate with the least squares method the 5 parameters of the model.
# b) On the basis of model (a), is there statistical evidence that the mean 
#    height at birth (t=0) is different between males and females?
# c) On the basis of model (a), is there statistical evidence that the mean 
#    height at adulthood (t=+inf) is different between males and females?
# d) On the basis of tests (b) and (c), implement an appropriate reduced model and 
#    estimate its parameters.
# e) On the basis of the reduced model (d), provide confidence intervals of global 
#    confidence 90% for the mean height at birth and at adulthood of males and females.

goat <- read.table('goat.txt', header=T)
head(goat)
dim(goat)

plot(goat$age, goat$height, col=as.numeric(goat$gender)+1, pch=20)

n <- dim(goat)[1]

# a)
fit <- lm(height ~ gender + I(1 - exp(-age)) + I(1 - exp(-age)):gender, data=goat)
summary(fit)

# b)
library(car)
linearHypothesis(fit,rbind(c(0,1,0,0)),0)

# c)
linearHypothesis(fit,rbind(c(0,1,0,1)),0)

# d)
fit.red <- lm(height ~ I(1 - exp(-age)) + I(1 - exp(-age)):gender, data=goat)
summary(fit.red)

# e)
C <- rbind(c(1,0,0), c(1,1,0), c(1,1,1))

Bf <- rbind(
  c((C %*% coefficients(fit.red))[1] - sqrt((C %*% vcov(fit.red) %*% t(C))[1,1]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[1] + sqrt((C %*% vcov(fit.red) %*% t(C))[1,1]) * qt(1 - 0.10/6, n-3)),
  c((C %*% coefficients(fit.red))[2] - sqrt((C %*% vcov(fit.red) %*% t(C))[2,2]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[2] + sqrt((C %*% vcov(fit.red) %*% t(C))[2,2]) * qt(1 - 0.10/6, n-3)),
  c((C %*% coefficients(fit.red))[3] - sqrt((C %*% vcov(fit.red) %*% t(C))[3,3]) * qt(1 - 0.10/6, n-3),
    (C %*% coefficients(fit.red))[3] + sqrt((C %*% vcov(fit.red) %*% t(C))[3,3]) * qt(1 - 0.10/6, n-3))
)
