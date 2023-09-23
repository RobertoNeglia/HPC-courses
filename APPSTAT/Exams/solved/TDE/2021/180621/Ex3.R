### Ex 3


library(MASS)
library(car)
library(rgl)

dataset = read.table("students.txt", header = T)
dim(dataset)

dummy_sex = ifelse(dataset$gender == 'female', 1, 0)   


y = dataset$watchtv 

attach(dataset)


fit <- lm(y ~ age + height + dummy_sex + distance + siblings + computertime + exercisehours + musiccds + playgames , data = dataset)
summary(fit)

b = fit$coefficients
sum(residuals(fit)^2)/fit$df

x11()
par(mfrow=c(2,2))
plot(fit)

shapiro.test(rstandard(fit))


# b) 
library(rgl)

library(glmnet)

Z <- model.matrix(y ~ age + height + dummy_sex + distance + siblings + computertime + exercisehours + musiccds + playgames)[,-1]

fit.lasso <- glmnet(Z, y, lambda = 0.3) # default: alpha=1 -> lasso 

fit.lasso$beta    # senza intercetta

coef(fit.lasso)   # con intercetta



# c)
lambda.grid <- 10^seq(0,-2,length=100)
fit.lasso_grid<- glmnet(Z, y, lambda = lambda.grid) # default: alpha=1 -> lasso 


x11()
plot(fit.lasso_grid , xvar='lambda', label=TRUE, col = rainbow(dim(Z)[2]))
legend('topright', dimnames(Z)[[2]], col =  rainbow(dim(Z)[2]), lty=1, cex=1)


# CV
cv.lasso <- cv.glmnet(Z, y,lambda=lambda.grid) # default: 10-fold CV
# NB: if we have a low number of data we have to choose a low k for k-fold CV
# eg: n = 13 -> k=3 and we write in the argument before lambda: "nfolds=3"

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

x11()
plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)


# Get the coefficients for the optimal lambda
coef.lasso <- predict(fit.lasso_grid, s=bestlam.lasso, type = 'coefficients')[1:10,]
coef.lasso 



# no scelgo lambda = 1
fit.lasso_opt <- glmnet(Z, y, lambda = 1)

coef(fit.lasso_opt)

# sono un po' diversi se li estraggo con predict oppure se rifitto il modello con bestlam.lasso
# però poco, forse approssimazione?

fm = lm(y ~ computertime)
summary(fm)
fm$coefficients

# d)
5.7498287211 + 0.0002905625*100 + 0.3365513849*1 + 0.1231026862*10 + 0.0003359460*35



# Prova per CI and PI ma non va

Z0.new <- data.frame(age = 0, height = 0, dummy_sex = 0, 
                     distance = 100, siblings = 0.3365513849, computertime = 0.1231026862, 
                     exercisehours = 0, musiccds = 0.0003359460, playgames = 0 )


fit.lasso_opt <- glmnet(Z, y, lambda = bestlam.lasso)

coef(fit.lasso_opt)


# Conf. int. for the mean
alpha = 0.05
Conf <- predict(fit.lasso_opt, Z0.new, interval='confidence', level=1-alpha)  
Conf
# Pred. int. for a new obs
Pred <- predict(fit.lasso_opt, Z0.new, interval='prediction', level=1-alpha)  
Pred



detach(dataset)














