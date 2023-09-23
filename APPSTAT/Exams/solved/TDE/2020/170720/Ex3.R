### Ex3

library(MASS)
library(car)
library(rgl)


dataset = read.table("toxicity.txt", header = T)

n = dim(dataset)[1]

attach(dataset)

fm <- lm(tox ~ C1 + C2 + C3 + C4 + C5 + C6)  # response ~ regressors

summary(fm) 

b = fm$coefficients
b

sum(residuals(fm)^2)/fm$df 

x11()
par(mfrow=c(2,2))
plot(fm)

shapiro.test(residuals(fm))

# b)
Z0.new <- data.frame(C1=100, C2=0.7 , C3 = 2, C4 = 4, C5 = 1.4 , C6=3)

alpha = 0.05
Pred <- predict(fm, Z0.new, interval='prediction', level=1-alpha)  
Pred

# c)
library(glmnet)

Z <- model.matrix(tox ~ C1 + C2 + C3 + C4 + C5 + C6)[,-1]
y = tox

## CV: Set lambda by cross-validation

# Let's set a grid of candidate lambda's for the estimate
lambda.grid <- 10^seq(0,-2,length=100)

fit.lasso <- glmnet(Z, y, lambda = lambda.grid) # default: alpha=1 -> lasso 

# Plot of the values of the coefficients
x11()
plot(fit.lasso , xvar='lambda', label=TRUE, col = rainbow(dim(Z)[2]))
legend('topright', dimnames(Z)[[2]], col =  rainbow(dim(Z)[2]), lty=1, cex=1)
# the coeff are going to 0 and at a certain point coeff reach 0 and then they
# remains 0 
# sometimes it's possible that we reach 0 then go up and then return to 0
# so if we choose a proper lambda we are also discarding some covariate

# CV
cv.lasso <- cv.glmnet(Z, y,lambda=lambda.grid) # default: 10-fold CV
# NB: if we have a low number of data we have to choose a low k for k-fold CV
# eg: n = 13 -> k=3 and we write in the argument before lambda: "nfolds=3"

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso

# Plot number of regressors that we keep
x11()
plot(cv.lasso)
abline(v=log(bestlam.lasso), lty=1)

coef.lasso <- predict(fit.lasso, s=bestlam.lasso, type = 'coefficients')[1: 7 ,]
coef.lasso 


# c) 

fit.lasso <- glmnet(Z, y, lambda = bestlam.lasso)

Z0.lasso = c(C1=100, C2=0.7 , C3 = 2, C4 = 4, C5 = 1.4 , C6=3)

alpha = 0.05
Pred <- predict(fit.lasso, Z0.lasso, interval='predicition', level=1-alpha)  
Pred


z0 = c(100, 0.7 ,2, 4,  1.4 , 3)
pred2 = predict(fit.lasso, s = bestlam.lasso, newx = z0) 
pred2





detach(dataset)



















