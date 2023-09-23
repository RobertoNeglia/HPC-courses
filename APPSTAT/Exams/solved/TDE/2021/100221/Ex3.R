### Ex3


dataset = read.table("landslides.txt", header = T)

n          <- dim(dataset)[[1]]
y          <- dataset$rate         # target
z1         <- dataset$rain        # regressor z1  
z2         <- dataset$hardness     # regressor z2
z3 = dataset$coarse
z4 = dataset$fine

fm <- lm(y ~ z1 + z2 + z3 + z4)
summary(fm)

b = fm$coefficients


# b)

fm1 <- lm(y ~ z1 + z3 + z4)
summary(fm1)

b = fm1$coefficients

# c)

library(MASS)
library(car)
library(rgl)

C <- rbind(c(0,0,1,-2))
p = dim(C)[1]
linearHypothesis(fm1, C, 0)


# CONSTRAINED MODEL  -> VEDI WORD


# fm2 <- lm(y ~ z1 + z3)
# summary(fm2)
# 
# b_new = fm2$coefficients

# fm3 <- lm(y ~ z1 + z4)
# summary(fm3)
# # worse 



# d)  (calcoli fermi alla sol di prima)
Z0.new <- data.frame(z1=700, z5=10 + 0.5*8)

# Conf. int. for the mean
alpha = 0.01
Conf <- predict(fm2, Z0.new, interval='confidence', level=1-alpha)  
Conf




