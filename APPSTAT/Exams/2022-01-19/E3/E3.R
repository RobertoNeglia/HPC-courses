setwd("~/shared-folder/HPC/APPSTAT/Exams/2022-01-19/E3")
library(car)


rm(list = ls())
tattoo <- read.csv("tattoo.txt", header = TRUE, sep = " ")
head(tattoo)
dim(tattoo)

attach(tattoo)
# linear model
dummy.method <- ifelse(method == "handmade", 1, 0)
lm <- lm(price ~ dimension + ncolors + dummy.method + dimension:dummy.method + ncolors:dummy.method)
summary(lm)
sum(residuals(lm)^2) / lm$df

coeffs <- lm$coefficients
coeffs
params <- rbind(handmade = c(coeffs[1] + coeffs[4], coeffs[2] + coeffs[5], coeffs[3] + coeffs[6]), machine = c(coeffs[1], coeffs[2], coeffs[3]))
print(params)

# question b)
# Verify assumptions (normality, homoscedasticity)
par(mfrow = c(2, 2))
plot(lm)

help(rstandard)
rstandard(lm)
# rstandard(lm) is the vector of standardized residuals, that is,
# residuals divided by their standard deviation.
shapiro.test(rstandard(lm))
C <- rbind(c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 1))
mu0 <- c(0, 0, 0)
linearHypothesis(lm, C, mu0)
# factor method certainly has an effect on price

C <- rbind(c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 1))
mu0 <- c(0, 0)
linearHypothesis(lm, C, mu0)
