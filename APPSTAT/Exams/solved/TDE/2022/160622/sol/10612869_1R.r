# The files discomaniac.txt and lipsticks.txt contain the prices (in e) and the conditions of sales (a real number
# from 0 to 12) of the same 20 LPs sold from two stores.
# a) Perform a statistical test of level 95% to verify if the mean prices and conditions of the two stores differ.
# b) Which are the assumptions of the previous test? Are they met?
#   c) Provide the plot of the confidence region at level 95% for the differences of the mean of prices and the conditions
# between the two stores.
# d) Provide four Bonferroni simultaneous confidence intervals of global level 95% for the means and the variances
# of the differences in prices and conditions. Interpret the results of the test at point (a) through the intervals.

load("~/Desktop/Applied Statistics/Labs/Lab 5/mcshapiro.test.RData")

discomaniac <- read.table('discomaniac.txt', header=T)
lipsticks <- read.table('lipsticks.txt', header=T)

library(car)


# a) Perform a statistical test of level 95% to verify if the mean prices and conditions of the two stores differ.

#Test on the means of two groups (assumptions: normality, independence and same covariance)
#the test is: H0: mu.1=mu.2 vs H1: mu.1 != mu.2
M1 <- colMeans(discomaniac[,3:4])
M2 <- colMeans(lipsticks[,3:4])
n1 <- 20
n2 <- 20
p <- 2
S1 <- cov(discomaniac[,3:4])
S2 <- cov(lipsticks[,3:4])
Sp <- ( (n1-1)*S1 + (n2-1)*S2 ) / (n1+n2-2)
Sp.inv <- solve(Sp)
alpha <- 0.05
cfr.fisher <- (n1+n2-2)*p/(n1+n2-1-p) * qf(1-alpha,p,n1+n2-2)   
T2 <- (1/n1 + 1/n2)^-1 * (M1-M2) %*% Sp.inv %*% (M1-M2)  
T2 < cfr.fisher #we accept null hypothesis


# b) Which are the assumptions of the previous test? Are they met?

mcshapiro.test(discomaniac[,3:4])$p
mcshapiro.test(lipsticks[,3:4])$p

par(mfrow=c(2,1))
image(S1)
image(S2)


# c) Provide the plot of the confidence region at level 95% for the differences of the mean of prices and the conditions
# between the two stores.

data <- discomaniac[,3:4] - lipsticks[,3:4]


plot(data,xlim=c(-2,3.5),xlab='Price difference',ylab='Condition difference')
ellipse(center=M1-M2, shape=Sp/n1, radius=sqrt(cfr.fisher), lwd=2)

# d) Provide four Bonferroni simultaneous confidence intervals of global level 95% for the means and the variances
# of the differences in prices and conditions. Interpret the results of the test at point (a) through the intervals.

k <- 4
alpha <- 0.05

A  <- rbind(c(1,0), c(0,1)) 
k  <- dim(A)[1]
A.s2 <- diag(A%*%Sp%*%t(A))
A.dm <- A%*%(M1-M2)
Bonf <- cbind(inf=A.dm - qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ), 
              center=A.dm, 
              sup=A.dm + qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ))
Bonf

ICVar1 <- cbind(inf= Sp[1,1]*(n1+n2-2)/qchisq(1-alpha/(2*k),n1+n2-2),
               center= Sp[1,1],
               sup= Sp[1,1]*(n1+n2-2)/qchisq(alpha/(2*k),n1+n2-2))
ICVar1
#price: 14.12882 22.36638 40.1193


ICVar2 <- cbind(inf= Sp[2,2]*(n1+n2-2)/qchisq(1-alpha/(2*k),n1+n2-2),
                center= Sp[2,2],
                sup= Sp[2,2]*(n1+n2-2)/qchisq(alpha/(2*k),n1+n2-2))
ICVar2
#condition: 1.901421 3.01001 5.399152
