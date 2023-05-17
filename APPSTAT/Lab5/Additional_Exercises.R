setwd("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 6 - 29032021")
load("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 5 - 23042021/mcshapiro.test.RData")

library(car)

#_______________________________________________________________________________
##### Pb 1 of 10/02/10
# The file exchange.txt collects the daily exchange rates dollar/euro and
# pound/euro of Jan 2010. Assume that the 30 daily increments are independent
# realization of a bivariate Gaussian distribution.
# a) Is there statistical evidence to state that during January the exchange 
#    rate has changed in mean?
# b) Using the Bonferroni inequality, provide four confidence intervals of
#    global confidence 90% for the mean of the increments and for their 
#    variances.

rates <- read.table(file='exchange.txt', header=T)
head(rates)
dim(rates)

x11()
matplot(rates, type='l', ylim=c(0,3), lwd = 2, xlab='days')
legend(1,3,legend=c('Dollar','Pound'),col=1:2,lwd=2,lty=1:2)


### question a)

# we need to build the correct data matrix:
# we know that the increments are independent realizations of a
# bivariate Gaussian

diffrates <- matrix(NA, 30, 2)
for(i in 1:30)
  diffrates[i,] <- as.numeric(rates[i+1,] - rates[i,])

diffrates

# we assume that the sample is now iid from a bivariate normal distr.
# H0: rates did not change => mean of increments == zero

# we first need to verify the Gaussian assumption
plot(diffrates, asp=1, pch=1)
mcshapiro.test(diffrates)

mu0      <- c(0, 0)
x.mean   <- colMeans(diffrates)
x.cov    <- cov(diffrates)
x.invcov <- solve(x.cov)
n <- 30
p <- 2

x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
Pb <- 1-pf(x.T2*(n-p)/(p*(n-1)), p, n-p)
Pb

# mean under H0 (blue)
points(mu0[1], mu0[2], col='blue', pch=16)

# sample mean (black)
points(x.mean[1], x.mean[2], col='black', pch=16)

# we represent the confidence region of level 95%: where does mu0 fall?
alpha <- .05
cfr.fisher <- (p*(n-1)/(n-p))*qf(1-alpha,p,n-p)
ellipse(center=x.mean, shape=x.cov/n, radius=sqrt(cfr.fisher), lwd=2)

# what about the region of level 99%?
alpha <- .01
cfr.fisher <- (p*(n-1)/(n-p))*qf(1-alpha,p,n-p)
ellipse(center=x.mean, shape=x.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='orange', add=TRUE)
dev.off()

### question b): Bonferroni method
###           (intervals for the components of the mean AND on the variances)

k <- 4
alpha <- 0.1

ICmean <- cbind(inf=x.mean - sqrt(diag(x.cov)/n) * qt(1 - alpha/(2*k), n-1),
                center= x.mean,
                sup= x.mean + sqrt(diag(x.cov)/n) * qt(1 - alpha/(2*k), n-1))

ICvar <- cbind(inf=diag(x.cov)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
               center=diag(x.cov),
               sup=diag(x.cov)*(n-1) / qchisq(alpha/(2*k), n-1))

ICmean
ICvar


#_______________________________________________________________________________
##### Pb 3 of 10/02/10
# An association of Milanese technologists collected in the file mobile.txt
# the prices of iphone and Nokia 5800 recorded for 100 shops in the province
# of Milan
# a) Build three simultaneous T2 intervals with global confidence 90%
#    for the mean of the two prices and for their difference
# b) Nokia states that - at a worldwide level - Nokia 5800 costs in mean
#    one third of the price for an iphone. Based on the data, can we
#    deny this statement for the province of Milan?

mobile <- read.table('mobile.txt', header=T)
plot(mobile, asp=1)

dev.off()

mcshapiro.test(mobile)

### question a)
n <- dim(mobile)[1]
p <- dim(mobile)[2]
x.mean   <- sapply(mobile,mean)
x.cov    <- cov(mobile)

alpha <- 0.10
cfr.fisher <- (p*(n-1)/(n-p))*qf(1-alpha,p,n-p)

A <- rbind(c(1,0), c(0,1), c(-1,1))

ICT2 <- cbind(A %*% x.mean - sqrt(diag(A %*% x.cov %*% t(A))/n * cfr.fisher), 
              A %*% x.mean,
              A %*% x.mean + sqrt(diag(A %*% x.cov %*% t(A))/n * cfr.fisher)) 
ICT2

plot(mobile, asp=1)
ellipse(x.mean, x.cov/n, cfr.fisher, add=T)

### question b)

# Test:
# H0: mu_nokia == mu_iphone/3 vs H1: mu_nokia != mu_iphone/3
# i.e.
# H0: mu_iphone - 3* mu_nokia == 0 vs H1: mu_iphone - 3* mu_nokia != 0
# (test a particular linear combination)

M   <- mobile[,2] - 3*mobile[,1]
t.M <- abs(mean(M))/sqrt(var(M)/n)
t.M
P   <- (1 - pt(t.M, n-1))*2
P

# or: (it is exactly analogous)
a <- c(-3, 1)
t.M.bis <- abs(x.mean%*%a)/sqrt((t(a)%*%x.cov%*%a)/n)
t.M.bis
P   <- (1 - pt(t.M.bis, n-1))*2
P

#_______________________________________________________________________________
##### Pb 3 of 24/09/10
# The PoliPharma has given a drug to 10 students in order to increase the levels
# of Matina and decrease the level of Fisina in their blood. In the files before.txt
# and after.txt the levels of Matina, Fisina, Chimina and Elettrina are reported
# for 10 students, before and after the administration of the drug.
# a) Having introduced and tested the appropriate hypotheses of Gaussianity,
#    perform a test of level 1% to prove the existence of an effect of the 
#    drug on the mean levels of the four enzymes.
# b) Provide four simultaneous T2 intervals for the mean of the four
#    increments
# c) Perform a test of level 1% to confirm/deny what stated by the Polipharma,
#    that is: the ingestion of the drug causes a mean increment of 2 units of
#    the Matina, a decrease of 1 unit in the Fisina and a mean increment in the
#    Chimina equal to the mean decrease in the Elettrina

before <- read.table('before.txt', header=T)
after <- read.table('after.txt', header=T)

### question a)
D <- after - before

# verify normality
mcshapiro.test(D)

# test
alpha <- 0.01
n <- dim(D)[1]
p <- dim(D)[2]

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

delta.0 <- c(0,0,0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)

cfr.fisher <- (p*(n-1)/(n-p))*qf(1-alpha,p,n-p)
D.T2 < cfr.fisher

# pvalue (not requested)
P <- 1-pf(D.T2 * (n-p)/(p*(n-1)), p, n-p)
P


### question b)

# Simultaneous T2 intervals for the components
T2 <- cbind( "Inf"= D.mean-sqrt(cfr.fisher*diag(D.cov)/n) , D.mean, Sup=D.mean+sqrt(cfr.fisher*diag(D.cov)/n) )
T2

### question c)

# we want to perform a test of global level 1% for the three linear 
# combination of the mean of D (datasets of the differences)
# H0: delta1 = 2, delta2 = -1, delta3 = - delta4
# i.e.,
# H0: delta1 = 2, delta2 = -1, delta3 + delta4 = 0

A <- rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,1))

## Way 1: we build a new dataset
D.new <- data.frame(as.matrix(D) %*% t(A))

p.new <- dim(D.new)[2]
n <- dim(D.new)[1]

alpha <- 0.01

D.new.mean   <- sapply(D.new,mean)
D.new.cov    <-  cov(D.new)
D.new.invcov <- solve(D.new.cov)

delta.0 <- c(2,-1,0)

D.T2 <- n * (D.new.mean-delta.0) %*% D.new.invcov %*% (D.new.mean-delta.0)
D.T2

cfr.fisher.new <- (p.new*(n-1)/(n-p.new))*qf(1-alpha,p.new,n-p.new)
D.T2 < cfr.fisher.new

# pvalue (not requested)
P <- 1-pf(D.T2 * (n-p.new)/(p.new*(n-1)), p.new, n-p.new)
P

## Way 2 (totally analogous): 
## We perform a test for 
## H0: A%*%mu.D == delta.0 vs H1: A%*%mu.D != delta.0
A <- rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,1))
delta.0 <- c(2,-1,0)

T2.A <- n * t(A %*% D.mean - delta.0) %*% solve(A %*% D.cov %*% t(A)) %*% (A %*% D.mean - delta.0)
T2.A

T2.A < cfr.fisher.new
# WARNING: here, the dimension 'p' that we use for cfr.fisher 
#          is p.new=3, it isn't the initial dimension of the data
#          (that was p=4) [ if I kept p=4 I would have obtained a
#          global level lower than 1% ]

#_______________________________________________________________________________
##### Pb 2 of 10/09/10
# The file pound.txt collects the exchange rates pound/euro used by 24 banks in 
# UK during the week between 30th Aug and 5th Sept 2010.
# a) Having framed the problem in the context of repeated measures, perform
#    a test of level 5% to verify the hypothesis that the mean of the exchange
#    rate is constant along time
# b) Provide 6 confidence intervals of global level 95% for the daily increments
#    in the exchange rate
# c) Comment the 6 intervals computed at point b).

pound <- read.table('pound.txt', header=T)
head(pound)

matplot(t(pound), type='l')

### question a)
n <- dim(pound)[1]
q <- dim(pound)[2]

M <- sapply(pound,mean)
S <- cov(pound)

# matrix of contrasts that looks at the daily increments
C   <-  matrix(c(-1, 1, 0, 0, 0, 0, 0,
                 0, -1, 1, 0, 0, 0, 0,
                 0, 0, -1, 1, 0, 0, 0,
                 0, 0, 0, -1, 1, 0, 0,
                 0, 0, 0, 0, -1, 1, 0,
                 0, 0, 0, 0, 0, -1, 1), 6, 7, byrow=T)
C
# I could choose any contrast matrix but if I look at the following
# questions I can save energy and time!
mcshapiro.test(pound)

# Test: H0: C%*%mu=0 vs H1: C%*%mu!=0
alpha   <- .05
delta.0 <- c(0, 0, 0, 0, 0, 0)

Md <- C %*% M 
Sd <- C %*% S %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) 
# attention to the multiplying factor!
T2 < cfr.fisher

P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1))
P

### questions b) and c)
k <- q-1
ICmean <- cbind(Md - sqrt(diag(Sd)/n) * qt(1 - alpha/(2*k), n-1),
                 Md,
                 Md + sqrt(diag(Sd)/n) * qt(1 - alpha/(2*k), n-1))
ICmean

### questions b) e c)
for (i in 1:(q-1))
  print(paste('Reject H0 in direction ',i,': ', !(delta.0[i]>ICmean[i,1] & delta.0[i]<ICmean[i,3]),sep=''))

#_______________________________________________________________________________
##### Pb 3 of 10/09/10
# The file extra.txt reports the representation expenses [$] of the English
# first minister and of his vice during the first 12 months of 2009. Assume
# those data to be independent realizations of a bivariate Gaussian.
# a) Build an ellipsoidal region of confidence 90% for the mean of the 
#    representation expenses
# b) Is there evidence of the fact that the prime minister spends in mean
#    more than twice the expenses of its vice?
# c) build a confidence interval of level 90% for the mean of the sum of
#    the expenses.

extra <- read.table('extra.txt', header=T)

plot(extra, asp=1)
extra

### question a)
mcshapiro.test(extra)

n <- dim(extra)[1]
p <- dim(extra)[2]

x.mean <- sapply(extra,mean)
x.cov <- cov(extra)
x.inv <- solve(x.cov)

cfr.fisher <- (n-1)*p/(n-p)*qf(1-alpha,p,n-p)

ellipse(center=x.mean, shape=x.cov/n, radius=sqrt(cfr.fisher), lwd=2)
# centre:
x.mean
# direction of the axes
eigen(x.cov)$vector[,1]
eigen(x.cov)$vector[,2]
# radius
sqrt(cfr.fisher)
# lenght of the semi-axes
sqrt(eigen(x.cov/n)$values)*sqrt(cfr.fisher)

### question b)
# Test:   H0: mu1<=2*mu2   vs H1: mu1>2*mu2
# i.e. H0: mu1-2*mu2<=0 vs H1: mu1-2*mu2>0
# i.e H0: a'mu<=0      vs H1: a'mu>0 con a=c(1,-2)

a <- c(1,-2)
delta.0 <- 0

extra <- as.matrix(extra)
t.stat <- (mean(extra %*% a) - delta.0 ) / sqrt( var(extra %*% a) / n ) # t-statistics (statistica 1!)

# Reject for large values of t 
# => compute the p-value as the probability of the right tail (i.e., of values >tstat)
P <- 1-pt(t.stat, n-1)
P

### question c)
a2 <- c(1,1)
alpha <- .1
cfr.t <- qt(1-alpha/2, n-1)

c(inf = mean(extra %*% a2) - cfr.t * sqrt( var(extra %*% a2) / n ),
  center = mean(extra %*% a2),
  sup = mean(extra %*% a2) + cfr.t * sqrt( var(extra %*% a2) / n ))

# Otherwise one can use the function t.test()
lc <- extra[,1] + extra[,2]
t.test(lc, alternative = 'two.sided', mu = 0, conf.level = 0.90)


#_______________________________________________________________________________
##### Pb 1 of 29/06/11
# The Spanish Authority for the paella measured the content of clams (vongole) [g] and 
# shrimps (gamberetti) [g] in a number of packs of "Paella of Cantabro" of 200g 
# (file cantabro.txt). They aim to verify if the mean content of clams and shrimps
# isn't significantly different from the nominal one: 30g of clams and 50g 
# of shrimps.
# a) Perform a test to prove the former hypothesis.
# b) Comment the test at point a) using three simultaneous T2 intervals (global
#    level 90%) for the mean content of clams, of shrimps and of their sum.
# c) Do you deem necessary an action of the Authority? Motivate the response.

cantabro <- read.table('cantabro.txt', header=T)
head(cantabro)

### question a)
mcshapiro.test(cantabro)

p <- dim(cantabro)[2]
n <- dim(cantabro)[1]

# Test: H0: mu=c(30,50) vs H1: mu!=c(30,50)
alpha    <- 0.01
mu0      <- c(30,50)
x.mean   <- colMeans(cantabro)
x.cov    <- cov(cantabro)
x.invcov <- solve(x.cov)

x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
x.T2 < cfr.fisher

P <- 1-pf(x.T2*(n-p)/(p*(n-1)), p, n-p)
P

### question b)
a1<-c(1,0)
a2<-c(0,1)
a3<-c(1,1)

ICT21<-data.frame(L=t(a1)%*%x.mean-sqrt(t(a1)%*%x.cov%*%a1/n)*sqrt(cfr.fisher),C=t(a1)%*%x.mean,U=t(a1)%*%x.mean+sqrt(t(a1)%*%x.cov%*%a1/n)*sqrt(cfr.fisher))
ICT22<-data.frame(L=t(a2)%*%x.mean-sqrt(t(a2)%*%x.cov%*%a2/n)*sqrt(cfr.fisher),C=t(a2)%*%x.mean,U=t(a2)%*%x.mean+sqrt(t(a2)%*%x.cov%*%a2/n)*sqrt(cfr.fisher))
ICT23<-data.frame(L=t(a3)%*%x.mean-sqrt(t(a3)%*%x.cov%*%a3/n)*sqrt(cfr.fisher),C=t(a3)%*%x.mean,U=t(a3)%*%x.mean+sqrt(t(a3)%*%x.cov%*%a3/n)*sqrt(cfr.fisher))
ICT2<-data.frame(rbind(ICT21,ICT22,ICT23))
ICT2

### question c)
h0 <- c(30,50,30+50)
for (i in 1:3)
  print(paste('Reject H0 for a',i,': ', !(h0[i]>ICT2[i,1] & h0[i]<ICT2[i,3]),sep=''))

#_______________________________________________________________________________
##### Pb 2 of 18/07/11
# The file uranium.txt reports the quotations of uranium [$/kg] in 32 
# stock exchange, for the days between 1st and 10th July 2011. 
# Having frames the problem in the context of repeated measures,
# answer the following questions.
# a) Is there evidence at 90% that the mean price was not constant
#    along the 10 days?
# b) Build 9 Bonferroni confidence intervals (global confidence 90%) for the
#    increments/decrements of the mean daily prices.
# c) Comment the conclusions deduced from the analyses (a) and (b)
# d) The IMF states that between the 8th and 9th July, the prices had a 
#    mean decrease of 1$/kg. Perform a test to confirm/deny this statement.

uranium <- read.table('uranium.txt', header=T)
uranium 

### question a) (with a look to b)!)
C <- matrix (0, nrow=9, ncol=10)
for(i in 1:9)
  C[i,c(i,i+1)]<-c(-1,1)

mcshapiro.test(uranium)

n <- dim(uranium)[1]
q <- dim(uranium)[2]

alpha    <- 0.10
delta.0      <- rep(0,9)
x.mean   <- colMeans(uranium)
x.cov    <- cov(uranium)
x.invcov <- solve(x.cov)

Md <- C %*% x.mean 
Sd <- C %*% x.cov %*% t(C)
Sdinv <- solve(Sd)

T2 <- n * t( Md - delta.0 ) %*% Sdinv %*% ( Md - delta.0 )

cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1))
T2 < cfr.fisher

P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)),(q-1),n-(q-1))
P

### question b)
alpha <- 0.10
k<-q-1
cfr.t <- qt(1-alpha/(2*k),n-1)

ICB <- cbind(L=Md-cfr.t*sqrt(diag(Sd)/n),C=Md,U=Md+cfr.t*sqrt(diag(Sd)/n))
ICB

### question c)
for (i in 1:k)
  print(paste('Reject H0 for the increment ',i,': ', !(0>ICB[i,1] & 0<ICB[i,3]),sep=''))

### question d)
a <- c(0,0,0,0,0,0,0,-1,1,0)  # combination that gives the increment between day9 and day8
a

# Test: a'mu==-1 vs a'mu!=-1  (decrease of 1 == increment of -1)
alpha <- 0.10
delta.0 <- -1

uranium <- as.matrix(uranium)
t2 <- (mean(uranium %*% a) - delta.0 )^2 / ( var(uranium %*% a) / n ) 
cfr.fisher <- qf(1-alpha,1,n-1)
t2<cfr.fisher

P <- 1-pf(t2,1,n-1)
P

## N.B. this test is univariate!
##      it is totally equivalent to:
t.test(uranium%*%a, alternative = 'two.sided', mu = delta.0, conf.level = 1-alpha)
t2
(-1.3699)^2

#_______________________________________________________________________________
##### Pb 1 of 4/07/2004
#####--------------------
# A lab is conducing a clinical trial to test the effectiveness of a new drug
# for diabetes. To each patient, the following quantities are measured 2h before
# and 2h after the administration of the drug: glycemia, body temperature,
# min and max blood pressure. The data are collected in the file diabete.txt.
# The pharmaceutical company that produces the drug declares that the drug
# is able to decrease of 60 units the glycemia, without any side effect on
# the body temperature, min and max pressure.
# a) perform a test to verify the statement of the producer.
# b) Using Bonferroni intervals of global confidence 95%, analyse the results
#    of the test at point (a).
# c) Verify the assumptions needed to execute the test at point (a).

diabetes <- read.table('diabetes.txt', header=TRUE)
diabetes

# question a)
D <- diabetes[,5:8]-diabetes[,1:4]
names(D) <- c('glycemia_diff', 'temperature_diff',
              'pressure_min_diff', 'pressure_max_diff')

# Test: H0: mu_D=c(-60,0,0,0) vs H1: mu_D!=0
n <- dim(D)[1]
p <- dim(D)[2]
M <- sapply(D,mean)
S <- cov(D)
Sinv <- solve(S)

delta0 <- c(-60,0,0,0)

T2 <- n*t(M-delta0)%*%Sinv%*%(M-delta0)

pvalue <- 1-pf(T2*(n-p)/(p*(n-1)), p, n-p)
pvalue

# question b)
alpha <- .05
k <- 4

cfr.t <- qt(1-alpha/(2*k),n-1)

ICB<-data.frame(L=M-sqrt(diag(S)/n)*cfr.t, C=M, U=M+sqrt(diag(S)/n)*cfr.t)
ICB

# question c)
mcshapiro.test(D)




#_______________________________________________________________________________
##### Pb 1 of 26/02/2008
#####--------------------
# Let X=(X1 X2 X3)'~N(mu, Sigma) a Gaussian random vector with
# mu=(1 1 1)' and Sigma=cbind(c(5,3,1),c(3,5,1),c(1,1,1)).
# a) Identify a region A such that P((X1 X2)' \in A)=0.9
# b) Identify a region A2 such that P((X1 X2)' \in A2 | X3=1)=0.9
# c) Having reported in a plot the graphs of the two regions, order in increasing
#    order the following probabilities:
#    P((X1 X2)' \in A )
#    P((X1 X2)' \in A2)
#    P((X1 X2)' \in A  | X3=1)
#    P((X1 X2)' \in A2 | X3=1)

mu=c(1,1,1)
Sigma=cbind(c(5,3,1),c(3,5,1),c(1,1,1))

### a) Consider only (X1 X2)'
eigen(Sigma[1:2,1:2])

# Direction of the axes:
eigen(Sigma[1:2,1:2])$vectors

# Center:
M <- mu[1:2]

# Radius of the ellipse:
r <- sqrt(qchisq(0.9,2))

# Length of the semi-axes:
r*sqrt(eigen(Sigma[1:2,1:2])$values)

# Plot
plot(M[1],M[2],xlim=c(-10,15),ylim=c(-10,15),col='blue',pch=19,xlab='X.1',ylab='X.2',asp=1)
ellipse(center=M, shape=cbind(Sigma[1:2,1:2]), radius=r, col = 'blue')
abline(h=0, v=0, lty=2, col='grey')
abline(a=0,b=1,col='grey',lty=2)
abline(a=2,b=-1,col='grey',lty=2)

### b) Consider the conditional distribution (X1 X2)'|X3=1

# Functions to compute the mean and the covariance matrix of the conditional
# distribution
mu.cond <- function(mu1,mu2,Sig11,Sig12,Sig22,x2)
{
  return(mu1+Sig12%*%solve(Sig22)%*%(x2-mu2))
}

Sig.cond <- function(Sig11,Sig12,Sig22)
{
  Sig21=t(Sig12)
  return(Sig11-Sig12%*%solve(Sig22)%*%Sig21)
}

M.c <- mu.cond(mu1=mu[1:2],mu2=mu[3],Sig11=Sigma[1:2,1:2],Sig12=Sigma[1:2,3],Sig22=Sigma[3,3],x2=1)
M.c
Sigma.c <- Sig.cond(Sig11=Sigma[1:2,1:2],Sig12=Sigma[1:2,3],Sig22=Sigma[3,3])
Sigma.c

eigen(Sigma.c[1:2,1:2])

# Direction of the axes:
eigen(Sigma.c[1:2,1:2])$vectors

# Center:
M.c

# Radius of the ellipse:
r <- sqrt(qchisq(0.9,2))

# Length of the semi-axes:
r*sqrt(eigen(Sigma.c[1:2,1:2])$values)

# Plot
ellipse(center=as.vector(M.c), shape=Sigma.c, radius=r, col = 'red')
