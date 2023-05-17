###--------------------###
### LAB 5 (24/03/2023) ###
###--------------------###

### TOPICS:
### Test for the mean of paired multivariate Gaussian observations
### Test for repeated measures
### Test for two independent Gaussian populations
### Exercises

library(car)
library(MVN)
quartz.options(height=6.5, width=11, reset=FALSE)

#_______________________________________________________________________________
##### Test for the mean of paired multivariate Gaussian observations
#####----------------------------------------------------------------

###   Example 6.1 JW + exercise 6.1 JW
###   (n=11, p=2)
###------------------------------------
# Municipal wastewater treatment plants are required by law to monitor
# their discharges into rivers and streams on a regular basis. Concern
# about the reliability of data from one of these self-monitoring 
# programs led to a study in which samples of effluent were divided 
# and sent to two laboratories for testing. One-half of each sample
# was sent to the Wisconsin State Laboratory of Hygiene, and one-half 
# was sent to a private commercial laboratory routinely used in the 
# monitoring program. Measurements of biochemical oxygen demand (BOD)
# and suspended solids (SS) were obtained, for n=11 sample splits, 
# from the two laboratories.
# (The experimenter divided each sample by first shaking it and then 
# pouring it rapidly into two bottles in order to avoid difference 
# in the suspended solids contained in the two half-samples)

# Do the two laboratories' chemical analyses agree?

effluent <- read.table('effluent.dat')
effluent

# We change the column names to have the correct label
colnames(effluent) <- c('BOD_Lab1', 'SS_Lab1', 'BOD_Lab2', 'SS_Lab2')

quartz()
pairs(effluent, pch=19, main='Dataset effluent')

dev.off()

# we compute the sample of differences
D <- data.frame(
  DBOD = effluent$BOD_Lab1 - effluent$BOD_Lab2,
  DSS  = effluent$SS_Lab1 - effluent$SS_Lab2
) 
D

# Scatter plot of the dataset of differences with the four quadrants
quartz()
plot(D, asp=1, pch=19, main='Dataset of Differences')
abline(h=0, v=0, col='grey35')
points(0,0, pch=19, col='grey35')
# DBOD: difference in 'biochemical oxygen demand'
# DSS:  difference in 'suspended solids'
#       measured by the two laboratories

dev.off()

# Now we can proceed as we already know, but working on D

### T2 Hotelling Test 
# H0: delta == delta.0 vs H1: delta != delta.0
# with delta.0 = c(0,0)

# Test the Gaussian assumption (on D!)
result <- mvn(data = D)
result$multivariateNormality

# The p-value isn't very high (but I don't reject for levels 5%, 1%).
# There might be outliers, but we don't remove them because we only have
# very few data

n <- dim(D)[1]  # 11
p <- dim(D)[2]  #  2

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

alpha   <- .05
delta.0 <- c(0,0)

D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)
D.T2

cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher

D.T2 < cfr.fisher # FALSE: we reject H0 at level 5%

# we compute the p-value
P <- 1 - pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P
# reject H0 at 5% (don't reject at 1%)

# Ellipsoidal confidence region with confidence level 95%
quartz()
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2)

# Adding delta.0 and the quadrants
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='grey35')

# Ellipsoidal confidence region with confidence level 99%
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt((n-1)*p/(n-p)*qf(1-0.01,p,n-p)),
        lty=2,col='grey',lwd=2)

# What if we set the radius as the quantile of order 1-pval?
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt((n-1)*p/(n-p)*qf(1-as.numeric(P),p,n-p)),
        lty=1, col='dark grey', lwd=2)

dev.off()

# Now, let's communicate our results to the client (because we need to take our salary!)
# Let's build confidence intervals for linear combination of the
# components of the mean vector

### Simultaneous T2 intervals in the direction of DBOD and DSS
IC.T2.DBOD <- c( D.mean[1] - sqrt(cfr.fisher*D.cov[1,1]/n), D.mean[1],
                 D.mean[1] + sqrt(cfr.fisher*D.cov[1,1]/n))
IC.T2.DSS  <- c( D.mean[2] - sqrt(cfr.fisher*D.cov[2,2]/n), D.mean[2],
                 D.mean[2] + sqrt(cfr.fisher*D.cov[2,2]/n))

T2 <- rbind(IC.T2.DBOD, IC.T2.DSS)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2

# Plot of the simulatenous T2 intervals
quartz()
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey')
abline(v = T2[1,1], col='red', lwd=1, lty=2)
abline(v = T2[1,3], col='red', lwd=1, lty=2)
abline(h = T2[2,1], col='red', lwd=1, lty=2)
abline(h = T2[2,3], col='red', lwd=1, lty=2)

points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.5)
abline(h=delta.0[1], v=delta.0[2], col='grey35')

segments(IC.T2.DBOD[1], 0, IC.T2.DBOD[3], 0, lty=1, lwd=2, col='red')
segments(0, IC.T2.DSS[1], 0, IC.T2.DSS[3], lty=1,lwd=2, col='red')

# We don't have enough evidence to reject H0 in any coordinate axes-direction
# But we can reject the global test of (global) level 5%
# => there exists at least one direction along which we are allowed
#    to reject the univariate test

# Recall:
# We reject the global H0 if in at least one direction we observe a 'high'
# value of the statistics T2 (univariate), i.e., we reject the global
# H0 if we reject the univariate test at least in direction (max(T2))

# Hence:
# Worst direction: direction along which the T2 statistics (univariate)
# is maximized

# => From the theory
#    - the maximum is realized (Hotelling T2-statistics)
D.T2
#    - the distribution of the maximum is known
#    - the direction along which the maximum is realized is known
worst <- D.invcov %*% (D.mean-delta.0)
worst <- worst/sqrt(sum(worst^2)) # Normalization
worst

# Confidence interval along the worst direction:
IC.worst  <- c(D.mean %*% worst - sqrt(cfr.fisher*(t(worst) %*% D.cov %*% worst) / n),
               D.mean %*% worst,
               D.mean %*% worst + sqrt(cfr.fisher*(t(worst) %*% D.cov %*% worst) / n))
IC.worst
delta.0 %*% worst
(IC.worst[1] < delta.0 %*% worst) & (delta.0 %*% worst < IC.worst[3])   
# Reject H0: a'mu == a'delta.0 in direction a = worst

# Extremes of IC.worst in the coordinate system (x,y):
x.min <- IC.worst[1] * worst # (x,y) coords of the lower bound of the interval
x.max <- IC.worst[3] * worst # (x,y) coords of the upper bound of the interval
m1.ort <- - worst[1] / worst[2] # Slope of the line orthogonal to worst
q.min.ort <- x.min[2] - m1.ort * x.min[1] # Intercept of line of slope m1.ort
# passing by x.min
q.max.ort <- x.max[2] - m1.ort * x.max[1] # Intercept of line of slope m1.ort
# passing by x.max
abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)

m1 = worst[2] / worst[1] # worst direction
abline(0, m1, col='grey35') # Intercept 0 because the line has to pass
# by delta.0 which is (0, 0)
segments(x.min[1], x.min[2], x.max[1], x.max[2], lty=1, lwd=2, col='forestgreen')

dev.off()

# If we are not convinced yet, let's look at all the directions:
# we compute confidence intervals for a'x where a varies in all the 
# directions between 0 and pi, with step pi/180. For each direction
# we compute the T2 statistic (univariate)
D <- as.matrix(D)
theta   <- seq(0, pi - pi/180, by = pi/180)
T2.d     <- NULL
Centerf  <- NULL
Maxf     <- NULL
Minf     <- NULL

for(i in 1:length(theta))
{
  a   <- c(cos(theta[i]), sin(theta[i]))
  t2  <- ( mean(D %*% a) - (delta.0 %*% a) )^2 / ( var(D %*% a) / n )
  T2.d  <- c(T2.d, t2)
  centerf  <- D.mean %*% a
  maxf     <- D.mean %*% a + sqrt( t(a) %*% D.cov%*% a / n) * sqrt(cfr.fisher)
  minf     <- D.mean %*% a - sqrt( t(a) %*% D.cov%*% a / n) * sqrt(cfr.fisher)
  Centerf  <- c(Centerf, centerf)
  Maxf     <- c(Maxf, maxf)
  Minf     <- c(Minf, minf)
}

quartz()
par(mfrow=c(1,3))

# Scatter plot of DSS and DBOD
plot(D, asp=1, pch=1, main='Dataset of the Differences', ylim=c(-15,60))

# Adding the quadrants
abline(h=delta.0[1], v=delta.0[2], col='red', lty=3)

# Adding the 95% confidence region for the true mean
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey')

# Adding the simultaneous confidence intervals in the direction of DSS and DBOD
segments(IC.T2.DBOD[1],0,IC.T2.DBOD[3],0,lty=1,lwd=2,col='red')
segments(0,IC.T2.DSS[1],0,IC.T2.DSS[3],lty=1,lwd=2,col='red')

# Adding the simultaneous confidence interval in the "worst" direction
x.min <- IC.worst[1] * worst
x.max <- IC.worst[3] * worst
segments(x.min[1], x.min[2], x.max[1], x.max[2], lty=1, lwd=2, col='forestgreen')
abline(0, m1, col='forestgreen', lty=3)
points(delta.0[1], delta.0[2], pch=16, col='black')

# Second plot with the confidence intervals for each direction tested
plot(theta, Centerf, main = 'Simultaneous T2 confidence intervals', ylim = c(-30,35), col = 'grey25', type='l',ylab='IC')
for(i in 1:length(theta))
  lines(c(theta[i], theta[i]), c(Minf[i], Maxf[i]), col = 'grey75')
abline(h=0, col='black')
lines(theta, Minf, col = 'red', lty = 2)
lines(theta, Maxf, col = 'red', lty = 2)

# Highlighting the directions of DSS and DBOD (theta = 0 and theta = pi/2)
lines(c(theta[1], theta[1]), c(Minf[1], Maxf[1]), col = 'red', lwd=2)
lines(c(theta[91], theta[91]), c(Minf[91], Maxf[91]), col = 'red', lwd=2)

# Highlighting the "worst" direction for the rejection of H0
lines(c(theta[which.max(T2.d)], theta[which.max(T2.d)]),
      c(Minf[which.max(T2.d)], Maxf[which.max(T2.d)]), col = 'forestgreen', lwd=2)

# Third plot with the value of T2 statistic
plot(theta, T2.d, main = 'T2 statistics', ylim = c(0,15), col = 'blue', type='l')
abline(v=c(0,pi/2), col = 'red', lty = 3)
abline(v=theta[which.max(T2.d)], col = 'forestgreen', lty = 3)

# Drawing the threshold given by the Fisher quantile at level .95
abline(h=cfr.fisher, col = 'grey', lty = 1, lwd=2)

dev.off()

### Bonferroni intervals
k <- p  # Number of Bonferroni intervals
cfr.t <- qt(1-alpha/(2*k), n-1) # Student quantile

# Bonferroni confidence intervals in the direction of DBOD and DSS
IC.BF.DBOD <- c(D.mean[1] - cfr.t*sqrt(D.cov[1,1]/n),
                D.mean[1],
                D.mean[1] + cfr.t*sqrt(D.cov[1,1]/n))
IC.BF.DSS  <- c(D.mean[2] - cfr.t*sqrt(D.cov[2,2]/n),
                D.mean[2],
                D.mean[2] + cfr.t*sqrt(D.cov[2,2]/n))
# Bonferroni region defined by the cartesian product of the Bf intervals
Bf <- rbind(IC.BF.DBOD, IC.BF.DSS)
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf

quartz()
# Scatter plot of DSS and DBOD
plot(D, asp=1, pch=1, main='Dataset of the Differences',ylim=c(-15,60))

# Adding the 95% confidence region for the true mean of the differences
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=2, col='grey',
        center.cex=1.25)

# Adding quadrant lines and delta.0
abline(h=0, v=0, col='grey', lty=1, lwd=2)
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

# Simultaneous confidence intervals in the direction of DSS and DBOD
abline(v = T2[1,1], col='red', lwd=1, lty=2)
abline(v = T2[1,3], col='red', lwd=1, lty=2)
abline(h = T2[2,1], col='red', lwd=1, lty=2)
abline(h = T2[2,3], col='red', lwd=1, lty=2)
segments(IC.T2.DBOD[1], 0, IC.T2.DBOD[3], 0, lty=1, lwd=2, col='red')
segments(0, IC.T2.DSS[1], 0, IC.T2.DSS[3], lty=1, lwd=2, col='red')

# Bonferroni intervals in the direction of DSS and DBOD
abline(v = Bf[1,1], col='blue', lwd=1, lty=2)
abline(v = Bf[1,3], col='blue', lwd=1, lty=2)
abline(h = Bf[2,1], col='blue', lwd=1, lty=2)
abline(h = Bf[2,3], col='blue', lwd=1, lty=2)
segments(IC.BF.DBOD[1], 0, IC.BF.DBOD[3], 0, lty=1, lwd=2, col='blue')
segments(0, IC.BF.DSS[1], 0, IC.BF.DSS[3], lty=1, lwd=2, col='blue')


dev.off()
# Note: The simultaneous confidence intervals and the Bonferroni intervals
# are built through different methods
# => they can lead to different conclusions (although in this particular 
#    case the conclusion is the same)

# Note 2: Here we have three confidence regions:
# (1) ellipse, (2) T2 rectangle, (3) Bonferroni rectangle
# In general, univariate confidence intervals for k given directions are
# associated with polygonal regions with 2*k sides.
# In this example, for k=2 the Bonferroni rectangle is smaller than the
# T2 rectangle. For increasing values of k, the Bonferroni region will 
# become larger and larger

### For instance, let's consider 3 linear combinations:
k <- 3

quartz()
par(mfrow=c(1,3))

# Scatter plot of DSS and DBOD
plot(D, asp=1, main='Confidence regions (k=3)')

# Adding the 95% confidence region for the true mean and delta.0
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=1, col='grey')
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

# Let's consider the following directions (DSS, DBOD and DSS + DBOD)
theta <- c(0, pi/4, pi/2)

# We draw the lines defining the Bonferroni confidence region (which is no
# longer a rectangle but a polygon)
for(i in 1:length(theta)) {
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='blue', lty=1,lwd=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='blue', lty=1,lwd=1)
}

# Same for T2 simultaneous
for(i in 1:length(theta)) {
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth),
         col='red', lwd=1, lty=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth), 
         col='red', lwd=1, lty=1)
}
# We add the legend
legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('blue', 'red'), lty=1)

### let's add another linear combination
k <- 4

# Same plot as before with the fourth linear combination added
plot(D, asp=1, main='Confidence regions (k=4)')
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=1, col='grey')
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

theta   <- c(theta,3*pi/4)

for(i in 1:length(theta))
{
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='orange', lty=1,lwd=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='orange', lty=1,lwd=1)
}

for(i in 1:length(theta))
{
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth),
         col='red', lwd=1, lty=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth), 
         col='red', lwd=1, lty=1)
}
legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('orange', 'red'), lty=1)

### let's add another linear combination
k <- 5

# Same plot as before with the fifth linear combination added
plot(D, asp=1, main='Confidence regions (k=5)')
ellipse(center=D.mean, shape=D.cov/n, radius=sqrt(cfr.fisher), lwd=1, col='grey')
points(delta.0[1], delta.0[2], pch=16, col='grey35', cex=1.25)

theta   <- c(theta,pi/6)

for(i in 1:length(theta))
{
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='forestgreen', lty=1,lwd=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * qt(1 - alpha/(2*k), n-1)) * a - 100*a.orth), 
         col='forestgreen', lty=1,lwd=1)
}

for(i in 1:length(theta))
{
  a       <- c( cos(theta[i]), sin(theta[i]))
  a.orth  <- c(-sin(theta[i]), cos(theta[i]))
  lines( rbind(D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean + as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth),
         col='red', lwd=1, lty=1)
  lines( rbind(D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a + 100*a.orth,
               D.mean - as.vector(sqrt( var(as.matrix(D) %*% a) / n ) * sqrt(cfr.fisher)) * a - 100*a.orth), 
         col='red', lwd=1, lty=1)
}
legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('forestgreen', 'red'), lty=1)

# The Bonferroni region is now bigger than the T2 region

graphics.off()

#_______________________________________________________________________________
##### Test for repeated measures 
#####----------------------------

### Exercise similar to Pb 3 of 14/09/2006
###------------------------------------------
# A pharmaceutical company performed a clinical test on 50 rats to investigate
# the effect of a new drug on the blood pressure.
# The blood pressure was measured to each rat four times: just before 
# giving the drug, 8, 16, 24 hours after the drug was given.
# (a) Perform a test at level 5% to prove that the drug has influence on 
#     the blood pressure during the 24 hours
# (b) Highlight the effect of the drug on the blood pressure

pressure <- read.table('pressure.txt',
                       col.names = c('h.0', 'h.8', 'h.16', 'h.24'))
head(pressure)
dim(pressure)

# Test for multivariate normality (Henze-Zirkler test by default)
mvn(pressure)$multivariateNormality

quartz()

# Plotting each observation as a line with 'h.0', 'h.8', 'h.16', 'h.24' on x-axis
matplot(t(pressure), type='l', lty = 1)

### Question (a) ###

n <- dim(pressure)[1]
q <- dim(pressure)[2]

M <- sapply(pressure, mean) # sample mean
M
S <- cov(pressure) # covariance matrix
S

# we build one of the possible contrast matrices to answer
# the question
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
C
# here we are looking at the effects on the pressure
# after 8, 16 and 24 hours from the instant the drug was
# given

# Test: H0: C %*% mu == c(0, 0, 0) vs H1: C %*% mu != c(0, 0, 0)
alpha   <- .05
delta.0 <- c(0, 0, 0)

Md <- C %*% M # Sample mean of the "constrated" observations
Sd <- C %*% S %*% t(C) # Sample covariance of the constrasted observations
Sdinv <- solve(Sd)

# Hotelling T2 statistics
T2 <- n * t(Md - delta.0) %*% Sdinv %*% (Md - delta.0)

# (q-1)*(n-1)/(n-(q-1)) times the 1-alpha Fisher quantile with q-1 and n-q+1 df 
cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha, (q-1), n-(q-1)) 

T2 < cfr.fisher # Testing if we are in the rejection region
T2
cfr.fisher

# T2 is much higher than cfr.fisher => the p-value will be very small
P <- 1 - pf(T2*(n-(q-1))/((q-1)*(n-1)), (q-1), n-(q-1))
P

### question (b)

# It is implicitly asking for confidence intervals on the components
# (for the mean of the increments after 8 hours, 16 hours and 24 hours)

# Simultaneous T2 intervals in the direction of the contrasts h.* - h.0
IC.T2 <- cbind(Md - sqrt(cfr.fisher*diag(Sd)/n),
               Md,
               Md + sqrt(cfr.fisher*diag(Sd)/n))
IC.T2

# Bonferroni intervals 
k <- q-1   # number of increments (i.e., dim(C)[1])
cfr.t <- qt(1-alpha/(2*k), n-1)

IC.BF <- cbind(Md - cfr.t*sqrt(diag(Sd)/n),
               Md,
               Md + cfr.t*sqrt(diag(Sd)/n))
IC.BF


quartz()
matplot(t(matrix(1:3, 3, 3)), t(IC.BF), type='b', pch='', xlim=c(0,4), xlab='',
        ylab='', main='Confidence intervals')

# Plotting the Bonferroni intervals
segments(matrix(1:3, 3, 1), IC.BF[,1], matrix(1:3, 3, 1), IC.BF[,3],
         col='orange', lwd=2)
points(1:3, IC.BF[,2], col='orange', pch=16)

# Plotting delta.0 under H0 (delta.0 == c(0, 0, 0))
points(1:3+.05, delta.0, col='black', pch=16)

# Plotting the simultaneous T2
segments(matrix(1:3+.1,3,1),IC.T2[,1],matrix(1:3+.1,3,1),IC.T2[,3], col='blue', lwd=2)
points(1:3+.1,IC.T2[,2], col='blue', pch=16)

legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('orange', 'blue'), lty=1, lwd=2)


### what happens if we change the constrast matrix?
Cbis <- matrix(c(-1, 1, 0, 0,
                 0, -1, 1, 0,
                 0, 0, -1, 1), 3, 4, byrow=T)
Cbis
# in this way we are looking at the mean increment of the pressure every 8 hours

Mdbis <- Cbis %*% M 
Sdbis <- Cbis %*% S %*% t(Cbis)
Sdinvbis <- solve(Sdbis)

T2bis <- n * t(Mdbis) %*% Sdinvbis %*% Mdbis

T2bis < cfr.fisher

# compare the T2 test statistics associated with C and Cbis
T2bis
T2


# What is changed?
# The confidence intervals on the contrasts
# (because we are looking at different contrasts!)

# Bonferroni intervals for Cbis
IC.BFbis <- cbind(Mdbis - cfr.t*sqrt(diag(Sdbis)/n),
                  Mdbis,
                  Mdbis + cfr.t*sqrt(diag(Sdbis)/n))

# Sim. T2 intervals for Cbis
IC.T2bis <- cbind(Mdbis - sqrt(cfr.fisher*diag(Sdbis)/n),
                  Mdbis,
                  Mdbis + sqrt(cfr.fisher*diag(Sdbis)/n))

IC.BFbis
IC.BF

IC.T2bis
IC.T2

### what if we want to verify the following hypothesis:
### "the drug decreases the pressure of two units with respect to
### the baseline at both 8 and 16 hours, and its effect vanishes in 24 hours
### from the drug administration"

C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)

delta.0 <- c(-2, -2, 0)

# or

C <- matrix(c(-1, 1, 0, 0,
              0, -1, 1, 0,
              0, 0, -1, 1), 3, 4, byrow=T)

delta.0 <- c(-2, 0, 2)

Md <- C %*% M # sample mean for the contrasted observations
Sd <- C %*% S %*% t(C) # sample covariance matrix
Sdinv <- solve(Sd)

# Hotelling T2 statistic
T2 <- n * t(Md - delta.0) %*% Sdinv %*% (Md - delta.0)

# (q-1)*(n-1)/(n-q+1) times the 1-alpha Fisher quantile with q-1 and n-q+1 df
cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1))

T2 < cfr.fisher # Do we accept H0?
T2
cfr.fisher

# p-value
P <- 1 - pf(T2*(n-(q-1))/((q-1)*(n-1)), (q-1), n-(q-1))
P

#_______________________________________________________________________________
##### Test for the mean of two independent Gaussian populations
#####------------------------------------------------------------

# we build the data: two bivariate samples of resp. 3 and 4 obs
t1 <- matrix(c(3, 3, 1, 6, 2, 3), 2)
t1 <- data.frame(t(t1))
t2 <- matrix(c(2, 3, 5, 1, 3, 1, 2, 3),2)
t2 <- data.frame(t(t2))

t1
t2

n1 <- dim(t1)[1] # n1 = 3
n2 <- dim(t2)[1] # n2 = 4
p  <- dim(t1)[2] # p = 2

# Test for multivariate normality -> useless regarding the sample sizes

# we compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)

# We compare the matrices -> here, using rule of thumb:
# we don't reject equality of covariance matrices if s1_ii and s2_ii differ from
# less than a factor ~10 (see J-W p.291)
list(S1=t1.cov, S2=t2.cov, Spooled=Sp)

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1 - mu2 == c(0, 0)  vs  H1: mu1 - mu2 != c(0, 0)

alpha   <- .01
delta.0 <- c(0, 0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # TRUE: no statistical evidence to reject H0 at level 1%

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  
# P-value high (we don't reject at 1%,5%,10%)

# Simultaneous T2 intervals
IC.T2.X1 <- c(t1.mean[1]-t2.mean[1] - sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)),
              t1.mean[1]-t2.mean[1] + sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)))
IC.T2.X2 <- c(t1.mean[2]-t2.mean[2] - sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)),
              t1.mean[2]-t2.mean[2] + sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)))

IC.T2 <- rbind(IC.T2.X1, IC.T2.X2)
dimnames(IC.T2)[[2]] <- c('inf','sup')                        
IC.T2

graphics.off()

#_______________________________________________________________________________
##### Pb 2 of 4/07/2013
#####--------------------
# The Chinese Institute of Genomics has reconstructed the genetic map of 200
# Chinese individuals, among which 100 were allergic to almonds and 100 were not
# allergic to almonds. The files hatingalmonds.txt e lovingalmonds.txt collect
# the presence (1) or absence (0) of 520 mutations, for the two groups respectively.
# The geneticist suspects that some of these mutations might be positively associated
# with the allergy to almonds, and for this reason they decide to compare the
# incidence of the 520 mutations in the two populations.

# a) Report the significant mutations, imposing a probability of at most 1%
#    that the single mutation is judged as influential if it isn't.
# b) Report the significant mutations, imposing a probability of at most 1%
#    that at least one of the non-influential mutations is judged as influential.
# c) Report the significant mutations, imposing that the expected proportion 
#    of mutations erroneously judged to be influential among all those judged to be
#    influential is lower than or equal to 1%

# You'll have to correct the univariate p-values of multiple tests.
# R function to correct p-values: p.adjust

allergy <- read.table('hatingalmonds.txt')
head(allergy)
dim(allergy)

noallergy <- read.table('lovingalmonds.txt')
head(noallergy)
dim(noallergy)

n1 <- dim(allergy)[1]
n2 <- dim(noallergy)[1]
p <- dim(noallergy)[2]

x.mean1 <- sapply(allergy, mean)
x.mean2 <- sapply(noallergy, mean)

p.hat <- (x.mean1*n1+x.mean2*n2)/(n1+n2)
x.var <- (p.hat*(1-p.hat))

# Test: H0.i: mu.i1 == mu.i2  vs H1.i: mu.i1 != mu.i2
# Asymptotic Z-test for the comparison of proportions (univariate and n=100 -> OK)

z.i <- (x.mean1-x.mean2)/sqrt(x.var*(1/n1+1/n2))
p.i <- ifelse(z.i<0, 2*pnorm(z.i), 2*(1-pnorm(z.i)))

which(p.i<.01)

# Bonferroni test
k <- 520

which(p.i*k<.01)

# or
p.Bf <- p.adjust(p.i, method='bonferroni')

which(p.Bf<.01)  

# Benjamini-Hockberg (control the false discovery rate)  
p.BH <- p.adjust(p.i, method='BH')

which(p.BH<.01)


quartz()
par(mfrow=c(1,3))
plot(p.i, main='Univariate')
abline(h=.01, lwd=2, col='red')

plot(p.Bf, main='Corrected - Bonferroni')
abline(h=.01, lwd=2, col='red')

plot(p.BH, main='Corrected - BH')
abline(h=.01, lwd=2, col='red')



### Additional execise on mean differences between independent Gaussian: 
### Pb 1 of 19/01/2022

# The Galician Food Association has launched an award for the Best Pulpo a la Gallega (meaning 
# Galician-style octopus). As part of the challenge, two tasters are sent to evaluate the 
# 30 finalist octopus dishes in A Coruna and the 30 finalist octopus dishes in Pontevedra. 
# Files "acoruna.txt" and "pontevedra.txt" collect the evaluations on each of the finalist 
# dishes given by the tasters in A Coruna and Pontevedra, respectively. Assume the evaluations 
# on different dishes to be independent, and the evaluations of the two tasters on the same 
# octopus dish to come from a bivariate Gaussian distribution.
# a) Perform a statistical test of level 99\% to verify if the mean evaluations in the two 
#    cities differ. State and verify the model assumptions. 
# b) Interpret the results of the test at point (a) through two Bonferroni intervals of
#    global level 99% for appropriate differences in the mean. Comment the result. 
# c) Is there statistical evidence to state that, at level 99%, the average evaluations of
#    A Coruna's octopus dishes are in mean higher than those of Pontevedra's octopus dishes? 
#
#    [by average evaluation of a octopus dishes is meant the one obtained by averaging the 
#    evaluations of the two tasters on that dishes]

