# The file musicCountry.txt contains the price (in e) and the average song length (in minutes) of 188 albums
# released in Germany and US.
# a) Knowing that, on average, 90% of albums are released in US, build a classifier to characterize the country of
# release based on the price and the average song length. Report the model for the data, the estimates of its
# parameters (means and covariances) and verify the model assumptions. Report the plot of the classification
# regions.
# b) Estimate the AER of the classifier through leave-one-out cross-validation.
# c) Using the classifier built at point (a), what is the estimated probability that a new album is classified as US?
#   d) How would you classify a new album with with a price of 50e, and an average song length of 3.5 minutes?
#   e) Use a support vector machine with linear kernel to classify albums. Tune the cost parameter with a 10-fold
# cross validation choosing between the values: 0.001, 0.01, 0.1, 1, 10, 100. Report the chosen cost and a plot of
# the classification regions. How would you classify the album at point d) with this classifier?

library(MASS) #lda and qda
library(e1071) #svm

#Normality test
load("~/Desktop/Applied Statistics/Labs/Lab 5/mcshapiro.test.RData")

musicCountry <- read.table('musicCountry.txt', header=T)


# a) Knowing that, on average, 90% of albums are released in US, build a classifier to characterize the country of
# release based on the price and the average song length. Report the model for the data, the estimates of its
# parameters (means and covariances) and verify the model assumptions. Report the plot of the classification
# regions.

germany = musicCountry[musicCountry$release.country=="Germany",1:2]
usa = musicCountry[musicCountry$release.country=="US",1:2]

mcshapiro.test(germany)$p
mcshapiro.test(usa)$p

S1 <- cov (germany)
S2 <- cov (usa)

par(mfrow=c(2,1))
image(S1)
image(S2)

fact <- as.factor(musicCountry[,3])
prior.c <- c(0.1,0.9)


lda.m <- lda(musicCountry[,1:2], fact, prior=prior.c)
#mean
lda.m$means
# price average.length
# Germany 59.82982       5.891951
# US      29.74888       4.377556

n1 <- dim(germany)[1]
n2 <- dim(usa)[1]
n <- n1+n2
Sp <- 1/(n-2) * ((n1-1)*S1 + (n2-1)*S2)
Sp
# price average.length
# price          161.37898      -2.671480
# average.length  -2.67148       2.853387

plot(musicCountry[,1:2], main='Music', xlab='price', ylab='avarage.length', pch=20)
points(germany, col='red', pch=20)
points(usa, col='blue', pch=20)
legend('topleft', legend=c('Germany','Usa'), fill=c('red','blue'), cex=.7)
x  <- seq(min(musicCountry[,1]), max(musicCountry[,1]), length=200)
y  <- seq(min(musicCountry[,2]), max(musicCountry[,2]), length=200)
xy <- expand.grid(price=x, average.length=y)

z  <- predict(lda.m, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# b) Estimate the AER of the classifier through leave-one-out cross-validation.

lda.CV <- lda(musicCountry[,1:2], fact, prior=prior.c,CV=T) 
misc <- table(class.true=fact, class.assignedCV=lda.CV$class)
G <- 2 #number of group
AER <- 0
for(g in 1:G)
  AER <- AER + sum(misc[g,-g])/sum(misc[g,]) * prior.c[g] 
AER
#0.04795322

# c) Using the classifier built at point (a), what is the estimated probability that a new album is classified as US?

lda.z <- predict(lda.m, xy)
colMeans(lda.z$posterior) #classification probability

#   d) How would you classify a new album with with a price of 50e, and an average song length of 3.5 minutes?


xy2 <- expand.grid(price=50, average.length=3.5)
lda.z <- predict(lda.m, xy2)
lda.z$class #classification class
#US



#   e) Use a support vector machine with linear kernel to classify albums. Tune the cost parameter with a 10-fold
# cross validation choosing between the values: 0.001, 0.01, 0.1, 1, 10, 100. Report the chosen cost and a plot of
# the classification regions. How would you classify the album at point d) with this classifier?

dat <- data.frame(x=musicCountry[,1:2], y=fact)

cost.c <- c(0.001,0.01,0.1,1,10,100)
best.c <- 0
best.err <- 1

for (i in cost.c) {
  svmfit <- svm( y ~ . , data=dat , kernel ='linear', cost =i , scale =F, CV=T)
  ygrid <- predict(svmfit,dat)
  misc <- table(class.true=fact, class.assignedCV=ygrid)
  err <- (misc[2]+misc[3])/188
  if (err<best.err){
    best.err <- err
    best.c <- i
  }
}
best.c
svmfit <- svm( y ~ . , data=dat , kernel ='linear', cost =best.c , scale =F, CV=T)

par(mfrow=c(1,2))
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19, asp=1)

#prediction
xy <- expand.grid(x.price=50, x.average.length=3.5)
ygrid <- predict(svmfit,xy)

#US
