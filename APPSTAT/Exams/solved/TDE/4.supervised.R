
#Libraries: 
library(MASS) #lda and qda
library(class) #knn
library(e1071) #svm

#Normality test
load("C:/Users/silvi/Desktop/4.APPLIED STATISTICS/labs/lab5/mcshapiro.test.RData")

### LDA 
###------------------
# Assumptions:
# 1) normality in each group
# 2) equal misclassification costs
# 3) same covariance in all groups (strong assumption)

### QDA 
###------------------
# Assumptions:
# 1) normality in each group
# 2) equal misclassification costs

#Example
data <- read.table(' .txt', header=T)

group = data$... #name of the class with the labels
tf <- factor(group=="NO") #group no-> true, the other is false
data_true = data[group=="NO",1:2] #name of "true" label
data_false = data[group=="AB",1:2] #name of "false" label

#see the data...
x11()
plot(data[,1:2],pch = 19)
points(data_true,col = 2, pch = 19)
points(data_false,col = 3, pch = 19)

#gaussianity within groups
mcshapiro.test(data_true)$p
mcshapiro.test(data_false)$p
#first look at the variance !!
#do groups have same covariance?
S1 <- cov(data_true)
S2 <- cov(data_false)
S1
S2
x11()
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
#same image -> I use LDA
#different image -> QDA
### per fare dei test sulla varianza quando ho solo una variabile:
var.test(sardinaa, sardinai)


#prior probability
prior_prob = c(pf = ..., pt = ...)

#if you have misclass costs: adjust the prior with the misclassification cost
c.true_false <- ....  #we predict true, but it is false
c.false_true <- ....
prior_prob <- c(pf*c.true_false, pt*c.false_true) /(pf*true_false + pt*c.false_true)

#model
lda.m <- lda(data[,1:2], tf, prior=prior_prob)
lda.m
#or
qda.m <- qda(data[,1:2], tf, prior=prior_prob)
qda.m
#rmk: ALWAYS CHECK THAT THE PRIOR PROBABILITY IS NOT EXCHANGED

#estimate of the paramenters
#mean
lda.m$means 
#covariance
S1 <- cov(data_true)
S2 <- cov(data_false)
#moreover, with lda
n1 <- dim(data_true)[1]
n2 <- dim(data_false)[1]
n <- n1+n2
Sp <- 1/(n-2) * ((n1-1)*S1 + (n2-1)*S2)
Sp #estimate of the covariance
#if more than 2 groups check the note

#plot
plot(data[,1:2], xlab='...', ylab='...', pch=20)
points(data_true, col='red', pch=20)
points(data_false, col='blue', pch=20)
legend('bottomleft', legend=levels(tf), fill=c('red','blue'), cex=.7)
x  <- seq(min(data[,1]), max(data[,1]), length=200)
y  <- seq(min(data[,2]), max(data[,2]), length=200)
xy <- expand.grid(...=x, ...=y) #metti nome delle variabili

z  <- predict(lda.m, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

#AER (by cross validation) (actual error rate)
lda.CV <- lda(data[,1:2], tf, prior=prior_prob,CV=T) #'same' for qda
misc <- table(class.true=tf, class.assignedCV = lda.CV$class)
G <- 2 #number of group
AER <- 0
for(g in 1:G)
  AER <- AER + sum(misc[g,-g])/sum(misc[g,]) * prior_prob[g] 
#AER per la knn -> cambio la misc e uso class.assignedCV = knn.m

#APER (approssimate error rate)
misc <- table(class.true=tf, class.assigned=predict(lda.m)$class)
G <- 2 #number of group
APER <- 0
for(g in 1:G)
  APER <- APER + sum(misc[g,-g])/sum(misc[g,]) * prior_prob[g] 

#Example
class.assignedCV
class.true true false
true     5    45		=> AER  = (45*pt+3*pf)/50
false    3    47   



#prediction
x <- 2
y <- 4
xy <- expand.grid(Altezza=x, Peso=y)
lda.z <- predict(lda.m, xy)
lda.z$class #classification class
lda.z$posterior #classification probability



####SE HO 3 GRUPPI: ####

data <- read.table('activity.txt', header=T)

group = factor(data$activity)
data_wal = data[1:150,1:2]
data_sit = data[151:300,1:2]
data_lay = data[301:450,1:2]

#see the data...
x11()
plot(data[,1:2],pch = 19)
points(data_wal,col = 2, pch = 19)
points(data_sit,col = 3, pch = 19)
points(data_lay,col = 4, pch = 19)

mcshapiro.test(data_wal)$p
mcshapiro.test(data_sit)$p
mcshapiro.test(data_lay)$p

#cov
S1 <- cov(data_wal)
S2 <- cov(data_sit)
S3 <- cov(data_lay)
S1
S2
S3
x11()
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))

image(S1)
image(S2)
image(S3)


#QDA
prior_prob = c(l=9/24, s=12/24, w=3/24) #da mettere in ordine alfabetico di factor

qda.m <- qda(data[,1:2], group, prior=prior_prob)
qda.m


#plot
x11()
plot(data[,1:2], xlab='...', ylab='...', pch=20)
points(data_wal, col='red', pch=20)
points(data_sit, col='green', pch=20)
points(data_lay, col='blue', pch=20)
legend('bottomright', legend=levels(group), fill=c('blue','green','red'), cex=.7)
#attenzione all'ordine dei nomi in group!!!!!! li mette in ordine alfabetico
x  <- seq(min(data[,1]), max(data[,1]), length=200)
y  <- seq(min(data[,2]), max(data[,2]), length=200)
xy <- expand.grid(accel=x, gyro=y) #metti nome delle variabili

z  <- predict(lda.iris, xy)$post  # these are P_i*f_i(x,y)  
z1 <- z[,1] - pmax(z[,2], z[,3])  # P_1*f_1(x,y)-max{P_j*f_j(x,y)}  
z2 <- z[,2] - pmax(z[,1], z[,3])  # P_2*f_2(x,y)-max{P_j*f_j(x,y)}    
z3 <- z[,3] - pmax(z[,1], z[,2])  # P_3*f_3(x,y)-max{P_j*f_j(x,y)}

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)

#### ####


#########################KNN classification #####################
#general command
knn.m <- knn(train = data[,1:2], test = data[,1:2], cl = tf, k = ..., prob=T)
#L1o cross validation (based on misclassification error)
n <- dim(data)[1]
best.k <- 0
best.err <- 100
best.model <- 0
#se chiede seed metti qui
for (k in 10:30) {
  tmp <- knn.cv(train = data[,1:2] , cl = tf,  k = k)  
  misc <- table(tf,tmp)
  if ((misc[2]+misc[3])/n < best.err){
    best.err <- (misc[2]+misc[3])/n
    best.k <- k
    best.model <- tmp
  }
}

#plot
plot(data[,1:2], main='', xlab='x', ylab='y', pch=20)
points(data_true, col='red', pch=20)
points(data_false, col='blue', pch=20)
legend('bottomleft', legend=levels(tf), fill=c('red','blue'), cex=.7)
x  <- seq(min(data[,1]), max(data[,1]), length=200)
y  <- seq(min(data[,2]), max(data[,2]), length=200)
xy <- expand.grid(...=x, ...=y)

knn.mt <- knn(train = data[,1:2], test = xy, cl = tf, k = best.k)
z <- as.numeric(knn.mt)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

#prediction
xy <- expand.grid(x=10.8, y=39.4)
knn.m <- knn(train = data[,1:2], test = xy, cl = tf, k = best.k, prob=T)
knn.m


########################## SVM ######################
dat <- data.frame(x=data[,1:2], y=tf)

svmfit <- svm( y ~ . , data=dat , kernel='linear', cost=0.1 , scale=F)
summary(svmfit) #osserva il numero di support vector
#the lower the cost, the higher the number of support vector!!!
#(more SV -> higher bias, lower variance)

#plot
par(mfrow=c(1,2))
plot(svmfit , dat, col =c('salmon', 'light blue'), pch=19)
#nel plot, le x rappresentano i support vector
#NB ascisse e ordinate sono invertite

#prediction
xy <- expand.grid(x.incidence=60, x.tilt=0)
ygrid <- predict(svmfit,xy)

#cross-validation for the cost:
set.seed (1)
tune.out <- tune(svm,y~.,data=dat ,kernel = 'linear',
                 ranges =list(cost=c(0.001 , 0.01, 0.1,10,100) ))
summary(tune.out)
# Extract the best model from the result of tune
bestmod <- tune.out$best.model
summary(bestmod)
plot(bestmod , dat, col =c('salmon', 'light blue'), pch=19, asp=1)




