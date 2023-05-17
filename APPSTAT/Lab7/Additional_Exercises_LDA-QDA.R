setwd("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 8 - 19042021")
load("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 5 - 23042021/mcshapiro.test.RData")

library(MASS)

#_______________________________________________________________________________
##### Problem 2 of 1/07/2009
#####--------------------------
# An art historian requires your help to identify a criterion of classification 
# to discriminate the sculptures created by Gian Lorenzo Bernini from those of 
# other contemporary sculptors, based on the weight [tons] and height [m] 
# of 100 sculptures of undoubted attribution (sculptures.txt files). Taking into
# account that Bernini's sculptures are about 25% of the sculptures which have 
# to be classified and that the purpose of the historian is to minimize the expected
# number of misclassifications:
# a) build two classifiers C1 and C2, respectively, assuming for C1 that the data
#    come from two normal populations with equal covariance matrix, and for C2 that
#    the data come from two normal populations with different covariance matrix;
# b) estimate by cross-validation the AER of the two classifiers and comment their
#    values;
# c) how will be classified by the two classifiers a 2 meter high and 4 tons heavy
#    statue?

sculpt <- read.table('sculptures.txt', header=T)
head(sculpt)

autore <- factor(sculpt[,3], levels=c('Bernini', 'Altro'))
autore

bernini <- sculpt[1:50,1:2]
altro <- sculpt[50:100,1:2]

# question a)
mcshapiro.test(bernini)
mcshapiro.test(altro)

# LDA
lda.s <- lda(sculpt[,1:2], autore, prior=c(0.25, 0.75))
lda.s

# QDA
qda.s <- qda(sculpt[,1:2], autore, prior=c(0.25, 0.75))
qda.s

x11()
plot(sculpt[,1:2], main='Sculptures', xlab='Height', ylab='Weight', pch=20)
points(bernini, col='red', pch=20)
points(altro, col='blue', pch=20)
legend('bottomleft', legend=levels(autore), fill=c('red','blue'), cex=.7)

points(lda.s$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)

x  <- seq(min(sculpt[,1]), max(sculpt[,1]), length=200)
y  <- seq(min(sculpt[,2]), max(sculpt[,2]), length=200)
xy <- expand.grid(Altezza=x, Peso=y)

z  <- predict(lda.s, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

z.q  <- predict(qda.s, xy)$post  
z1.q <- z.q[,1] - z.q[,2] 
z2.q <- z.q[,2] - z.q[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

contour(x, y, matrix(z1.q, 200), levels=0, drawlabels=F, add=T, lty=2)  
contour(x, y, matrix(z2.q, 200), levels=0, drawlabels=F, add=T, lty=2)

dev.off()

# question b)
# LDA
LdaCV.s <- lda(sculpt[,1:2], autore, prior=c(0.25, 0.75), CV=T)
table(class.true=autore, class.assignedCV=LdaCV.s$class)

AER.CV.l <- 45/50*0.25+3/50*0.75
AER.CV.l

# QDA
QdaCV.s <- qda(sculpt[,1:2], autore, prior=c(0.25, 0.75), CV=T)
table(class.true=autore, class.assignedCV=QdaCV.s$class)

AER.CV.q <- 20/50*0.25+6/50*0.75
AER.CV.q

# question c)
predict(lda.s, c(Altezza=2, Peso=4))
predict(qda.s, c(Altezza=2, Peso=4))

points(2,4, pch=3, col='springgreen', lwd=2)

graphics.off()

Qda.m <- predict(qda.s)
table(class.true=autore, class.assigned=Qda.m$class)

#_______________________________________________________________________________
##### Problem 2 of 16/07/2010
#####--------------------------
# A young Portuguese engineer wants to build a machine able to automatically
# distinguish between two species of sardines (the Atlantic sardines and the
# Iberian sardine) on the basis of the length [cm] of the sardine. To
# do this, it aims to build a classifier based on the measurements of 500 
# Atlantic sardines and 500 Iberian sardine (sardine.txt file).
# a) Perform two tests to verify the hypothesis of normality of the two 
#    populations, a test for equality of means, a test for equality of 
#    variances.
# b) On the basis of the previous tests and knowing that 75% of fished sardines
#    belong to the Atlantic species while 25% to the Iberian species, build a
#    classifier that minimizes the number of misclassified sardines and report
#    the parameters.
# c) Estimate the AER of the classifier analytically using the estimated 
#    probability densities.

sardine <- read.table('sardine.txt', header=T)
head(sardine)

# question a) 
sardinaa <- sardine[[1]]
sardinai <- sardine[[2]]

shapiro.test(sardinaa)
shapiro.test(sardinai)
var.test(sardinaa, sardinai)
t.test(sardinaa, sardinai, var.eq=T)

# question b) and c)

MA <- mean(sardinaa)
MI <- mean(sardinai)
SD <- sqrt((var(sardinaa) + var(sardinai))/2)

# Analitically
misclass <- function(x){
  0.25*(1 - pnorm(x, MI, SD)) + 0.75*pnorm(x, MA, SD)
  }
optimize(f=misclass, lower=min(sardine), upper=max(sardine)) 

R   <- optimize(f=misclass, lower=min(sardine), upper=max(sardine))$minimum
AER <- optimize(f=misclass, lower=min(sardine), upper=max(sardine))$objective

# check
AER == 0.25*(1 - pnorm(R, mean(sardinai), SD)) + 0.75*pnorm(R, mean(sardinaa), SD)

# some plots
specie <- factor(rep(c('Atlantica', 'Iberica'), each=500))
lunghezza <- c(sardine[[1]], sardine[[2]])

fit <- lda(specie ~ lunghezza, prior=c(0.75, 0.25))
x <- data.frame(lunghezza=seq(min(sardine), max(sardine), 0.05))

LDA <- predict(fit, x)$posterior[,2] # classe iberica

x11()
par(mfrow=c(2,1))
plot(x[,1], 0.25*(dnorm(x[,1], MI, SD)), type='l', col='red', lty=1, xlab='x', ylab='density * prior', ylim=c(0,0.6))
lines(x[,1], 0.75*(dnorm(x[,1], MA, SD)), type='l', col='blue', lty=1, xlab='x', ylab='density * prior')
abline(v=R, lty=2)
points(sardinaa, rep(0, 500), pch=3, col='blue')
points(sardinai, rep(0, 500), pch=3, col='red')
legend(7,.5,legend=c('Atlantic','Iberian'),fill=c('blue','red'),cex=.7)

plot(x[,1], LDA, type='l', col='red', lty=1, xlab='x', ylab='estimated posterior')# rosso = iberica
lines(x[,1], 1 - LDA, type='l', col='blue', lty=1, xlab='x', ylab='estimated posterior')# blu = atlantica
abline(h = 0.5)
abline(v=R, lty=2)

points(sardinaa, rep(0, 500), pch=3, col='blue')
points(sardinai, rep(0, 500), pch=3, col='red')

dev.off()

# Compute the APER
prior <- c(0.75,0.25)
G <- 2
misc <- table(classe.vera=specie, classe.allocata=predict(fit)$class)
misc

APER <- 0
for(g in 1:G)
  APER <- APER + sum(misc[g,-g])/sum(sum(misc[g,])) * prior[g]  
APER
AER

#_______________________________________________________________________________
##### Problem 2 of 14/02/2011
#####--------------------------
# The ATM is considering the possibility of including in the turnstiles optical
# readers able to measure the size [mm] of tickets.
# In order to build a suitable software, it asks you to build a classification 
# rule which minimizes the expected number of errors.
# Starting with the measures relating to 100 regular tickets (true.txt file)
# and 100 counterfeit tickets (false.txt files) and knowing that the
# 0.5% of the banknotes in circulation is counterfeit:
# a) construct an appropriate classification rule (in particular, verify the 
#    assumptions when possibile, and provide a qualitative graph of the 
#    classification regions);
# b) compute the APER and discuss its value;
# c) on the basis of the rule identified at point (a), how will be classified
#    a ticket long 85.5 mm and wide 55.0 mm? What is the probability that the 
#    ticket will be false?

true <- read.table('true.txt', header=T)
false <- read.table('false.txt', header=T)
head(true)
head(false)

fv <- as.factor(rep(c('vero','falso'), c(100,100)))
fv

biglietti <- rbind(true, false) 

### question a)
# normality
mcshapiro.test(true)
mcshapiro.test(false)

# homogeneity of covariances
S1<-cov(true)
S2<-cov(false)
x11()
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))

x11()
plot(biglietti, main='Tickets', xlab='Length', ylab='Width', pch=20)
points(false, col='red', pch=20)
points(true, col='blue', pch=20)
legend('bottomleft', legend=levels(fv), fill=c('red','blue'), cex=.7)

# The assumption of homogeneity of covariances doesn't seem to be satisfied
dev.off()
dev.off()

# QDA
qda.bigl <- qda(biglietti, fv, prior = c(0.005, 0.995))
qda.bigl # look at the order in levels to set the priors!!
Qda.bigl <- predict(qda.bigl, biglietti)
Qda.bigl

x11()
plot(biglietti, main='Tickets', xlab='Length', ylab='Width', pch=20)
points(false, col='red', pch=20)
points(true, col='blue', pch=20)
legend('bottomleft', legend=levels(fv), fill=c('red','blue'), cex=.7)
points(qda.bigl$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)

x  <- seq(min(biglietti[,1]), max(biglietti[,1]), length=200)
y  <- seq(min(biglietti[,2]), max(biglietti[,2]), length=200)
xy <- expand.grid(Lunghezza=x, Larghezza=y)

z  <- predict(qda.bigl, xy)$post   
z1 <- z[,1] - z[,2]   
z2 <- z[,2] - z[,1]      

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

dev.off()

### question b)

MC <- table(classe.vera=fv, classe.allocata=Qda.bigl$class)
MC
APER <- 66/100 * 0.005 
APER

### domanda c)

predict(qda.bigl, cbind(Lunghezza = 85.5, Larghezza = 55))

#_______________________________________________________________________________
##### Problem 3 of 12/02/2014
#####--------------------------

# To cope with the economic crisis, the Hotel Mont Blanc in Courmayeur has 
# decided to apply special discounts on their Carnival rates in case the 
# snow at skiing facilities is predicted to be of poor quality. The file neve.txt 
# reports, for the last 60 years, the data on total snowfall [cm] and the medium
# temperature [°C] recorded in the two months from December to January, together 
# with the judgment on the quality of the Carnival snow provided by the Alpine Guides 
# Society of Courmayeur.
# a) Build a classifier for the quality of the Carnival snow that minimizes the expected
#    cost of misclassification (display a qualitative graph of the classification regions)
#    when assuming that:
#    - There is no economic loss in the case in which the snow is good and no discount
#      is applied; there is no economic loss in case the snow is bad and the Hotel
#      operates the discount; there is a loss of 3000 euros in the case the snow is bad
#      and no discounts are applied; there is a loss of 2000 euros in the case the snow is 
#      good and discounts are applied;
#    - A season characterized by good snow is associated with a higher variability in
#      temperatures and in the amount of snow.
# b) Compute the APER of the classifier.
# c) Based on the estimates at point (b), estimate the expected economic loss of the 
#    classifier.
# d) Based on the classifier build at point (a) and knowing that the last two months
#    December-January a total of 200 cm of snow have fallen and the average temperature
#    was -4 °C, would you recommend to the hotel to apply the special discount of Carnival?

# question a)
neve <- read.table('neve.txt', header=T)
good<-neve[neve[,3]=='good',1:2]
bad<-neve[neve[,3]=='bad',1:2]

mcshapiro.test(good)$pvalue
mcshapiro.test(bad)$pvalue

prior <- c(dim(bad)[1]/(sum(dim(bad)[1],dim(good)[1])),dim(good)[1]/(sum(dim(bad)[1],dim(good)[1])))
pb <- prior[1]
pg <- prior[2]

c.bg <- 2000
c.gb <- 3000

# Modified prior to account for the misclassification costs
prior.c <- c(bad=pb*c.gb/(c.bg*pg+c.gb*pb),good=pg*c.bg/(c.bg*pg+c.gb*pb))
prior.c

qda.m <- qda(giudizio ~ quantita + temperatura, data=neve, prior=prior.c)
qda.m

x11()
plot(neve[,1:2], main='Snow', xlab='V1', ylab='V2', pch=20)
points(bad, col='red', pch=20)
points(good, col='blue', pch=20)
legend('bottomleft', legend=levels(neve[,3]), fill=c('red','blue'), cex=.7)

points(qda.m$means, pch=4,col=c('red','blue') , lwd=2, cex=1.5)

x  <- seq(min(neve[,1]), max(neve[,1]), length=200)
y  <- seq(min(neve[,2]), max(neve[,2]), length=200)
xy <- expand.grid(quantita=x, temperatura=y)

z  <- predict(qda.m, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

# question b) [APER]
Qda.m <- predict(qda.m)
table(classe.vera=neve[,3], classe.allocata=Qda.m$class)

APER  <- (2+3)/(46+14)
APER

# question c) [Expected economic loss]
(3*c.bg+2*c.gb)/60

# question d) [Classification of 2012-2013]
z0.new <- data.frame(quantita=200, temperatura=-4)
points(z0.new, pch=4, col='springgreen', lwd=2, cex=1.5)

predict(qda.m,z0.new)$class

graphics.off()

#_______________________________________________________________________________
##### Problem 2 of 9/09/2009
#####--------------------------

# A young statistician from Tromso wants to build a classifier able to
# predict the presence or absence of snowfall on the day D + 1 using 
# the average temperature [°C] and humidity of day G. Using 
# the data of the last 30 days (file snow.txt) and knowing that in the
# last 40 years in Tromso it has snowed an average of 207 days a year:
# a) build a classifier [report the model assumptions and provide a 
#    qualitative graph of the classification regions].
# b) Estimate the APER of classifier (a), and compare it with that of the
#    trivial classifier.
# c) The temperature and the humidity measured yesterday in Tromso are 
#    -10 Â° C and 0.75, respectively. Estimate the likelihood of snow 
#    for today by using both the classifier (a) that the trivial classifier.


snow <- read.table('snow.txt', header=TRUE)
snow

attach(snow)

i1 <- which(Snow.G.1=='no-snow')
i2 <- which(Snow.G.1=='snow')

# question a)
mcshapiro.test(snow[i1,1:2])$p
mcshapiro.test(snow[i2,1:2])$p

x11()
plot(Temperature.G,Humidity.G, pch=19, col=ifelse(Snow.G.1=='no-snow','blue','lightblue'))

dev.off()

# QDA
library(MASS)
qda.s <- qda(snow[,1:2], snow[,3], prior=c(1-207/365,207/365))
qda.s

x11()
plot(snow[,1:2], main='Snow', pch=20, col=ifelse(Snow.G.1=='no-snow','blue','lightblue'))
legend('bottomleft', legend=levels(snow[,3]), fill=c('blue','steelblue2'), cex=.7)

points(qda.s$means, pch=4,col=c('blue','steelblue2') , lwd=2, cex=1.5)

x  <- seq(min(snow[,1]), max(snow[,1]), length=200)
y  <- seq(min(snow[,2]), max(snow[,2]), length=200)
xy <- expand.grid(Temperature.G=x, Humidity.G=y)

z  <- predict(qda.s, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]  

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)

dev.off()

# question b)
Qda.s <- predict(qda.s)
table(classe.vera=snow[,3], classe.allocata=Qda.s$class)

prior <- c(1-207/365,207/365)

APER  <- 2/10*prior[1]+1/20*prior[2]
APER

# Trivial classifier; classifies always as the most likely class a priori
APER.banale <- prior[1]
APER.banale

# question c)
new.day <- c(Temperature.G=-10, Humidity.G=0.75)
predict(qda.s, new.day)$posterior[2]

prior.snow=207/365
prior.snow

graphics.off()
