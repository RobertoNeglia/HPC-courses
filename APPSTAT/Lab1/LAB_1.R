###--------------------###
### LAB 1 (02/03/2022) ###
###--------------------###

### TOPICS:
### Basic commands (scalars, vectors, matrices, and operations)
### Import/Export of dataframes
### Examples of univariate statistical analyses with plots
### Visualization of Multivariate Data
### Visualization of Categorical Data
### 3d plots, functions, "for" cycles
### Save plots


#_______________________________________________________________________________
##### Getting started

# To download R:   http://www.r-project.org/

# For further info on R:
# - material on the R website
# - "An Introduction to R" (Venable and Smith)
# - "help" button in the R interface


# for a "user friendly" interface --> RStudio
# (you can download it from rstudio.org)

# Comments: to comment a line insert "#" 

# To execute a command directly from the script (i.e., without copy & paste
# in the Console)
# With R: ctrl+r
# With RStudio: ctrl+enter 

# the command:
setwd("/home/rubuntu/HPC/APPSTAT/Lab1")

# selects the working directory, i.e., the directory in which R will 
# look for or save the files (data, plots, ..)
# NB: Not necessary if the working directory was already opened as a "Project"

# the command:

getwd()

# returns the current working directory

# R coding style guidelines
# https://rstudio-pubs-static.s3.amazonaws.com/390511_286f47c578694d3dbd35b6a71f3af4d6.html

#_______________________________________________________________________________
##### VECTORS/MATRICES & LINEAR ALGEBRA IN R

# scalars
a <- 1
a

# In the classical R syntax, values are assigned by using "<-". 
# In recent versions of R the assignment can be equivalently done 
# through "="
a = 1
a

4 -> b

a <- b
b = 5

# vectors

v <- c(2, 3, 5, 4)
v

u <- seq(2, 5, len=4)
u

u <- seq(2, 5, by=1)
u

z <- 2:5
z


# matrices
W <- rbind(c(11, 13, 15), c(12, 14, 16))
W

W <- cbind(c(11, 12), c(13, 14), c(15, 16))
W

W <- matrix(data = c(11, 12, 13, 14, 15, 16), nrow = 2, ncol = 3, byrow = F)
W

W <- matrix(c(11, 12, 13, 14, 15, 16), 2, 3)
W

################################################################################
#_______________________________________________________________________________
##### help

# if you remember the name of the command use
# help(COMMANDNAME)
# e.g.
help(matrix)

help.start()

################################################################################

#_______________________________________________________________________________
# Extraction of elements from a vector or a matrix

v
v[2]
v[c(2, 3)]
v[-3]
v[-c(1, 4)]

W
W[2, 3]
W[2, c(2, 3)]
W[2, ]
W[, c(2, 3)]

# Remark: in R vectors "are not" matrices n*1 o 1*n:
# vectors have only one index, whereas matrices have
# two indices (for rows and columns)

v
rbind(v)
cbind(v)

# Algebraic operations in R
# Remark: by default, operations are done component-wise

a <- 1
b <- 2

c <- c(2, 3, 4)
d <- c(10, 10, 10)

Z <- matrix(c(1, 10, 3, 10, 5, 10), nrow = 2, ncol = 3, byrow = F)

# Sum and multiplication (component-wise).
# (this default is different from that of matlab!)

a + b # scalar + scalar

c + d # vector + vector

a * b # scalar * scalar

c * d # vector * vector (component-wise)

c + a # vector + scalar

c^2 # attention: operations are always component-wise!

exp(c) 

sum(c) # sums the components of c

prod(c) # returns the product of the components of c

# Operations on matrices

V <- t(W) # transpose of a matrix

Z + W # matrix + matrix (component-wise)

Z * W # matrix * matrix (component-wise) 

V * W # matrix * matrix (component-wise) (error!)

V %*% W # Matrix multiplication
W %*% V 

W + a # matrix + scalar

W + c # matrix + vector

W + 2:5
# Remark: R uses the "recycling", i.e., it tries to make the
# terms dimensions compatible by recycling data if missing
W + 2:4 # recycling without warning!

# Inverse of a matrix (square and invertible)
A <- matrix(c(11, 13, 12, 14), ncol=2, nrow=2, byrow=TRUE)
det(A)
solve(A)

# Solution of a linear system Ax=b
b <- c(1, 1)
solve(A, b)

#_______________________________________________________________________________
##### Categorical data
# The command 'factor' converts the argument (vector of numbers or strings)
# in a vector of realizations of a categorical random variable, whose possible
# values are collected in 'Levels'

district <- c('MI', 'MI', 'VA', 'BG', 'LO', 'LO', 'CR', 'Alt', 'CR', 'MI',  
              'Alt', 'CR', 'LO', 'VA', 'MI', 'Alt', 'LO', 'MI')
district <- factor(district, levels=c('MI', 'LO', 'BG', 'CR', 'VA', 'Alt'))
district

resass <- table(district) # table of absolute frequencies
resass
resrel <- table(district) / length(district) # table of relative frequencies
resrel

#_______________________________________________________________________________
##### lists: objects made of objects (objects can be of different type)

exam <- list (course = 'Applied Statistics',  
              date = '27/09/2022',
              enrolled = 7,
              corrected = 6,
              student_id = as.character(c(45020, 45679, 46789, 43126, 42345,
                                          47568, 45674)),
              evaluation = c(30, 29, 30, NA, 25, 26, 27) 
)
exam

exam$evaluation
exam[[6]]

#_______________________________________________________________________________
##### data.frame: objects made of vectors of the same lengths,
##### possibly of different types.
##### (Remark: they look like matrices by they aren't!)

exam <- data.frame(
  student_id = factor(as.character(c(45020, 45679, 46789, 43126, 42345, 47568,
                                     45674))),
  evaluation_W = c(30, 29, 30, NA, 25, 26, 17), 
  evaluation_O = c(28, 30, 30, NA, 28, 27, NA), 
  evaluation_P = c(30, 30, 30, 30, 28, 28, 28),
  outcome  = factor(c('Passed', 'Passed', 'Passed', 'To be repeated', 'Passed',
                      'Passed', 'To be repeated')))
exam

exam$evaluation_W    # a data.frame is a particular kind of list!
exam[[2]]
exam[2, ]

evaluation_W

attach(exam)
evaluation_W
detach(exam)
evaluation_W

#_______________________________________________________________________________
##### Data Import/Export

record <- read.table('record.txt', header=T)
record

head(record)                                                           
dim(record)
dimnames(record)

# Transform times in seconds
record[, 4:7] <- record[, 4:7] * 60
record

# to save a data frame (or a matrix)
write.table(record, file = 'record_mod.txt')

# Remark. The file containing 'record_mod.txt' will be saved in the working directory 

# to save several objects in the workspace
W <- matrix(data = c(11, 12, 13, 14, 15, 16), nrow = 2, ncol = 3, byrow = F)
V <- t(W)
a <- 1

save(W, V, a, file = 'variousobjects.RData')

# to save the entire workspace: save.image('FILENAME.RData')
save.image("myworkspace.RData")

# this command remove all the variable of the workspace
ls()
rm(list=ls())

# to load a workspace (i.e., .RData)
load("variousobjects.RData")

ls()

#_______________________________________________________________________________
##### Example: analysis of quantitative data (with plots) 

record <- read.table('record_mod.txt', header=T)
record

# some synthetic indices
# mean of the columns
colMeans(record)
sapply(record, mean)
# standard deviation
sapply(record, sd)
# variance
sapply(record, var)
# covariance matrix
cov(record)
# correlation matrix
cor(record)

##### Descriptive/inferential analysis on the variable m100 ('very basic'!)
attach(record)

### t-test for the mean value of the quantity 
### H0: mean == 11.5 vs H1: mean != 11.5

# Recall: qqplot to verify (qualitatively) the Gaussian assumption on the
# distribution generating sample
qqnorm(m100) # quantile-quantile plot
qqline(m100, col='red') # theoretical line
# Recall: Shapiro-Wilk test to verify (quantitatively) the Gaussian assumption on the
# distribution generating sample
shapiro.test(m100)
# p-value = 0.2569
# this is greater than any significative value of alpha, which means that we can't 
# reject the assumption that the distribution is indeed Gaussian

# now move on testing about the mean of the distribution
alpha <- .05
mean.H0 <- 11.5

# automatically
t.test(m100, mu = mean.H0, alternative = 'two.sided', conf.level = 1-alpha)
# p-value > alpha --> cannot reject the null hypothesis
# accept H0 at 5 percent of significance
# reject H0 at 10% of significance

# manually
sample.mean <- mean(m100)
sample.sd <- sd(m100)
n <- length(m100)
tstat <- (sample.mean - mean.H0) / (sample.sd / sqrt(n))
cfr.t <- qt(1 - alpha/2, n-1)
abs(tstat) < cfr.t  # cannot reject H0 (accept H0)

pval  <- ifelse(tstat >= 0, (1 - pt(tstat, n-1))*2, pt(tstat, n-1)*2)
pval

IC <- c(inf     = sample.mean - sample.sd / sqrt(n) * qt(1 - alpha/2, n-1), 
        center  = sample.mean, 
        sup     = sample.mean + sample.sd / sqrt(n) * qt(1 - alpha/2, n-1))
IC

##### Simple linear regression (variable m200 vs m100)

# More than one plot in a unique device (commands par or layout)
# (command par)
quartz() # quartz() on PC
par(mfrow=c(2, 2))
hist(m100, prob=T, main="Histogram records 100m", xlab="sec")
hist(m200, prob=T, main="Histogram records 200m", xlab="sec")
boxplot(record[,1:2], main="Boxplot records 100m e 200m", xlab="sec")
plot(m100, m200, main='Scatter plot records 100m e 200m', xlab="Records 100m", ylab="Records 200m")

dev.off()

# command layout
quartz()
layout(cbind(c(1, 1), c(2, 3)), widths=c(2, 1), heights=c(1, 1))
plot(m100, m200)
hist(m100, prob=T)
hist(m200, prob=T)

dev.off()

# Fit of the linear model (command lm)
# Model: m200 = beta0 + beta1 * m100 + eps, eps ~ N(0, sigma^2)
regression <- lm(m200 ~ m100)
regression

summary(regression)

coef(regression)
vcov(regression)
residuals(regression)
fitted(regression)

quartz()
plot(m100, m200, asp=1, cex=0.75)
abline(coef(regression))
points(m100, fitted(regression), col='red', pch=19)

legend(
  'bottomright',
  c('Obs.', 'Fit', 'Reg. line'),
  col = c('black', 'red', 'black'),
  lwd = c(1, 1, 1),
  lty = c(-1, -1, 1),
  pch = c(c(1, 19, -1))
)

title(main='Linear regression (m200 vs m100)')

dev.off()

# Test F "by hand" (H0: beta1=0 vs H1: beta1!=0)
SSreg <- sum((fitted(regression) - mean(m200))^2)
SSres <- sum(residuals(regression)^2)
SStot <- sum((m200 - mean(m200))^2)

n <- length(m200)
Fstat <- (SSreg/1) / (SSres/(n-2))
P <- 1 - pf(Fstat, 1, n-2)
P # reject H0

# Confidence and prediction intervals (command predict)
newdata <- data.frame(m100=c(10, 11, 12))
pred_nd <- predict(regression, newdata)
pred_nd

IC_nd <- predict(regression, newdata, interval='confidence', level=.99)
IC_nd
IP_nd <- predict(regression, newdata, interval='prediction', level=.99)
IP_nd

quartz()
plot(m100, m200, asp=1, ylim=c(18.5, 27.5), cex=0.5)
abline(coef(regression))
points(m100, fitted(regression), col='red', pch=20)
points(c(10, 11, 12), pred_nd, col='blue', pch=16)

matlines(rbind(c(10, 11, 12), c(10, 11, 12)), t(IP_nd[, -1]), type="l", lty=2,
         col='dark grey', lwd=2)
matpoints(rbind(c(10, 11, 12), c(10, 11, 12)), t(IP_nd[, -1]), pch="-", lty=2,
          col='dark grey', lwd=2, cex=1.5)
matlines(rbind(c(10, 11, 12), c(10, 11, 12)), t(IC_nd[, -1]), type="l", lty=1,
         col='black', lwd=2)
matpoints(rbind(c(10, 11, 12), c(10, 11, 12)), t(IC_nd[, -1]), pch="-", lty=1,
          col='black', lwd=2, cex=1.5)

legend(
  'bottomright',
  c('Obs.', 'Fit', 'Reg. line', 'Pred. new', 'IC', 'IP'),
  col = c('black', 'red', 'black', 'blue', 'black', 'dark grey'),
  lwd = c(1, 1, 1, 1, 2, 2),
  lty = c(-1, -1, 1, -1, 1, 2),
  pch = c(c(1, 19, -1, 19, -1, -1))
)

title(main='Linear regression (m200 vs m100)')

dev.off()

# diagnostic of residuals
quartz()
par(mfrow=c(2, 2))
boxplot(residuals(regression), main='Boxplot of residuals')
qqnorm(residuals(regression))
plot(m100, residuals(regression), main='Residuals vs m100')
abline(h=0, lwd=2)
plot(fitted(regression), residuals(regression), main='Residuals vs fitted m200')
abline(h=0, lwd=2)

quartz()
par(mfrow=c(2, 2))
plot(regression)

dev.off()

detach(record)
graphics.off()



#_______________________________________________________________________________
##### Visualization of Multivariate Data
##### Example 1: dataset record (all the variables)

record <- read.table('record_mod.txt', header=T)
record

# Scatterplot
quartz()
pairs(record)  # or plot(record)

# Boxplot
boxplot(record, col='gold')

boxplot(log(record), col='gold')

# Starplot
quartz()
stars(record, col.stars=rep('gold',55))

# Radarplot
stars(record, draw.segments=T)

# Chernoff faces
quartz()
source('faces.R')
faces(record)

##### Example 2: cerebral aneurysm
aneurysm <- read.table('aneurysm.txt', header=T, sep=',')
head(aneurysm)
dim(aneurysm)

aneurysm.geometry <- aneurysm[, 1:4]
aneurysm.position <- factor(aneurysm[, 5])

head(aneurysm.geometry)

color.position <- ifelse(aneurysm.position == '1', 'red', 'blue')

attach(aneurysm.geometry)

quartz()
layout(cbind(c(1, 1), c(2, 3)), widths=c(2, 1), heights=c(1, 1))
plot(R1, R2, asp=1, col=color.position, pch=16)
hist(R1, prob=T, xlim=c(-10, 15))
hist(R2, prob=T, xlim=c(-10, 15))

quartz()
layout(cbind(c(1, 1), c(2, 3)), widths=c(2, 1), heights=c(1, 1))
plot(C1, C2, asp=1, col=color.position, pch=16)
hist(C1, prob=T, xlim=c(-5, 5))
hist(C2, prob=T, xlim=c(-5, 5))

graphics.off()

detach(aneurysm.geometry)

# some statistical indices
sapply(aneurysm.geometry, mean)
sapply(aneurysm.geometry, sd)
cov(aneurysm.geometry)
cor(aneurysm.geometry)

# Attention: rounded zeros!
round(sapply(aneurysm.geometry, mean), 1)
round(cov(aneurysm.geometry), 1)
round(cor(aneurysm.geometry), 1)

# Scatterplot
quartz()
pairs(aneurysm.geometry, col=color.position, pch=16)

# Boxplot
boxplot(aneurysm.geometry, col='gold')

# Stratified boxplots
quartz()
par(mfrow=c(1, 4))
boxplot(aneurysm.geometry$R1 ~ aneurysm.position, col=c('red', 'blue'), main='R1')
boxplot(aneurysm.geometry$R2 ~ aneurysm.position, col=c('red', 'blue'), main='R2')
boxplot(aneurysm.geometry$C1 ~ aneurysm.position, col=c('red', 'blue'), main='C1')
boxplot(aneurysm.geometry$C2 ~ aneurysm.position, col=c('red', 'blue'), main='C2')

# Stratified boxplots (same scale)
par(mfrow=c(1, 4))
boxplot(aneurysm.geometry$R1 ~ aneurysm.position, col=c('red', 'blue'), main='R1',
        ylim=range(aneurysm.geometry))
boxplot(aneurysm.geometry$R2 ~ aneurysm.position, col=c('red', 'blue'), main='R2',
        ylim=range(aneurysm.geometry))
boxplot(aneurysm.geometry$C1 ~ aneurysm.position, col=c('red', 'blue'), main='C1',
        ylim=range(aneurysm.geometry))
boxplot(aneurysm.geometry$C2 ~ aneurysm.position, col=c('red', 'blue'), main='C2',
        ylim=range(aneurysm.geometry))

dev.off()

# Chernoff faces
source('faces.R')
quartz()
faces(aneurysm.geometry)

# matplot
matplot(t(aneurysm.geometry), type='l')
matplot(t(aneurysm.geometry), type='l', col=color.position)


#_______________________________________________________________________________
##### Visualization of Categorical Data

district <- c('MI', 'MI', 'VA', 'BG', 'LO', 'LO', 'CR', 'Alt', 'CR', 'MI',  
              'Alt', 'CR', 'LO', 'VA', 'MI', 'Alt', 'LO', 'MI')
district <- factor(district, levels=c('MI', 'LO', 'BG', 'CR', 'VA', 'Alt'))
district

# Pie chart (no ordering of levels)
quartz()
pie(table(district), col=rainbow(length(levels(district)))) 

# 3D Pie chart (never use them!!)
library(plotrix)
quartz(width=14)
par(mfrow=c(1, 2))

pie3D(
  table(district)[1:length(levels(district))],
  labels = levels(district),
  explode = 0.1,
  main = "Pie Chart of Districts ",
  col = rainbow(length(levels(district)))
)

set.seed(020323)
shuffle = sample(1:length(levels(district)), size=length(levels(district)),
                 replace = F)
pie3D(table(district)[shuffle], labels=levels(district)[shuffle], explode=0.1,
      main="Pie Chart of Districts ", col=rainbow(length(levels(district)))[shuffle])

# Barplot (levels are ordered)
quartz()
barplot(table(district) / length(district))  

# or
plot(district)   # barplot of absolute frequences

# Remark: Some functions (e.g., the function plot()) may behave differently
# depending on the object it takes as input

is(district)[1]
plot(district)

# record is a data frame
is(record)[1]
plot(record) # scatterplot

# Remark 2: be careful to the scale of representation
quartz(width=11, height=5)
par(mfrow=c(1, 3))
barplot(table(district) / length(district), ylim=c(0, 1)); box()                       
barplot(table(district)/length(district),ylim=c(0, 10)); box() 
barplot(table(district)/length(district),ylim=c(0, 0.47)); box() 

graphics.off()

#_______________________________________________________________________________
### 3d plots (for instance, let's plot a bivariate Gaussian density)

x <- seq(-4, 4, 0.15)
y <- seq(-4, 4, 0.15)

# To build a function in R
gaussian <- function(x, y) {
  exp(-(x^2 + y^2 + x * y))
}

w <- matrix(NA, nrow = length(x), ncol=length(y))

# for
for(i in 1:length(x)) {
  for(j in 1:length(y)) {
    w[i, j] <- gaussian(x[i], y[j])
  }
}

# or
w <- outer(x, y, gaussian)
# help(outer)

quartz()
image(x, y, w)
contour(x, y, w, add=T)

persp(x, y, w, col='red')
persp(x, y, w, col='red', theta=30, phi=30, shade=.05, zlab='density')

# To download a package: 
# from RStudio: Tools -> Install Packages -> type PACKAGENAME 
#               and click install
# from R: Packages -> Install Packages -> Choose a CRAN mirror
#         (e.g., Italy (Milano)) -> Choose the package and click OK
# or install.packages('PACKAGENAME') in a R console

library(rgl)
options(rgl.printRglwidget = TRUE)
persp3d(x, y, w, col='red', alpha=1)
lines3d(x,x, gaussian(x,x), col='blue', lty=1)
lines3d(x,x, 0, col='blue', lty=2)

dev.off()

############################################################
####       More on graphical representation in R        ####
####           http://raw.densitydesign.org/            ####
####                     Ggobi                          ####
####               http://www.ggobi.org                 ####
############################################################


#_______________________________________________________________________________
##### Save plots

plot(rnorm(10), rnorm(10))

bmp(file="myplot.bmp")
plot(rnorm(10), rnorm(10))
dev.off()

jpeg(file="myplot.jpeg")
plot(rnorm(10), rnorm(10))
dev.off()

png(file="myplot.png")
plot(rnorm(10), rnorm(10))
dev.off()

pdf(file="myplot.pdf")
plot(rnorm(10), rnorm(10))
dev.off()

pdf(file="myplot2.pdf", onefile=T)
plot(rnorm(10), rnorm(10))
plot(rnorm(10), rnorm(10))
plot(rnorm(10), rnorm(10))
dev.off()
