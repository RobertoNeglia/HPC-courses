#_______________________________________________________________________________
##### Principal component analysis of the dataset 'running records'

runrec <- read.table('record_mod.txt', header=T)
n <- dim(runrec)[1]
p <- dim(runrec)[2]

# we make the units of measure homogeneous across the variables
runrec[,4:7] <- runrec[,4:7]*60

# Boxplot
x11()
boxplot(runrec)

# We observe that the variability increases non-linearly for increasing
# lengths; the record times of the marathon have much more variability
# than that of the other disciplines. This could significanlty influence 
# the PCA

# We perform the PCA on original data
pc.runrec <- princomp(runrec, scores=T)
pc.runrec
summary(pc.runrec)

# To obtain the rows of the summary:
# standard deviation of the components
pc.runrec$sd
# proportion of variance explained by each PC
pc.runrec$sd^2/sum(pc.runrec$sd^2)
# cumulative proportion of explained variance
cumsum(pc.runrec$sd^2)/sum(pc.runrec$sd^2)

# loadings (recall: coefficients of the linear combination of the original 
#           variables that defines each principal component)

load.rec <- pc.runrec$loadings
load.rec

load.rec[,1:7]

# graphical representation of the loadings of the first six principal components
x11()
par(mar = c(1,4,0,2), mfrow = c(6,1))
for(i in 1:6) barplot(load.rec[,i], ylim = c(-1, 1))

# Interpretation of the loadings:
# First PCs: long distances disciplines (I PC: variable 'Marathon')
# Last PCs: short distances disciplines 

# The loadings reflect the previous observation: the first PC is represented
# by the variable "Marathon", the second by the long distances, etc. The
# short distances appear in the last PCs

# Explained variance
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.runrec, las=2, main='Principal components', ylim=c(0,3.5e6))
barplot(sapply(runrec,sd)^2, las=2, main='Original Variables', ylim=c(0,3.5e6), ylab='Variances')
plot(cumsum(pc.runrec$sd^2)/sum(pc.runrec$sd^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(runrec),labels=1:ncol(runrec),las=2)

# The first PC (var. Marathon) explains more than 99.98% of the 
# total variability. This is due to the masking effect of that
# variable over the others

# scores
scores.runrec <- pc.runrec$scores  

layout(matrix(c(1,2),2))
boxplot(runrec, las=2, col='red', main='Original variables')
scores.runrec <- data.frame(scores.runrec)
boxplot(scores.runrec, las=2, col='red', main='Principal components')

# biplot
x11()
biplot(pc.runrec, scale=0, cex=.7)

graphics.off()

##### Principal component analysis of the dataset 'runningrecords',
##### but on the standardized variables

# We compute the standardized variables
runrec.sd <- scale(runrec)
runrec.sd <- data.frame(runrec.sd)

pc.runrec <- princomp(runrec.sd, scores=T)
pc.runrec
summary(pc.runrec)

# Explained variance
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.runrec, las=2, main='Principal Components', ylim=c(0,6))
abline(h=1, col='blue')
barplot(sapply(runrec.sd,sd)^2, las=2, main='Original Variables', ylim=c(0,6), ylab='Variances')
plot(cumsum(pc.runrec$sde^2)/sum(pc.runrec$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(runrec.sd),labels=1:ncol(runrec.sd),las=2)

# If we wanted to perform dimensionality reduction, we could keep
# 1 or 2 PCs

# loadings
load.rec <- pc.runrec$loadings
load.rec

x11()
par(mar = c(2,2,2,1), mfrow=c(3,1))
for(i in 1:3)barplot(load.rec[,i], ylim = c(-1, 1), main=paste('Loadings PC ',i,sep=''))

# Interpretation of the loadings:
# In this case, the first PC represents an average of the times of all the disciplines,
# taken with very similar (positive) weights. The sencond PC contrasts the short
# distances (m100, m200, m400) with the long distances (m800, m1500, m3000, Marathon)

# High PC1: long times in all the disciplines
# Low PC1: short times in all the disciplines 
# High PC2: long times in short distances, short times in long distances
# Low PC2: short times in short distances, long times in long distances

# scores
scores.runrec <- pc.runrec$scores
scores.runrec

x11()
plot(scores.runrec[,1:2])

x11()
layout(matrix(c(1,2),2))
boxplot(runrec.sd, las=2, col='red', main='Original variables')
scores.runrec <- data.frame(scores.runrec)
boxplot(scores.runrec, las=2, col='red', main='Principal components')

x11()
layout(matrix(c(1,2),1))
plot(runrec.sd[,'m100'],runrec.sd[,'Marathon'],type="n",xlab="m100",ylab="Marathon", asp=1)
text(runrec.sd[,'m100'],runrec.sd[,'Marathon'],dimnames(runrec)[[1]],cex=.7)
plot(scores.runrec[,1],scores.runrec[,2],type="n",xlab="pc1",ylab="pc2", asp=1)
text(scores.runrec[,1],scores.runrec[,2],dimnames(runrec)[[1]],cex=.7)

x11()
biplot(pc.runrec)

graphics.off()
