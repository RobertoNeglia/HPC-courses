#############################################################################
#                                Example PCA                                #
#############################################################################

# Along the ringroads of Milan four control units measure the concentration  
# of the pollutant NO in the air.
# The measures collected during the last year are reported in the file NO.txt.
# Perform a principal component analysis of the available data. In particular:
# (a) Compute the loadings
# (b) Compute the variances along the PCs
# (c) Comment and interpret the results at points (a) and (b)
# (d) On the 3rd July the control unites registered the values (13, 10, 11, 13).
#     Compute the corresponding scores

#############################################################################

# Import the data
NO <- read.table('NO.txt', header=T)
NO

dim(NO)
dimnames(NO)

NO <- data.frame(NO)
var.names <- c("I Control Unit","II Control Unit","III Control Unit","IV Control Unit")
dimnames(NO)[[2]] <- var.names

## DATA EXPLORATION ##
# Scatter plot 
pairs(NO, col=rainbow(dim(NO)[1]), pch=16, main='Scatter plot')

## COMMENT ## 
# Positive correlation between I-II, II-III, I-III, low correlation between I-IV,
# II-IV, very low and negative correlation between III-IV.

# This is confirmed by quantitative analyses
# we compute the sample mean, sample covariance matrix and sample correlation matrix
M <- sapply(NO,mean)
M
S <- cov(NO)
S
R <- cor(NO)
R

# Boxplot
x11()
boxplot(NO, las=1, col='red', main='Boxplot',grid=T)

# Matplot + boxplot
x11()
matplot(t(NO), type='l', axes=F)
box()
boxplot(NO, add=T, boxwex=0.1, col='red')

## COMMENT ## 
# The variability of the variables are similar, medians of the II and IV control
# units seem lower (as well as the means); almost symmetric distributions, with 
# little asymmetry in the measurements of the III control unit. 
# Presence of outliers in all the variables. In particular, there is an outlying
# measurement (NO[5,]) for the CUs I-II-III, not for the CU IV.

###############################################################################
# (a) Compute the loadings                                                    #
# (b) Compute the variances along the PCs                                     #
# (c) Comment and interpret the results at points (a) and (b)                 #
###############################################################################

## PRINCIPAL COMPONENT ANALYSIS ##
## COMMENT ## 
# The original variability along the 3 variables is similar; the units of 
# measure of the variables are homogeneous. We thus perform a PCA based on S
pca.NO <- princomp(NO, scores=T)
pca.NO
summary(pca.NO)

x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(pca.NO$sdev^2, las=2, main='Principal components', ylim=c(0,10), ylab='Variances')
barplot(sapply(NO,sd)^2, las=2, main='Original variables', ylim=c(0,5), ylab='Variances')
plot(cumsum(pca.NO$sdev^2)/sum(pca.NO$sde^2), type='b', axes=F, xlab='number of components', 
     ylab='contribution to the total variability', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(NO),labels=1:ncol(NO),las=2)

## COMMENT ## 
# Dimensionality reduction: the first two components explain 73.6%<80% of the total
# variability. However, we do not see an elbow in the proportion of explained variability
# in correspondence of the III or IV PC. It could be thus sufficient to analyse the 
# sample through the first 2 PCs.
# Note. The results of the PCA on the correlation matrix would not give very different results
NO.sd <- scale(NO)
pca.NO.std <- princomp(NO.sd, scores=T)
pca.NO.std

summary(pca.NO.std) # PCA on R
summary(pca.NO)     # PCA on S

# SCORES AND LOADINGS (PCA on S)#
# Scores
scores.NO <- pca.NO$scores
scores.NO
layout(matrix(c(1,2),1,2))
boxplot(NO, las=2, col='red', main='Variabili originarie')
scores.NO <- data.frame(scores.NO)
boxplot(scores.NO, las=2, col='red', main='Componenti principali')
# pairs(scores.NO, col=rainbow(dim(scores.NO)[1]), pch=16)
# cor(scores.NO)

## COMMENT ## 
# The variability of the first two PCs (looking at the size of the boxes)
# is greater than that of the components 3 and 4. These have very similar
# variability: we should avoid to keep 3 PCs, let's keep 2 or 4 PCs 
# instead.
# The proportion of variability explained by the second PC is influenced
# by the presence of 4 outliers. These represents 4 out of the 5 outliers
# of the IV CU (NO[3,4], NO[59,4], NO[125,4], NO[308,4]). The outlier of
# PC 1 corresponds to the outlying observation mentioned before and to the
# measurement NO[217,]

# We plot the outlying data for the first 2 PCs
color.outliers <- NULL
for (i in 1:365){
  color.outliers[i] <- 'blue'
}
color.outliers[5] <- 'red'
color.outliers[217] <- 'red'
color.outliers[3] <- 'green'
color.outliers[125] <- 'green'
color.outliers[59] <- 'green'
color.outliers[308] <- 'green'
pairs(NO, col=color.outliers, pch=16, main='Scatter plot - Outliers')

# Loadings
load.NO    <- pca.NO$loadings
load.NO

x11()
par(mar = c(1,4,0,2), mfrow = c(4,1))
for(i in 1:4)
  barplot(load.NO[,i], ylim = c(-1, 1))

## COMMENT ## 
# The first PC is a weighted mean (with very similar weights) of the measurements
# recorded by the CU I-II-III; the second PC corresponds to the IV CU. The last
# two PCs do not have a clear interpretations, although they are combinations of
# the CU I-II-III.
# Note: it makes sense that the outliers of the PC1 correspond to the ouliers of the
# CU I-II-III, and the outlier of the PC2 are those of CU IV.

# We can geometrically represent in a 3D space the dimensionality reduction, since 
# the loadings of the CU IV are almost 0 in that direction

library(rgl)

M <- sapply(NO[,1:3],mean)
S <- cov(NO[,1:3])

open3d()
points3d(NO[,1:3], col=rainbow(dim(NO)[1]), asp=1, size=5)
axes3d()
plot3d(ellipse3d(S, centre=M, level= 9/10), alpha=0.25, add = TRUE)
title3d(xlab="I CU",ylab="II CU",zlab="III CU")
pca.NO1 <- NULL
for(i in 1:365)
  pca.NO1 <- rbind(pca.NO1, pca.NO$loadings[1:3,1]*pca.NO$scores[i,1] + M)
points3d(pca.NO1, col='black', size=6)
lines3d(rbind(M + 4*pca.NO$sdev[1] * pca.NO$loadings[1:3,1], M - 4*pca.NO$sdev[1] * pca.NO$loadings[1:3,1]), col='black') 
color=rainbow(dim(NO)[1])
for(i in 1:365)
  lines3d(rbind(NO[i,], pca.NO1[i,]),col=color[i])

open3d()
points3d(NO[,1:3], col=color.outliers, asp=1, size=5)
axes3d()
plot3d(ellipse3d(S, centre=M, level= 9/10), alpha=0.25, add = TRUE)
title3d(xlab="I CU",ylab="II CU",zlab="III CU")
pca.NO1 <- NULL
for(i in 1:365)
  pca.NO1 <- rbind(pca.NO1, pca.NO$loadings[1:3,1]*pca.NO$scores[i,1] + M)
points3d(pca.NO1, col='black', size=5)

# We separately represent the PC2, and compare it with the CU IV variable
# (in this case, the loadings corresponding to the CUs I-II-IV are almost zero)

x11()
plot(NO[,4],scores.NO[,2],col=rainbow(n),pch=19,xlab='IV CU',ylab='Comp.2')

##TRY TO REMOVE THE OUTLIERS##
NO.outliers <- NO[which(color.outliers=='blue'),]
dim(NO.outliers)

pca.NO.outliers <- princomp(NO.outliers, scores=T)
pca.NO.outliers
summary(pca.NO.outliers)

# SCORES AND LOADINGS #
# Scores
scores.NO.outliers <- pca.NO.outliers$scores
scores.NO.outliers
layout(matrix(c(1,2),1,2))
boxplot(NO.outliers, las=2, col='red', main='Original Variables')
scores.NO.outliers <- data.frame(scores.NO.outliers)
boxplot(scores.NO.outliers, las=2, col='red', main='Principal Components')

# Loadings
load.NO.outliers    <- pca.NO.outliers$loadings
load.NO.outliers

x11()
par(mar = c(1,4,0,2), mfrow = c(4,1))
for(i in 1:4)
  barplot(load.NO.outliers[,i], ylim = c(-1, 1))

## COMMENT ## 
# The results obtained by removing the outliers are not very different from
# those obtained before. Both the analyses lead to the same choice to keep
# 2 PCs for the purpose of dimensionality reduction

# We finally note that if we keep the first 2 PCs we get that:
# - the first PC is roughly the mean of the observations of the CUs I-II-III
# - the second PC corresponds to the IV CU.
# From the application viewpoint, this could suggest a malfunctioning of the
# IV CU, whose measurements are actually uncorrelated with those of the 
# CU I-II-III, although they are geographically near.

#################################################################################
# (d) On the 3rd July the control unites registered the values (13, 10, 11, 13).#
#     Compute the corresponding scores                                          #
#################################################################################

data.3jul <- c(13, 10, 11, 13)
scores.3jul <- t(pca.NO$loadings)%*%(data.3jul-colMeans(NO))
scores.3jul

x11()
pairs(data.frame(rbind(NO,data.3jul)), col=c(rep(1,n),2), pch=16, main='Scatter plot - Data 3rd July')

x11()
plot(scores.NO[,1],scores.NO[,2],col='grey',pch=19,xlab='Comp.1',ylab='Comp.2')
points(scores.3jul[1],scores.3jul[2],col='black',pch=19) 
