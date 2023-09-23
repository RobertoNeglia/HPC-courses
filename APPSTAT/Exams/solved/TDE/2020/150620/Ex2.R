### Ex2


dataset = read.table("stoneflakes.txt", header = T)

n = dim(dataset)[1]

# Plot data
#2D
x11()
plot(dataset)

dataset.e <- dist(dataset, method='euclidean')

dataset.ew <- hclust(dataset.e, method='ward.D2')

x11()
plot(dataset.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')



k_cl = 3

x11()
plot(dataset.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ew, k=k_cl)


clusters <- cutree(dataset.ew, k=k_cl) # Compute the mean in each cluster
clust1 = dataset[which(clusters == 1), ]
clust2 = dataset[which(clusters == 2), ]
clust3 = dataset[which(clusters == 3), ]

group = rep(0,n)
group[which(clusters == 1)] = 1
group[which(clusters == 2)] = 2
group[which(clusters == 3)] = 3

group = factor(group, labels = c('1', '2', '3'))

n1 = dim(clust1)[1]
n2 = dim(clust2)[1]
n3 = dim(clust3)[1]

g= 3
p=2

fit <- manova(as.matrix(dataset) ~ group)  # the dataset should be a matrix-> for this we put as.matrix
s = summary.manova(fit,test="Wilks")
s

# verify 
# 1)  normality (multivariate) in each group (g tests)
load("mcshapiro.test.RData") 
P = c(mcshapiro.test(clust1)$pvalue, mcshapiro.test(clust2)$pvalue, mcshapiro.test(clust3)$pvalue)
# assume gaussianity
P
# 2) same covariance structure (= same covariance matrix Sigma)
S  <-  cov(dataset)
S1 <-  cov(clust1)   # S nel gruppo i = 1..g
S2 <-  cov(clust2)
S3 <-  cov(clust3)   # g righe

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)
round(S3,digits=1)

# to compare visually
x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))




### Via Bonferroni
alpha <- 0.1
p = 2
g = 3


k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)

W <- summary.manova(fit)$SS$Residuals   # save the residuals
m  <- sapply(dataset,mean)         # estimates mu
m1 <- sapply(clust1,mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(clust2,mean)    # estimates mu.2=mu+tau.2
m3 <- sapply(clust3,mean)    # estimates mu.3=mu+tau.3    # g righe

inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
inf13 <- m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
sup13 <- m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
inf23 <- m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
sup23 <- m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )

CI <- list(group1_group2=cbind(inf12, sup12), group1_group3=cbind(inf13, sup13), group2_group3=cbind(inf23, sup23))
CI







