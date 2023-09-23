## Ex1

load('mcshapiro.test.RData')

dataset = read.table("pinnanobilis.txt", header = T)
# only numeric

n = dim(dataset)[1]

# Plot data
#2D
x11()
plot(dataset)

dataset.e <- dist(dataset, method='euclidean')
dataset.ec <- hclust(dataset.e, method='complete')


k_cl = 2
x11()
plot(dataset.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ec, k=k_cl)

cluster.ec <- cutree(dataset.ec, k=k_cl) 





# b) 
# Verify assumptions
g=2

P = c(mcshapiro.test(clust1)$pvalue, mcshapiro.test(clust2)$pvalue)
P

S1 = cov(clust1)
S2 = cov(clust2)

x11()
par(mfrow=c(1,g))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))

round(S1,digits=1)
round(S2,digits=1)

group = rep(0,n)
group[which(cluster.ec == 1)] = 1
group[which(cluster.ec == 2)] = 2

group = factor(group, labels=c('1', '2'))

fit <- manova(as.matrix(dataset) ~ group)  # the dataset should be a matrix-> for this we put as.matrix
s = summary.manova(fit,test="Wilks")
s




# strong diff
clust1 = dataset[which(cluster.ec == 1), ]
clust2 = dataset[which(cluster.ec == 2), ]

m  <- sapply(dataset,mean)         # estimates mu
m1 <- sapply(clust1,mean)    # estimates mu.1=mu+tau.1
m2 <- sapply(clust2,mean)    # estimates mu.2=mu+tau.2


tau1 = m1-m
tau2 = m2-m

DF <- fit$df # fit è il modello finale

Spooled <- sum(fit$res^2)/DF
Spooled  # estimate for the variance 

### Via Bonferroni
alpha <- 0.1
p= 2 
k <- p*g*(g-1)/2
qT <- qt(1-alpha/(2*k), n-g)

W <- summary.manova(fit)$SS$Residuals   # save the residuals
m  <- sapply(dataset,mean)         # estimates mu

n1 = dim(clust1)[1]
n2 = dim(clust2)[1]
inf12 <- m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )


CI <- list(group1_group2=cbind(inf12, sup12))
CI





