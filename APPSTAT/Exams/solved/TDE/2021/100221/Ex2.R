### Ex2


dataset = read.table("rice.txt", header = T)


# compute the dissimilarity matrix
dataset.e <- dist(dataset, method='euclidean')

dataset.ec <- hclust(dataset.e, method='complete')

k_cl = 3

x11()
plot(dataset.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.ec, k=k_cl)

cluster.ec <- cutree(dataset.ec, k=3)
table(cluster.ec)       # tabella che mi dice quanti samples sono in ogni cluster

x11()
plot(dataset, col= cluster.ec + 1, pch=19)

clust1 = dataset[which(cluster.ec == 1), ]
clust2 = dataset[which(cluster.ec == 2), ]
clust3 = dataset[which(cluster.ec == 3), ]

mean1 = colMeans(clust1)
mean2 = sapply(clust2, mean)
mean3 = sapply(clust3, mean)

m = rbind(clust1 = mean1, clust2 = mean2, clust3 = mean3)
m

# b) 
# coph

coph.ec <- cophenetic(dataset.ec)
ec <- cor(dataset.e, coph.ec)
ec

# k-means

n = dim(dataset)[1]

k = 3       # number of clusters 

result.k <- kmeans(dataset, centers=k) # Centers: fixed number of clusters

x11()
plot(dataset, col = result.k$cluster+1)  # we have one color for each cluster

result.k$centers      # centers of the clusters
result.k$size 



# c) 
x.mean = colMeans(clust1)[1]
x.cov = var(clust1[,1])
n = dim(clust1)[1]

shapiro.test(clust1[,1])

## Bonf
alpha = 0.05
k <- 2
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = x.mean - cfr.t*sqrt((x.cov)/n),
            center = x.mean, 
            sup = x.mean + cfr.t*sqrt((x.cov)/n))
Bf

Bf_var <- cbind(inf= (x.cov*(n-1) )/ qchisq(1 - alpha/(2*k), n-1),
                center=x.cov,
                sup=(x.cov*(n-1)) / qchisq(alpha/(2*k), n-1))
Bf_var





















#### cambio linkage

# compute the dissimilarity matrix
dataset.e <- dist(dataset, method='euclidean')

dataset.es <- hclust(dataset.e, method='single')

k_cl = 3

x11()
plot(dataset.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(dataset.es, k=k_cl)

cluster.es <- cutree(dataset.es, k=3)
table(cluster.es)       # tabella che mi dice quanti samples sono in ogni cluster

x11()
plot(dataset, col= cluster.es + 1, pch=19)

coph.es <- cophenetic(dataset.es)
es <- cor(dataset.e, coph.es)
es


