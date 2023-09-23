## Ex2

dataset = read.table("streaming.txt", header = T)
# only numeric

n = dim(dataset)[1]

# Plot data
#2D
x11()
plot(dataset)


dataset.e <- dist(dataset, method='euclidean')
dataset.es <- hclust(dataset.e, method='single')

x11()
plot(dataset.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')


k_cl = 3
cluster.es <- cutree(dataset.es, k=k_cl) 
cluster.es
cluster.es <- cutree(dataset.es, k= k_cl)
table(cluster.es)       # tabella che mi dice quanti samples sono in ogni cluster
which(cluster.ec==2)    # estraggo gli indici delle righe nel cluster 2



clust1 = dataset[which(cluster.es == 1), ]
clust2 = dataset[which(cluster.es == 2), ]
clust3 = dataset[which(cluster.es == 3), ]


mean1 = colMeans(clust1)
mean2 = colMeans(clust2)
mean3 = colMeans(clust3)

mean1
mean2
mean3

m = rbind(clust1 = mean1, clust2 = mean2, clust3 = mean3)
m

coph.es <- cophenetic(dataset.es)

es <- cor(dataset.e, coph.es)
es

## Bonf

# clust 1 

cov1    <- cov(clust1)
n = dim(clust1)[1]

k <- 2
alpha = 0.05
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = mean1 - cfr.t*sqrt(diag(cov1)/n),
            center = mean1, 
            sup = mean1 + cfr.t*sqrt(diag(cov1)/n))
Bf


# clust 2 
cov2    <- cov(clust2)
n = dim(clust2)[1]

k <- 2
alpha = 0.05
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = mean2 - cfr.t*sqrt(diag(cov2)/n),
            center = mean2, 
            sup = mean2 + cfr.t*sqrt(diag(cov2)/n))
Bf


# clust 3
cov3    <- cov(clust3)
n = dim(clust3)[1]

k <- 2
alpha = 0.05
cfr.t <- qt(1 - alpha/(k*2), n-1)

Bf <- cbind(inf = mean3 - cfr.t*sqrt(diag(cov3)/n),
            center = mean3, 
            sup = mean3 + cfr.t*sqrt(diag(cov3)/n))
Bf





