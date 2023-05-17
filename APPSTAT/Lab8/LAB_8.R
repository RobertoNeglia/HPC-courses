###--------------------###
### LAB 8 (28/04/2023) ###
###--------------------###

### TOPICS:
### Hierarchical clustering
### K-means clustering
### DBSCAN
### Exercises
### Multidimensional Scaling

quartz.options(height=6.5, width=11, reset=FALSE)
options(rgl.printRglwidget = TRUE)

library(mvtnorm)
library(MVN)
library(rgl)
library(car)
library(dbscan)
library(cluster)
library(fields)

#_______________________________________________________________________________
##### Hierarchical Clustering

##### Example 1: iris dataset
#####-------------------------

# let's forget the labels (we perform a cluster analysis, not a discriminant one)

species.name <- iris[,5]
iris4        <- iris[,1:4]

quartz()
pairs(iris4, pch=19)

dev.off()

# compute the dissimilarity matrix of the data
# we choose the Euclidean metric (and then we look at other metrics)
help(dist)
iris.e <- dist(iris4, method='euclidean')

quartz()
image(1:150,1:150,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j')

# with other metrics:
iris.m <- dist(iris4, method='manhattan')
iris.c <- dist(iris4, method='canberra')

quartz()
par(mfrow=c(1,3))
image(1:150,1:150,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

# actually, the data are never ordered according to (unknown) labels
misc <- sample(150)
iris4 <- iris4[misc,]

iris.e <- dist(iris4, method='euclidean')
iris.m <- dist(iris4, method='manhattan')
iris.c <- dist(iris4, method='canberra')

quartz()
par(mfrow=c(1,3))
image(1:150,1:150,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:150,1:150,as.matrix(iris.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

graphics.off()

# we now aim to perform hierarchical clustering of the dataset iris
# we limit to Euclidean distance 

# Command hclust()
help(hclust)

iris.es <- hclust(iris.e, method='single')
iris.ea <- hclust(iris.e, method='average')
iris.ec <- hclust(iris.e, method='complete')

# if we want more detailed information on euclidean-complete
# clustering:
names(iris.ec)
iris.ec$merge  # order of aggregation of statistical units / clusters
iris.ec$height # distance at which we have aggregations
iris.ec$order  # ordering that allows to avoid intersections in the dendrogram

# plot of the dendrograms
quartz()
par(mfrow=c(1,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

dev.off()

# plot dendrograms (2 clusters)
quartz()
par(mfrow=c(1,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.es, k=2)
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ec, k=2)
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ea, k=2)

dev.off()

# plot dendrograms (3 clusters)
quartz()
par(mfrow=c(1,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
#rect.hclust(iris.es, k=2)
rect.hclust(iris.es, k=3)
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
#rect.hclust(iris.ec, k=2)
rect.hclust(iris.ec, k=3)
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
#rect.hclust(iris.ea, k=2)
rect.hclust(iris.ea, k=3)

dev.off()

# How to cut a dendrogram?
# We generate vectors of labels through the command cutree()
help(cutree)

# Fix k=2 clusters:
cluster.ec <- cutree(iris.ec, k=2) # euclidean-complete:
cluster.ec

cluster.es <- cutree(iris.es, k=2) # euclidean-single
cluster.ea <- cutree(iris.ea, k=2) # euclidean-average

# interpret the clusters
table(label.true = species.name[misc], label.cluster = cluster.es)
table(label.true = species.name[misc], label.cluster = cluster.ec)
table(label.true = species.name[misc], label.cluster = cluster.ea)

quartz()
plot(iris4, col=ifelse(cluster.es==1,'red','blue'), pch=19)
quartz()
plot(iris4, col=ifelse(cluster.ec==1,'red','blue'), pch=19)
quartz()
plot(iris4, col=ifelse(cluster.ea==1,'red','blue'), pch=19)

graphics.off()

# Let's give a mark to the algorithms: did they aggregate coherently with
# the dissimilarity matrix or not?

# compute the cophenetic matrices 
coph.es <- cophenetic(iris.es)
coph.ec <- cophenetic(iris.ec)
coph.ea <- cophenetic(iris.ea)

# compare with dissimilarity matrix (Euclidean distance)
quartz()
layout(rbind(c(0,1,0),c(2,3,4)))
image(as.matrix(iris.e), main='Euclidean', asp=1 )
image(as.matrix(coph.es), main='Single', asp=1 )
image(as.matrix(coph.ec), main='Complete', asp=1 )
image(as.matrix(coph.ea), main='Average', asp=1 )

dev.off()

# compute cophenetic coefficients
es <- cor(iris.e, coph.es)
ec <- cor(iris.e, coph.ec)
ea <- cor(iris.e, coph.ea)

c("Eucl-Single"=es,"Eucl-Compl."=ec,"Eucl-Ave."=ea)

# Exercise: repeat the analysis with other metrics

#_______________________________________________________________________________
##### Example 2: simulated data
#####---------------------------------------------------------------

### Univariate case
###-----------------
set.seed(123)
# x : vector of NON clustered data
x <- 0:124/124 + rnorm(125, sd=0.01)
# y : vector of clustered data (5 clusters)
y <- c(rnorm(25, mean=0, sd=0.01), rnorm(25, mean=0.25, sd=0.01),
       rnorm(25, mean=0.5, sd=0.01), rnorm(25, mean=0.75, sd=0.01),
       rnorm(25, mean=1, sd=0.01))

x <- sample(x)
y <- sample(y)

quartz()
par(mfrow=c(1,2))
plot(rep(0,125),x, main='data: no clust',xlab='')
plot(rep(0,125),y, main='data: clust',xlab='')

dx <- dist(x)
dy <- dist(y)
hcx<- hclust(dx, method='single')
hcy<- hclust(dy, method='single')

quartz()
par(mfrow=c(1,2))
plot(hcx, labels=F, cex=0.5, hang=-0.1, xlab='', sub='x')
plot(hcy, labels=F, cex=0.5, hang=-0.1, xlab='', sub='y')

quartz()
par(mfrow=c(2,2))
image(as.matrix(dx), asp=1, main='not ordered x' )
image(as.matrix(dy), asp=1, main='not ordered y' )
image(as.matrix(dx)[hcx$order,hcx$order], asp=1, main='reordered x' )
image(as.matrix(dy)[hcy$order,hcy$order], asp=1, main='reordered y' )

graphics.off()


### Bivariate case (example of chaining effect)
###--------------------------------------------
p <- 2
n <- 100

mu1 <- c(0,1)
mu2 <- c(5,1.2)
sig <- diag(rep(1,p))

set.seed(1)
X1 <- rmvnorm(n, mu1, sig)
X2 <- rmvnorm(n, mu2, sig)

X <- rbind(X1, X2)
# If we knew the labels:
quartz()
plot(X, xlab='Var 1', ylab='Var 2', col=rep(c('red','blue'),each=100), asp=1, pch=19)

# How we actually see the data:
plot(X, xlab='Var 1', ylab='Var 2', asp=1, pch=19)

dev.off()

# let's compare the clustering results with Euclidean distance and three 
# types of linkage
x.d <- dist(X, method = 'euclidean')

x.es <- hclust(x.d, method='single')
x.ea <- hclust(x.d, method='average')
x.ec <- hclust(x.d, method='complete')

quartz()
par(mfrow=c(1,3))
plot(x.es, main = 'Single linkage'  , hang=-0.1, xlab='', labels=F, sub='')
plot(x.ea, main = 'Average linkage' , hang=-0.1, xlab='', labels=F, sub='')
plot(x.ec, main = 'Complete linkage', hang=-0.1, xlab='', labels=F, sub='')


# let's cut the tree as to get 2 clusters
cluster.es <- cutree(x.es, k = 2)
cluster.ea <- cutree(x.ea, k = 2)
cluster.ec <- cutree(x.ec, k = 2)

quartz()
par(mfrow=c(1,3))

plot(X, xlab='Var 1', ylab='Var 2', main = 'Single linkage', col=ifelse(cluster.es==1,'red','blue'), pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Average linkage', col=ifelse(cluster.ea==1,'red','blue'), pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Complete linkage', col=ifelse(cluster.ec==1,'red','blue'), pch=16, asp=1)

graphics.off()

### Bivariate case (example of ellipsoidal clusters)
###-------------------------------------------------

p <- 2
n <- 100

mu1 <- c(0,1)
mu2 <- c(6.5,1)

e1 <- c(1,1)
e2 <- c(-1,1)
sig <- 5*cbind(e1)%*%rbind(e1)+.1*cbind(e2)%*%rbind(e2)

set.seed(2)
X1 <- rmvnorm(n, mu1, sig)
X2 <- rmvnorm(n, mu2, sig)

X <- rbind(X1, X2)

# If we knew the labels:
quartz()
plot(X, xlab='Var 1', ylab='Var 2', asp=1, col=rep(c('red','blue'),each=100), pch=19)

# How we actually see the data:
plot(X, xlab='Var 1', ylab='Var 2', asp=1, pch=19)

dev.off()

## compare results with different linkage, when Euclidean distances is used

x.d <- dist(X, method = 'euclidean')

x.es <- hclust(x.d, method='single')
x.ea <- hclust(x.d, method='average')
x.ec <- hclust(x.d, method='complete')

quartz()
par(mfrow=c(1,3))

plot(x.es, main = 'Single linkage', ylab='Euclidean distance', hang=-0.1, xlab='', labels=F, sub='')
plot(x.ea, main = 'Average linkage', hang=-0.1, xlab='', labels=F, sub='')
plot(x.ec, main = 'Complete linkage', hang=-0.1, xlab='', labels=F, sub='')


# let's cut the tree to get 2 clusters
quartz()
par(mfrow=c(1,3))

plot(x.es, main = 'Single linkage', ylab='Euclidean distance', hang=-0.1, xlab='', labels=F, sub='')
rect.hclust(x.es, k=2)
plot(x.ea, main = 'Average linkage', hang=-0.1, xlab='', labels=F, sub='')
rect.hclust(x.ea, k=2)
plot(x.ec, main = 'Complete linkage', hang=-0.1, xlab='', labels=F, sub='')
rect.hclust(x.ec, k=2)

cluster.es <- cutree(x.es, k = 2)
cluster.ea <- cutree(x.ea, k = 2)
cluster.ec <- cutree(x.ec, k = 2)

quartz()
par(mfrow=c(1,3))

plot(X, xlab='Var 1', ylab='Var 2', main = 'Euclidean, Single linkage', col=cluster.es+1, pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Euclidean, Average linkage', col=cluster.ea+1, pch=16, asp=1)
plot(X, xlab='Var 1', ylab='Var 2', main = 'Euclidean, Complete linkage', col=cluster.ec+1, pch=16, asp=1)

graphics.off()

# Change the mean mu2
# Example:
# mu2 <- c(4.5,1.2)


#_______________________________________________________________________________
##### Example 3: earthquakes dataset 
#####-------------------------------

help(quakes)

head(quakes)
dim(quakes)

Q <- cbind(quakes[,1:2], depth = -quakes[,3]/100)
head(Q)

plot3d(Q, size=3, col='orange', aspect = F)

# dissimilarity matrix (Euclidean metric)
d <- dist(Q)

quartz()
image(as.matrix(d))

quartz()
par(mfrow=c(2,2))

clusts <- hclust(d, method='single')
plot(clusts, hang=-0.1, labels=FALSE, main='single', xlab='', sub='')
# rect.hclust(clusts, k=2)
# rect.hclust(clusts, k=3)

clusta <- hclust(d, method='average')
plot(clusta, hang=-0.1, labels=FALSE, main='average', xlab='', sub='')
# rect.hclust(clusta, k=2)
# rect.hclust(clusta, k=3)

clustc <- hclust(d, method='complete')
plot(clustc, hang=-0.1, labels=FALSE, main='complete', xlab='', sub='')
# rect.hclust(clustc, k=2)
# rect.hclust(clustc, k=3)

clustw <- hclust(d, method='ward.D2')
plot(clustw, hang=-0.1, labels=FALSE, main='ward', xlab='', sub='')
# rect.hclust(clustw, k=2)
# rect.hclust(clustw, k=3)

## Ward-Linkage (see J-W p. 692-693): 
# Ward considered hierarchical clustering procedures based on minimizing the
# 'loss of information' from joining two groups. This method is usually implemented
# with loss of information taken to be an increase in an error sum of squares 
# criterion, ESS. First for a given cluster k, let ESS[k] be the sum of the
# squared deviations of every item in the cluster from the cluster mean (centroid).
# (ESS[k]=sum_{x.j in cluster k} t(x.j-x.mean[k])%*%(x.j-x.mean[k]), where 
# x.mean[k] is the centroid of cluster k).
# If there are currently K clusterts, define ESS as the sum of the ESS[k] 
# (ESS=ESS[1]+ESS[2]+...+ESS[K]). At each step in the analysis, the union of
# every possible pair of clusters is considered, and the two clusters whose 
# combination results in the smallest increase in ESS (minimum loss of information)
# are joined. 
# Initially, each cluster consists of a single item and, if there are N items, 
# ESS[k]=0, k=1,2,...,N, so ESS=0. At the other extreme, when all the clusters are 
# combined in a single group of N items, the value of ESS is given by
# ESS=sum_j(t(x.j-x.mean)%*%(x.j-x.mean)), where x.j is the multivariate measurement
# associated with the jth item and x.mean is the mean of all items.

open3d()

# single linkage
clusters <- cutree(clusts, 2)
plot3d(Q, size=3, col=clusters+1, aspect = F)

clusters <- cutree(clusts, 3)
plot3d(Q, size=3, col=clusters+1, aspect = F) 

# average linkage
clustera <- cutree(clusta, 2)
plot3d(Q, size=3, col=clustera+1, aspect = F) 

clustera <- cutree(clusta, 3)
plot3d(Q, size=3, col=clustera+1, aspect = F) 

# complete linkage
clusterc <- cutree(clustc, 2)
plot3d(Q, size=3, col=clusterc+1, aspect = F) 

clusterc <- cutree(clustc, 3)
plot3d(Q, size=3, col=clusterc+1, aspect = F) 

# ward linkage
clusterw <- cutree(clustw, 2)
plot3d(Q, size=3, col=clusterw+1, aspect = F) 

clusterw <- cutree(clustw, 3)
plot3d(Q, size=3, col=clusterw+1, aspect = F)


#_______________________________________________________________________________
##### K-means method

#################################################################
# Recall K-means algorithm:

# simulated data
n <- 100
set.seed(1)
x <- matrix(rnorm(n*2), ncol=2)
x[1:(n/2),1] <- x[1:(n/2),1]+2
x[1:(n/2),2] <- x[1:(n/2),2]-2

quartz()
plot(x,pch=20,cex=2,xlab='x1',ylab='x2')

# k-means algorithm

k <- 2
cluster <- sample(1:2, n, replace=TRUE)
iter.max <- 3

colplot <- c('royalblue','red')
colpoints <- c('blue4','red4')

quartz()
par(mfrow = c(iter.max,3))

for(i in 1:iter.max)
{
  C <- NULL
  for(l in 1:k)
    C <- rbind(C, colMeans(x[cluster == l,]))
  
  plot(x, col = colplot[cluster],pch=19)
  line <- readline()
  
  points(C, col = colpoints, pch = 4, cex = 2, lwd = 2)
  line <- readline()
  
  plot(x, col = 'grey',pch=19)
  points(C, col = colpoints, pch = 4, cex = 2, lwd = 2)
  line <- readline()
  
  QC <- rbind(C, x)
  Dist <- as.matrix(dist(QC, method = 'euclidean'))[(k+1):(k+n),1:k]
  for(j in 1:n)
    cluster[j] <- which.min(Dist[j,])
  
  plot(x, col = colplot[cluster],pch=19)
  points(C, col = colpoints, pch = 4, cex = 2, lwd = 2)
  line <- readline()
  
}

graphics.off() 

#################################################################

# back to the earthquakes dataset

### in automatic, command kmeans()
help(kmeans)

result.k <- kmeans(Q, centers=2) # Centers: fixed number of clusters

names(result.k)

result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares within cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimension of the clusters

quartz()
plot(Q, col = result.k$cluster+1)

open3d()
plot3d(Q, size=3, col=result.k$cluster+1, aspect = F) 
points3d(result.k$centers,size=10)

### How to choose k:
### 1) evaluate the variability between the groups with respect to 
###   the variability withing the groups
### 2) evaluate the result of hierarchical clustering (not recommended,
###    quite computationally expensive)

# (we've just seen method 2) and suggested k = 2)

# method 1)
b <- NULL
w <- NULL
for(k in 1:10){
  
  result.k <- kmeans(Q, k)
  w <- c(w, sum(result.k$wit))
  b <- c(b, result.k$bet)
  
}

quartz()
matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
lines(1:10, w/(w+b), type='b', lwd=2)

# this method seems to suggest k = 2 or 4
# let's try also k=4:
result.k <- kmeans(Q, 4)

quartz()
plot(Q, col = result.k$cluster+1)

open3d()
plot3d(Q, size=3, col=result.k$cluster+1, aspect = F) 
points3d(result.k$centers, size = 10)

#_______________________________________________________________________________
##### DBSCAN

# Simulate the data
set.seed(2)
n <- 400
x <- cbind(x = runif(4) + rnorm(n, sd = 0.1), y = runif(4) + rnorm(n, sd = 0.1))
true_clusters <- rep(1:4, time = 100)
quartz()
plot(x, col = true_clusters, pch = 19)

# Choice of hyperparameters for DBSCAN
# Rule of thumb, minPts = dimensionality + 1 = 3 here

# How to choose eps from minPts?
# Plot of the distances to the minPts nearest neighbor
kNNdistplot(x, k = 3)
# Taking eps = 0.05 seems to be a good threshold
abline(h = 0.05, col = "red", lty = 2)

# Run the dbscan
dbs <- dbscan(x, eps = 0.05, minPts = 3)
dbs

# Plot of the resulting clustering
plot(x, col = dbs$cluster + 1L, pch=19)


# How to tune the algorithm and find the "best" eps and minPts?
# We can consider grid search BUT we need an objective function

# In the absence of the true labels (which would be the case in practice),
# we have to use an "intrinsic" clustering quality measure

# Silhouette score (from the package "cluster")
help(silhouette)

# Let's compute the silhouette score on the clustering performed before
# WARNING (specific to DBSCAN): We need to remove the noise points as they do
# not belong to a cluster, before computing the silhouette score
clustered_index <- which(dbs$cluster != 0) # Index of non noise points
clustered_points <- x[clustered_index] # only clustered points
clustered_labels <- dbs$cluster[clustered_index] # corresponding labels

sil <- silhouette(clustered_labels, dist(clustered_points))
summary(sil)

sil_score <- function(labels, dist) {
  # Compute the average of the silhouette widths
  sil <- silhouette(labels, dist)
  sil_widths <- sil[,"sil_width"]
  mean(sil_widths)
}

sil_score(clustered_labels, dist(clustered_points))

# Grid Search
minPts_grid <- 1:20
eps_grid <- seq(0.01, 0.2, by = 0.01)

max_share_noise <- 0.2

dbscan_perf <- function(minPts, eps) {
  # Compute the silhouette score resulting from dbscan clustering
  dbs <- dbscan(x, eps, minPts) # Run dbscan
  
  clustered_index <- which(dbs$cluster != 0) # Index of non noise points
  clustered_points <- x[clustered_index] # only clustered points
  clustered_labels <- dbs$cluster[clustered_index] # corresponding labels
  nb_clusters <- length(unique(clustered_labels))
  
  if ((nb_clusters > 1 & nb_clusters < n) & (length(which(dbs$cluster == 0))/n < max_share_noise)) { 
    # Silhouette score is defined only if 2 <= nb_clusters <= n-1
    sil_score(clustered_labels, dist(clustered_points))
  }
  
  else {
    # otherwise we return 0 which would be the approx. value of the silhouette
    # score if the clusters were completely overlapping
    0
  }
}

# We compute the silhouette score for all combinations of minPts and eps
perf_grid <- outer(minPts_grid, eps_grid, FUN = Vectorize(dbscan_perf))
dimnames(perf_grid) <- list(minPts_grid, eps_grid)

# Histogram of the Silhouette scores
quartz()
hist(perf_grid, breaks = 20, xlab = "Silhouette score", xlim = c(-1, 1), main = NULL)

max_score <- max(perf_grid)
min_score <- min(perf_grid)
max_abs <- max(abs(max_score), abs(min_score))

quartz(height = 6.5, width = 6.5)
image.plot(x = eps_grid, y = minPts_grid, z = perf_grid, xlab = "eps", ylab = "minPts",
      main = 'Silhouette score', col = hcl.colors(64, palette = 'Blue-Red'),
      breaks = c(seq(-max_abs, 0, length=33)[-33], seq(0, max_abs, length=33)))

# Retrieve best parameter values
max_score <- max(perf_grid)
argmax_score <- which(perf_grid == max_score, arr.ind = TRUE)
best_eps <- eps_grid[argmax_score[2]]
best_minPts <- minPts_grid[argmax_score[1]]
best_eps
best_minPts
max_score

# Run the dbscan
dbs <- dbscan(x, best_eps, best_minPts)
dbs

quartz()
plot(x, col = dbs$cluster + 1L, pch=19)

# Let's try now with eps = 0.09 and minPts = 15
dbs <- dbscan(x, eps = 0.09, minPts = 15)
dbs
# Recovered the original clusters!
quartz()
plot(x, col = dbs$cluster + 1L, pch=19)

# Example of dataset where DBSCAN succeed and k-means & hierarchical fail

data("moons")
plot(moons, pch=19)

minPts = 3 # Dimensionality + 1
quartz()
kNNdistplot(moons, k = 3)
abline(h = 0.3, col = "red", lty = 2)

eps <- 0.3
minPts <- 3

dbs <- dbscan(moons, eps, minPts)
dbs

quartz()
plot(moons, col = dbs$cluster + 1L, pch=19)

# Let's try the other algorithms that we saw:
hclust.s <- hclust(dist(moons), method='single')
hclust.a <- hclust(dist(moons), method='average')
hclust.c <- hclust(dist(moons), method='complete')


# plot of the dendrograms
quartz()
par(mfrow=c(1,3))
plot(hclust.s, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(hclust.c, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(hclust.a, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

dev.off()

# Let's look at:
# - single linkage with three clusters
# - complete linkage with two clusters
# - average linkage with two clusters

cluster.c <- cutree(hclust.c, k=2) # euclidean-complete:
cluster.s <- cutree(hclust.s, k=3) # euclidean-single
cluster.a <- cutree(hclust.a, k=2) # euclidean-average


quartz()
par(mfrow=c(1,3))
plot(moons, col=cluster.c + 1L, pch=19, main='Complete')
plot(moons, col=cluster.a + 1L, pch=19, main='Average')
plot(moons, col=cluster.s + 1L, pch=19, main='Single')

# k-means

# Choose k -> elbow method
b <- NULL
w <- NULL
for(k in 1:10){
  
  k.means <- kmeans(moons, k)
  w <- c(w, sum(k.means$wit))
  b <- c(b, k.means$bet)
  
}

quartz()
matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
lines(1:10, w/(w+b), type='b', lwd=2)
# -> 4 clusters

quartz()
k.means <- kmeans(moons, centers = 4, nstart = 25)
plot(moons, col=k.means$cluster + 1L, pch=19, main='k-means')

#_______________________________________________________________________________
##### Exercises

##### Problem 3 of July 1, 2009
#####-----------------------------
# The Veritatis Diagram is a mysterious work attributed to Galileo. Some
# semiologists believe that some pages of the book hide a coded message;
# they also believe that these pages are characterized by an abnormal
# numerosity of some letters of the alphabet. The veritatis.txt file
# lists, for 132 pages of the book, the absolute frequencies of the five
# vowels of the Latin alphabet.
# a) By using an agglomerative clustering algorithm (Manhattan distance,
#    average linkage), identify two clusters and report suspicious pages;
# b) assuming that, for each of the two clusters, the five absolute 
#    frequencies are (approximately) normally distributed with the same 
#    covariance matrix perform a test to prove the existence of a difference
#    in the mean value of the two distributions;
# c) using five Bonferroni intervals of global confidence 90%, comment the
#    results of the test at point (b).

vowels <- read.table('veritatis.txt', header=T)
head(vowels)
dim(vowels)

### question a)
quartz()
plot(vowels)

HC <- hclust(dist(vowels, method='manhattan'), method = 'average')

quartz()
plot(HC, hang=-0.1, sub='', labels=F, xlab='')

# we cut the dendrogram at k=2 clusters
rect.hclust(HC, k=2)

pag <- cutree(HC, k=2)
table(pag)
which(pag==2)

quartz()
plot(vowels , col=pag+1, asp=1, pch=16, lwd=2)

### question b)
p  <- 5
n1 <- table(pag)[1]
n2 <- table(pag)[2]

# Verify gaussianity
mvn(vowels[pag=='1',])$multivariateNormality
mvn(vowels[pag=='2',])$multivariateNormality

# Test for independent Gaussian populations
t1.mean <- sapply(vowels[pag=='1',],mean)
t2.mean <- sapply(vowels[pag=='2',],mean)
t1.cov  <-  cov(vowels[pag=='1',])
t2.cov  <-  cov(vowels[pag=='2',])
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)

# Test: H0: mu.1-mu.2==0 vs H1: mu.1-mu.2!=0
delta.0 <- c(0,0,0,0,0)
Spinv   <- solve(Sp)
T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)
P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P

### question c)
alpha <- 0.1
IC <- cbind(t2.mean-t1.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(p*2), n1+n2-2),
            t2.mean-t1.mean,
            t2.mean-t1.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(p*2), n1+n2-2))
IC

graphics.off()

#_______________________________________________________________________________
##### Problem 2 of February 18, 2009
#####--------------------------------
# Friday, October 17, 2008, in Black Fortune skies a crash occurred between
# two artificial satellites. About a hundred debris were found on the
# ground (satellite.txt file). However, due to friction with the atmosphere it
# was not possible to define the origin of any debris. To clarify the 
# accident, the US Air Force requests to estimate the relative position of 
# the two points of impact on the ground.
# Assuming that the debris coming from the same satellite are scattered on
# the ground according to a normal law of unknown mean and covariance matrix:
# a) attribute the debris to the two satellites via a hierarchical clustering
#    algorithm that uses the Euclidean distance and the average linkage;
# b) report the numerosity of the two clusters and the value of the cophenetic
#    coefficient;
# c) introducing the appropriate hypotheses, identifying the points of impact
#    on the ground with the mean of the two normal distributions and assuming 
#    correct the allocation of the debris to the two satellites, provide an elliptical
#    confidence region (global level 99%) for the relative position of the two
#    points of impact on the ground.

satellite <- read.table('satellite.txt', header=T)
head(satellite)

quartz()
plot(satellite, asp = 1, pch=16)

dev.off()

### question a)

D.s <- dist(satellite)

HCa <- hclust(D.s, method = 'average')

quartz()
plot(HCa, hang=-0.1, sub='', xlab='', labels=F)

# we know that there are 2 clusters, so we cut the dendrograms 
# accordingly:
rect.hclust(HCa, k=2)

sata <- cutree(HCa, k=2)

quartz()
plot(satellite , col=sata+1, asp=1, pch=16, main='Average')

dev.off()

### question b)

table(sata)

coph.a <- cophenetic(HCa)

coph.sat <- cor(D.s, coph.a)
coph.sat


### question c)

p  <- 2
n1 <- table(sata)[1]
n2 <- table(sata)[2]

# Assumptions:
# - Normality
# - Independent populations
# - Homogeneity of covariance structures

# Verificy Gaussian assumption
mvn(satellite[sata=='1',])$multivariateNormality
mvn(satellite[sata=='2',])$multivariateNormality

t1.mean <- sapply(satellite[sata=='1',],mean)
t2.mean <- sapply(satellite[sata=='2',],mean)
t1.cov  <-  cov(satellite[sata=='1',])
t2.cov  <-  cov(satellite[sata=='2',])

# Homogeneity of covariances
t1.cov
t2.cov

Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)

# Elliptic confidence region at 99%
alpha <- 0.01
cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha, p, n1+n2-1-p)

# Characterize the ellipse:
# Directions of the axes
eigen(Sp)$vector

# Radius
r <- sqrt(cfr.fisher)

# Length of the semi-axes
r*sqrt(eigen(Sp)$values*(1/n1+1/n2))

quartz(width=14, height = 7)
par(mfrow=c(1,2))
plot(satellite , col=sata+1, asp=1, pch=16, main='Original data and groups')
plot(satellite, xlim=c(-23,-15), ylim=c(-25,-17), pch='', asp=1, 
     main='Elliptic region for the mean diff. (red - green)')
# confidence region and sample mean in blue
ellipse(center=t1.mean-t2.mean, shape=Sp*(1/n1+1/n2), radius=sqrt(cfr.fisher), 
        lwd=2, col='blue')

graphics.off()

#_______________________________________________________________________________
##### Problem 2 of February 28, 2013
#####--------------------------------
# In view of the week of haute couture in Paris, the fashion house La Boutin
# has decided to create a unique-edition pair of shoes decorated with 
# gemstones instead of crystals. During the processing of the pair of shoes,
# a bumbling craftsman caused the fall of the gemstones and of a
# box of crystals. The file preziosi.txt collects data related to the Cartesian
# coordinates of the found stones.
# a) Through a hierarchical clustering algorithm (Euclidean distance and
#    Ward linkage), identify the two clusters of stones.
# b) reporting the numerosity of the two groups identified at point a) and
#    the value of the cophenetic coefficient.
# c) Identify the gemstones with the smaller group. Having introduced
#    (and verified) the appropriate assumptions, estimate an elliptical region
#    containing 99% of gemstones (report the center, the length and the direction
#    of the principal axes of the ellipse).

stones <- read.table('preziosi.txt',header=TRUE)
plot(stones)

# question a)
p.e <- dist(stones, method='euclidean')
p.ew <- hclust(p.e, method='ward.D2')

plot(p.ew, hang=-0.1, sub='', xlab='', labels=F)
cl.ew <- cutree(p.ew, k=2) # euclidean-ward

# question b)
table(cl.ew)

coph.ew <- cophenetic(p.ew)
ew <- cor(p.e, coph.ew)
ew

# question c)
p.pr<-stones[which(cl.ew==2),]
n <- dim(p.pr)[1]
p <-dim(p.pr)[2]

plot(stones)
points(p.pr, pch=19)

# Verify gaussian assumptions
mvn(p.pr)$multivariateNormality

M <- sapply(p.pr, mean)
S <- cov(p.pr)
alpha <- 0.01
cfr.chisq <- qchisq(1-alpha,p)

# Characterize the ellipse:
# Axes directions:
eigen(S)$vectors
# Center:
M
# Radius of the ellipse:
r <- sqrt(cfr.chisq)
# Length of the semi-axes:
r*sqrt(eigen(S)$values)  

quartz()
plot(p.pr, asp = 1, col='gold', pch=19, xlim=c(-10,50))
points(M[1], M[2], pch = 4, cex = 1.5, lwd = 2)

ellipse(center=M, shape=S, radius=sqrt(cfr.chisq), col = 'black', lty = 2, center.pch = 4)
points(stones[which(cl.ew==1),])

graphics.off()

#_______________________________________________________________________________
##### Problem 3 of February 28, 2007
#####--------------------------------

# The dataset Pb3.txt reports Length (cm) Width (cm) of the chest of
# 50 sparrows. The biologist who has collected the measurements aims to 
# demonstrate that the sparrows can be divided into two distinct groups
# in terms of length and width of the chest. Help him to prove his theory 
# by implementing and commenting on the following analyses:
# (a) By means of an agglomerative hierarchical clustering algorithm that
#     uses the Manhattan distance and the Single linkage state if it is 
#     reasonable to cluster the data into two groups.
# (b) Implement a test to prove the difference of the means of the two groups.
# (c) Identify and comment the four Bonferroni intervals with global confidence
#     90% (lower bound, central value, upper bound) for:
#     - The difference of the mean of the variable length.
#     - The difference of the mean of the variable width.
#     - The difference of the mean of the sum of the variables length and width.
#     - The difference of the mean of the difference of variable length and
#       width.
# (d) verify the assumptions necessary to the implementation of the test.

sparrows <- read.table('Pb3.txt')
head(sparrows)
dim(sparrows)

### question (a)

quartz()
plot(sparrows, pch=16)

dev.off()

gruppi <- hclust(dist(sparrows, method='manhattan'), method='single')
plot(gruppi, hang=-0.1, sub='', xlab='', labels=F)

# cut in two groups
cluster <- cutree(gruppi, k=2)

plot(sparrows, pch=16, col = as.vector(cluster)+1)


### question (b)
g1 <- sparrows[cluster==1,]
g2 <- sparrows[cluster==2,]
g1
g2

# Test: H0: mu.1-mu.2==0 vs H1: mu.1-mu.2!=0
p  <- 2
n1 <- dim(g1)[1]
n2 <- dim(g2)[1]
alpha <- 0.10

mean1 <- sapply(g1,mean)
mean2 <- sapply(g2,mean)
cov1  <-  cov(g1)
cov2  <-  cov(g2)
Sp      <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)

delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (mean1-mean2-delta.0) %*% Spinv %*% (mean1-mean2-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)

pvalue <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
pvalue

### question (c)

dm <- (mean1-mean2)
A  <- rbind(c(1,0), c(0,1), c(1,1), c(1,-1))
k  <- dim(A)[1]

A.s2 <- diag(A%*%Sp%*%t(A))
A.dm <- A%*%(mean1-mean2)

Bonf <- cbind(inf=A.dm - qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ), 
              center=A.dm, 
              sup=A.dm + qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ))
Bonf

### question (d)
mvn(g1)$multivariateNormality
mvn(g2)$multivariateNormality

cov1
cov2


#_______________________________________________________________________________
##### Multidimensional Scaling (principal coordinate analysis)
#####---------------------------------------------------------

### Given the distances (dissimilarities) among n statistical units, look for 
### the k-dimensional representation (k small) of the n statistical units
### such that the distances (dissimilarities) among the representations
### of the n units are as close as possible to the original distances
### (dissimilarities) among the n units.


### Example 1: European cities
###---------------------------

help(eurodist)
eurodist
# road distances among European cities

# R function for multidimensional scaling: cmdscale
help(cmdscale)

location <- cmdscale(eurodist, k=2)
location

# I have to set asp=1 (equal scales on the two axes)
# to correctly represent Euclidean distances
quartz()
plot(location[,1], location[,2], type='n', asp=1, axes=FALSE, main="MDS of European cities",xlab='',ylab='')
text(location[,1], location[,2], labels=colnames(as.matrix(eurodist)), cex = 0.75, pos = 3)

# change the sign to get the North in the upper part of the plot
quartz()
plot(location[,1], -location[,2], type='n', asp=1, axes=FALSE, main="MDS of European cities",xlab='',ylab='')
text(location[,1], -location[,2], labels=colnames(as.matrix(eurodist)), cex = 0.75, pos = 3)

# compare the original matrix d_ij = d(x_i,x_j) and delta_ij = d(y_i,y_j) 
quartz()
plot(eurodist, dist(location))

# visualize the most different distances
quartz()
par(cex = 0.75, mar = c(10,10,2,2))
image(1:21, 1:21, asp=1, abs(as.matrix(dist(location)) - as.matrix(eurodist)), axes = F, xlab = '', ylab ='')
axis(1, at = 1:21, labels = colnames(as.matrix(eurodist)), las = 2, cex = 0.75)
axis(2, at = 1:21, labels = colnames(as.matrix(eurodist)), las = 1, cex = 0.75)
box()

# Rome-Athens
as.matrix(eurodist)[19,1]
as.matrix(dist(location))[19,1]

# Cologne-Geneve
as.matrix(eurodist)[6,8]
as.matrix(dist(location))[6,8]


# Compute the "stress": the higher it is, the worse
# the matching between original distances and their
# geometrical representation through MDS
Stressk <- NULL
for(k in 1:4)
{
  location.k <- cmdscale(eurodist, k)
  Stress <- (sum( (as.vector(eurodist) - as.vector(dist(location.k)))^2)  /
               sum( as.vector(location.k)^2))^(1/2)
  Stressk <- c(Stressk, Stress) 
}

quartz()
plot(1:4,Stressk,xlab='k',ylab='Stress',lwd=2)
# the stress increases for k>2 because of numerical problems 
# (the representation of minimal stress is not always found)

graphics.off()
