### Ex4


library(fda)
library(KernSmooth)


dataset = read.table("spectra.txt",header=T)

dataset = t(dataset)

spectrum_1 = dataset[ , 1]

Xobs0 <- spectrum_1         # valore
abscissa <- c(1:80)         # ascissa
N <- length(abscissa) # number of locations of observations
n_functions = dim(dataset)[2]

x11()
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
plot(abscissa,Xobs0,xlab="t",ylab="observed data", type = "l")


# Set parameters
m <- 4        # spline order
degree <- m-1    # spline degree = order -1



# We can choose the correct number of basis via GCV
nbasis_range <- 6:30
gcv <- numeric(length(nbasis_range))
for (i in 1:length(nbasis_range)){
  basis_temp <- create.bspline.basis(range(abscissa), nbasis_range[i], m)
  gcv[i] <- smooth.basis(abscissa, Xobs0, basis_temp)$gcv
}
x11()
par(mfrow=c(1,1))
plot(nbasis_range,gcv)
n_basis_opt = nbasis_range[which.min(gcv)]
abline(v = n_basis_opt, col = 2)

n_basis_opt
nbasis = 11

# Create the basis
basis <- create.bspline.basis(rangeval = range(abscissa), nbasis=nbasis, norder=m)
# If breaks are not provided, equally spaced knots are created





# Evaluate the basis on the grid of abscissa
basis_eval <- eval.basis(abscissa, basis)   
dim(basis_eval) # number of data x number of basis
head(basis_eval)
# each column is the evaluation of the basis function in the location where we have the data
# we see that basis 6,7,8,9 are 0 in the first nodes and this is because of the support

# Fit via LS
est_coef = lsfit(basis_eval, Xobs0, intercept=FALSE)$coef
est_coef
# remove the intercept since it is already considered 








basis$params


DATA =dataset


n_functions = dim(DATA)[2]



data.fd <- Data2fd(y = DATA ,argvals = abscissa ,basisobj = basis)
x11()
plot.fd(data.fd)


# c)

library(fdakma)

# y0 = t(DATA)

abscissa2 = t(abscissa)

n_cluster = 3

fdakma_align <- kma(
  x=abscissa2, y0= y0 , n.clust = n_cluster, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson', 
  center.method = 'k-means',
  # seeds = c(1,21) # you can give a little help to the algorithm...
)

x11()
kma.show.results(fdakma_align)



### ----------------

# Smoothing per n funzioni
#x11()
z.all.rs0 = matrix(nrow = N, ncol = n_functions)
for(i in 1:n_functions){
  z.obs0.i = dataset[, i]  # i-th functional datum 
  # Create the basis
  basis.i = create.bspline.basis(rangeval=c(min(abscissa), max(abscissa)), nbasis=nbasis, norder=m)
  z.all.rs0[, i] = eval.fd(abscissa, (smooth.basis(argvals=abscissa, y=z.obs0.i, fdParobj=basis.i))$fd) 
  # Add the plot of the i_th smoothed curve 
  
  # # forse togli questo if 
  # if(i != 1)
  #   points(abscissa, eval.fd(abscissa, (smooth.basis(argvals=abscissa, y=z.obs0.i, fdParobj=basis.i))$fd), 
  #          type="l", col=rainbow(n)[i], lwd=2) 
}
# ho popolato z.all.rs0 e fatto il plot


# Plot of all the smoothed curves 
x11()
plot(abscissa, z.all.rs0[, 1], type="l", col=rainbow(n_functions)[1], lwd=2) 
for(i in 2:n_functions){
  points(abscissa, z.all.rs0[, i], type="l", col=rainbow(n_functions)[i], lwd=2) 
}



## 
# From the smoothing 
z0 = t(z.all.rs0)     # all the smoothed curves  (matrix N x n_functions)

num.clusts = 3      # number of clusters 
warp.method = 'affine'   # warping method (e.g.: 'NOalignment', 'shift', 'dilation', 'affine')
simil.method = 'd0.pearson'  # similarity method (e.g.: 'd0.pearson', 'd1.pearson') 
## QUANDO DICE DI USARE LA CORRELAZIONE COME SIMILARITY INTENDE d0.pearson
library(fdakma)
set.seed(3)
fda_kma = kma(
  x=abscissa, y0=z0, n.clust = num.clusts, 
  warping.method = warp.method, 
  similarity.method = simil.method,   
  center.method = 'k-means'
) 
# Output the results 
kma.show.results(fda_kma)
# GLI VANNO BENE QUESTI PLOTS ?? 


# Labels 
fda_kma$labels
















