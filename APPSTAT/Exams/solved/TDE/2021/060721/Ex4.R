### Ex4

library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics

dataset = read.table("colours.txt", header = T)

resp = dataset$revenue  # resp (è il target): variabile per cui vogliamo studiare la dipendenza spaziale 

coordinates(dataset) <- c('x','y')

# a) 
# Stationarity
v <- variogram( resp ~ 1, dataset) 

x11()
plot(v, main = 'Sample Variogram',pch=19)

#                        vgm(sill, model, range, nugget)
v.fit = fit.variogram(v, vgm(45 , "Sph",500))   

# Plot of the fit
x11()
plot(v, v.fit, pch = 19)  



g.tr <- gstat(formula = revenue ~ 1, data = dataset, model = v.fit)
dataset[1,]
predict(g.tr, dataset[1,], BLUE = TRUE)
# i obtain the exact value and variance 0 

# b)
# create the dummy 
n = 84
dummy_red = rep(0,n)
dummy_orange = rep(0,n)
dummy_red[which(dataset$colour == 'red')] = 1
dummy_orange[which(dataset$colour == 'orange')] = 1

v2=variogram( resp ~ dummy_red + dummy_orange , data = dataset)

x11()
plot(v2,pch=19)


#                        vgm(sill, model, range, nugget)
v.fit2 = fit.variogram(v2, vgm(5 , "Sph",2000))   

# Plot of the fit
x11()
plot(v2, v.fit2, pch = 19)
v.fit2

dataset$dummy_red = dummy_red
dataset$dummy_orange = dummy_orange

g.tr <- gstat(formula = revenue ~  dummy_red + dummy_orange, data = dataset, model = v.fit2)

predict(g.tr, dataset[1 ,], BLUE = TRUE)
predict(g.tr, dataset[2 ,], BLUE = TRUE)
predict(g.tr, dataset[3 ,], BLUE = TRUE)

# c) 
# Choose the model
# Let's compare the variogram of the data and of the residuals:

                  # variogram under stationarity 
v.gls = v2             # variogram without stationarity

              # fit with stationarity 
v.gls.fit = v.fit2          # fit without stationarity 


x11()
plot(v$dist, v$gamma, xlab='distance' ,ylab='semivariance', pch=19, col='skyblue1', ylim=c(0,60))
points(v.gls$dist, v.gls$gamma, xlab='distance', ylab='semivariance', pch=19, col='red')

# drift significant -> choose the non-stationary model


##  Prediction at a single new location

#                               coordinata x , coordinata y           1,3 è la dimensione si s0.new
s0.new = as.data.frame(matrix(c(514811.55 ,  5037308.54  , 1 , 0),1,4))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','dummy_red', 'dummy_orange')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = revenue    ~  dummy_red + dummy_orange  , data = dataset, model = v.fit2)

# prediction
predict(g.t, s0.new, BLUE = FALSE)



s0.new = as.data.frame(matrix(c(514811.55 ,  5037308.54  , 0 , 1),1,4))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','dummy_red', 'dummy_orange')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = revenue    ~  dummy_red + dummy_orange  , data = dataset, model = v.fit2)

# prediction
predict(g.t, s0.new, BLUE = FALSE)




s0.new = as.data.frame(matrix(c(514811.55 ,  5037308.54  , 0 , 0),1,4))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','dummy_red', 'dummy_orange')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = revenue    ~  dummy_red + dummy_orange  , data = dataset, model = v.fit2)

# prediction
predict(g.t, s0.new, BLUE = FALSE)



s0.new = as.data.frame(matrix(c(514811.55 ,  5037308.54  , 0 , 0),1,4))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','dummy_red', 'dummy_orange')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = revenue    ~  dummy_red + dummy_orange  , data = dataset, model = v.fit2)

# prediction
predict(g.t, s0.new, BLUE = FALSE)














