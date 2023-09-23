### EX4



library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics


v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}

dataset = read.table("hotels.txt", header = T)

resp = dataset$price  # resp (è il target): variabile per cui vogliamo studiare la dipendenza spaziale 

coordinates(dataset) <- c('x','y')   # 'x', 'y' devono essere i nomi delle colonne

# a) under stationairty
# Hyp: second-order stationarity and isotropy
v <- variogram( resp ~ 1, dataset) 


x11()
plot(v, main = 'Sample Variogram',pch=19)

#                        vgm(sill, model, range, nugget)
v.fit = fit.variogram(v, vgm(4000 , "Sph", 500 , 1000 )) 
v.fit
x11()
plot(v, v.fit, pch = 19)  

# parameter a0
g.tr <- gstat(formula = price  ~ 1, data = dataset, model = v.fit)
predict(g.tr, dataset[1,], BLUE = TRUE)


# b)
n = dim(dataset)[1]
dummy_no = rep(0, n)
dummy_no[which(dataset$winter == 'no')] = 1 
dataset$dummy_no = dummy_no

v2 <- variogram( resp ~ dummy_no + distance + dummy_no:distance, dataset)

x11()
plot(v2, main = 'Sample Variogram',pch=19)

v.fit2 = fit.variogram(v2, vgm(1000 , "Sph", 1000 , 100 )) 
v.fit2
x11()
plot(v2, v.fit2 , pch = 19) 

# Parameters
g.tr <- gstat(formula = resp ~ dummy_no + distance + dummy_no:distance, data = dataset , model = v.fit)

### ------------
dum_1 = dataset[1,]  # dummy = 1
dum_0 = dataset[3,]  # dummy = 0

# Estimate a0_g
dum_1$distance = 0
obs_1 = dum_1

a0_1 = predict(g.tr, obs_1, BLUE = TRUE)
# a0_1 = 225.7315

dum_0$distance = 0
obs_0 = dum_0

a0_0 = predict(g.tr, obs_0, BLUE = TRUE)
# a0_0 = 444.4611

# Estimate a1_g
dumreg_1 = dataset[1,]
dumreg_0 = dataset[3,]

dumreg_1$distance = 1
y_1 = dumreg_1

predict(g.tr, y_1, BLUE = TRUE) 
a1_1 = 225.7196 - 225.7315
a1_1   # -0.0119

dumreg_0$distance = 1
y_0 = dumreg_0

predict(g.tr, y_0, BLUE = TRUE) 
a1_0 = 444.3928 - 444.4611
a1_0    # -0.0683


## -----------



# Estimates
#a0_2 = a0_noWin -> dummy = 1
s0.new = data.frame(x= 344076.8, y= 5072878, dummy_no = 1, distance=0)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 
# 225.7315

#a0_1 = a0_Win -> dummy = 0
s0.new = data.frame(x= 344076.8, y= 5072878, dummy_no = 0, distance=0)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 
# 444.4611

# a1_2 = a1_noWint -> dummy = 1
s0.new = data.frame(x= 344076.8, y= 5072878, dummy_no = 1, distance=1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 

a1_2 = 225.7196 - 225.7315 
# -0.0119


# a1_1 = a1_Wint -> dummy = 0
s0.new = data.frame(x= 344076.8, y= 5072878, dummy_no = 0, distance=1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 

a1_1 = 444.3928 - 444.4611 
# -0.0683






# c) 

v = v                   # variogram under stationarity 
v.gls = v2               # variogram without stationarity

v.fit = v.fit              # fit with stationarity 
v.gls.fit = v.fit2          # fit without stationarity 


# If you have summed more than 2 models, then change v.fit[ , ]
x11()
plot(v$dist, v$gamma, xlab='distance' ,ylab='semivariance', pch=19, col='skyblue1', ylim=c(0,8000))
points(v.gls$dist, v.gls$gamma, xlab='distance', ylab='semivariance', pch=19, col='red')

# they are different -> the drift is important
# we choose the model under non-stationary

# d) 
dummy.s0 = 0
d = sqrt((342362.58-342399.74)^2 + (5072518.24- 5072272.75)^2)
reg.s0 =  d    # valore del regressore nella nuova location

#                               coordinata x , coordinata y           1,3 è la dimensione si s0.new
s0.new = as.data.frame(matrix(c(342399.74, 5072272.75 , reg.s0, dummy.s0),1,4))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','distance', 'dummy_no')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = price ~  distance + dummy_no + dummy_no:distance  , data = dataset, model = v.fit2)

# prediction
predict(g.t, s0.new, BLUE = FALSE)










