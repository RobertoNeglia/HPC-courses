## Ex4

library(sp)           ## Data management
library(lattice)      ## Data management
# library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics


dataset = read.table("walesharks.txt", header = T)

resp = dataset$sights  # resp (è il target): variabile per cui vogliamo studiare la dipendenza spaziale 
resp = log(resp)

dataset$log.sights = resp

reg = dataset$log.chlorofill   # regressore (se ho un modello non stationary)

coordinates(dataset) <- c('x','y')   # 'x', 'y' devono essere i nomi delle colonne

v <- variogram( resp ~ reg, dataset)


x11()
plot(v, main = 'Sample Variogram',pch=19)

v.fit = fit.variogram(v, vgm(0.5 , "Exp", 110000))   
v.fit

x11()
plot(v, v.fit, pch = 19)  


# parameters
# Estimates
g.tr <- gstat(formula = log.sights ~ log.chlorofill, data = dataset , model = v.fit)

existent_loc = c(dataset$x[1], dataset$y[1]) 


# a0
s0.new = data.frame(x= existent_loc[1] , y= existent_loc[2], log.chlorofill=0)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE)  
# 2.745206

# a1
s0.new = data.frame(x= existent_loc[1] , y= existent_loc[2], log.chlorofill=1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 

# 3.48432 - 2.745206 = 0.739114


## model for the log of the chlorofill
resp_c = dataset$log.chlorofill

v_c <- variogram( resp_c  ~ 1, dataset) 

x11()
plot(v_c, main = 'Sample Variogram',pch=19)


v.fit_c = fit.variogram(v_c, vgm(3.5 , "Sph", 50000))   
# come range metti un valore intorno al quale il variogram si assesta sull'asintoto
v.fit_c


# Plot of the fit
x11()
plot(v_c, v.fit_c, pch = 19)  


# prediction 253844.8, 385997.7).
s0.new = as.data.frame(matrix(c(253844.8,385997.7),1,2))

#                ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y')

#                      ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) <- c('x','y')

# Create a gstat object gstat(g.obj, id, formula, data, model, set,...)
#                     = nome_colonna_resp ~ 1,
g.tr <- gstat(formula = log.chlorofill    ~ 1, data = dataset, model = v.fit_c)
# model is the variogram model we use for the data

# Prediction of resp(s_0):
predict(g.tr, s0.new, BLUE = FALSE)  # by default predict uses Kriging




reg.s0 = 3.983721     # valore del regressore nella nuova location

#                               coordinata x , coordinata y           1,3 è la dimensione si s0.new
s0.new = as.data.frame(matrix(c( 253844.8,385997.7, reg.s0),1,3))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','log.chlorofill')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = log.sights      ~ log.chlorofill   , data = dataset, model = v.fit)

# prediction
predict(g.t, s0.new, BLUE = FALSE)

# 5.739903






















