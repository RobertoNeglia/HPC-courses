### Ex4



library(sp)           ## Data management
library(lattice)      ## Data management
# library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics



dataset = read.table("revenues.txt", header = T)

resp = dataset$revenue  # resp (è il target): variabile per cui vogliamo studiare la dipendenza spaziale 

reg = dataset$population

## Define the sample coordinates
coordinates(dataset) <- c('x','y')   # 'x', 'y' devono essere i nomi delle colonne

v <- variogram( resp ~ reg, dataset)


x11()
plot(v, main = 'Sample Variogram',pch=19)


v.fit = fit.variogram(v, vgm( 500, "Sph", 1500))   
# come range metti un valore intorno al quale il variogram si assesta sull'asintoto
v.fit


# Plot of the fit
x11()
plot(v, v.fit, pch = 19) 

## 2. Caso con solo regresore numerico

#             formula = nome_colonna_resp ~ nome_colonna_reg
g.tr <- gstat(formula = revenue  ~ population , data = dataset , model = v.fit)

existent_loc = c(dataset$x[1], dataset$y[1]) 


# a0
#                                                            nome_colonna_reg = 0
s0.new = data.frame(x= existent_loc[1] , y= existent_loc[2], population = 0)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE)  
# -29.03012

# a1
#                                                           nome_colonna_reg = 1
s0.new = data.frame(x= existent_loc[1] , y= existent_loc[2], population = 1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) 
# a1 = ciò che predic qui - a0
# -29.0092 + 29.03012 = 0.02092




# b) FATTO CON LM
# 
# 
# fm <- lm(dataset$population ~ dataset$distance)
# summary(fm)
# 
# 
# d0 = sqrt((514711.6 - 514703.8)^2 + (5033903.0 - 5035569.3)^2)
# 
# p_s0 = 7555.30254-0.85395*d0 
# 
# 
# 
# reg.s0 = p_s0     # valore del regressore nella nuova location
# 
# #                               coordinata x , coordinata y           1,3 è la dimensione si s0.new
# s0.new = as.data.frame(matrix(c( 514703.8          , 5035569.3         , reg.s0),1,3))
# 
# #                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
# names(s0.new) = c('x','y','population')
# 
# #                 ('nome colonna coord x', 'nome colonna coord y')
# coordinates(s0.new) = c('x','y')
# 
# # create a g.stat object
# #           formula = nome_colonna_resp ~ nome_colonna_reg
# g.t <- gstat(formula = revenue   ~  population         , data = dataset, model = v.fit)
# 
# # prediction
# predict(g.t, s0.new, BLUE = FALSE)



# b) FATTO CON SPATIAL


resp_p = dataset$population  # resp (è il target): variabile per cui vogliamo studiare la dipendenza spaziale 

reg_p = dataset$distance



v_p <- variogram( resp_p ~ reg_p, dataset)


x11()
plot(v_p, main = 'Sample Variogram',pch=19)


v.fit_p = fit.variogram(v_p, vgm( 3500, "Sph", 2000, 1000))   
 
# come range metti un valore intorno al quale il variogram si assesta sull'asintoto
v.fit_p


# Plot of the fit
x11()
plot(v_p, v.fit_p, pch = 19) 

d0 = sqrt((514711.6 - 514703.8)^2 + (5033903.0 - 5035569.3)^2)

#                               coordinata x , coordinata y   1,2 è la dimensione si s0.new  -> DA CONTROLLARE
s0.new = as.data.frame(matrix(c( 514703.8          , 5035569.3         , d0),1,3))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','distance')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = population   ~  distance         , data = dataset, model = v.fit_p)

# prediction
predict(g.t, s0.new, BLUE = FALSE)

pop.s0 = 6091.035


#                               coordinata x , coordinata y           1,3 è la dimensione si s0.new
s0.new = as.data.frame(matrix(c( 514703.8          , 5035569.3         , pop.s0),1,3))

#                 ('nome colonna coord x', 'nome colonna coord y', 'nome colonna regressore')
names(s0.new) = c('x','y','population')

#                 ('nome colonna coord x', 'nome colonna coord y')
coordinates(s0.new) = c('x','y')

# create a g.stat object
#           formula = nome_colonna_resp ~ nome_colonna_reg
g.t <- gstat(formula = revenue   ~  population         , data = dataset, model = v.fit)

# prediction
predict(g.t, s0.new, BLUE = FALSE)

# 97.57383

