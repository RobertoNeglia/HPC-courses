
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics




data=read.table('hotels.txt')
attach(data)

dWin = ifelse(winter=='yes', 1, 0)

data <- data.frame(cbind(x, y, dWin, distance, price))
names(data) <- c('x','y','dWin','distance','price')
coordinates(data) <- c('x','y')
head(data)

#a)
v=variogram(price ~ 1, data=data)
plot(v,pch=19)

#fitting the variogram (with gls)
v.fit <- fit.variogram(v, vgm(1000, "Sph", 1500, 500))
plot(v, v.fit, pch = 3)
v.fit

g.tr <- gstat(formula = price ~ 1, data = data, model = v.fit)
predict(g.tr, data[1,], BLUE = TRUE)


#b)
v2=variogram(price ~ dWin + distance + distance:dWin, data=data) #??? -> 1 or something else
plot(v2,pch=19)

#fitting the variogram (with gls)
v2.fit <- fit.variogram(v2, vgm(1000, "Sph", 1500, 500))
plot(v2, v2.fit, pch = 3)
v2.fit

g.tr <- gstat(formula = price ~ dWin + distance + distance:dWin, data = data, model = v2.fit)

#a0_nowin
s0.new = data.frame(x=0, y=0, dWin = 0, distance=0)
coordinates(s0.new) = c('x','y')
predict(g.tr, data[1,], BLUE = TRUE) #a0_nowin = 204.2135

#a0_win + a0_nowin
s0.new = data.frame(x=0, y=0, dWin = 1, distance=0)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) #a0_win = 447.1602 - 204.2135 = 242.9467

#a1_nowin + a0_nowin
s0.new = data.frame(x=0, y=0, dWin = 0, distance=1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE) #a1_nowin = 220.7112 - 204.2135 = 16.4977

#a01 + a02 + a11 + a12
s0.new = data.frame(x=0, y=0, dWin = 1, distance=1)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE)
#a1_win = 447.0909  - 204.2135 - 242.9467 - 16.4977 = -16.567


#c)
s0_dist = sqrt((342399.74-342362.58)^2+(5072272.75-5072518.24)^2)

s0.new = data.frame(x=342399.74, y=5072272.75, dWin = 1, distance=s0_dist)
coordinates(s0.new) = c('x','y')
predict(g.tr, s0.new, BLUE = TRUE)
