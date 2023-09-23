### Ex1

dataset <- read.table('holiday.txt', header = T)  # se serve aggiungi il nome delle colonne col.names=c('',..)
dataset


group1   <- factor(dataset$location)   # qui species è il primo grouping
group2   <- factor(dataset$type)    # qui gender è il secondo grouping

x  <- dataset[,1]   # salvo la variabile numerica in x 

g <- length(levels(group1))
b <- length(levels(group2))
n <- length(x)/(g*b)   



group_12 <- group1
levels(group_12) <- c('C_H','C_B','C_A','R_H','R_B','R_A')  #variable which divides the g*b groups (only for plot)
group_12[group1=='Chiavari' & group2=='hotel'] <- 'C_H'
group_12[group1=='Chiavari' & group2=='bb'] <- 'C_B'
group_12[group1=='Chiavari' & group2=='apt'] <- 'C_A'
group_12[group1=='Rapallo' & group2=='hotel'] <- 'R_H'
group_12[group1=='Rapallo' & group2=='bb'] <- 'R_B'
group_12[group1=='Rapallo' & group2=='apt'] <- 'R_A'


load("mcshapiro.test.RData") 
Ps <- c(shapiro.test(x[ group_12==levels(group_12)[1]])$p,
        shapiro.test(x[ group_12==levels(group_12)[2]])$p,
        shapiro.test(x[ group_12==levels(group_12)[3]])$p,
        shapiro.test(x[ group_12==levels(group_12)[4]])$p,
        shapiro.test(x[ group_12==levels(group_12)[5]])$p,
        shapiro.test(x[ group_12==levels(group_12)[6]])$p)   
Ps


bartlett.test(x, group_12)

fit.aov2.int <- aov(x ~ group1 + group2 + group1:group2)  # aov(y ~ factor1 + factor2 + interaction)
summary.aov(fit.aov2.int)

# b) interaction is not significant
fit.aov2.ad <- aov(x ~ group1 + group2)
summary.aov(fit.aov2.ad)


fit<- aov(x ~ group2)   # qui sto rimuovendo group1
summary.aov(fit)

# c)
# variance:
DF <- fit$df # fit è il modello finale

Spooled <- sum(fit$res^2)/DF
Spooled  # estimate for the variance 


# Estimate the great mean mu:
m <- mean(x)

# Point-wise estimates of the means in each group
m1  <- mean(x[group2=='hotel'])   
m2 <- mean(x[group2=='bb'])  
m3 <- mean(x[group2=='apt'])  

# Estimate tau.i:
tau1  <- mean(x[group2=='hotel'])  - m  
tau2 <- mean(x[group2=='bb'])  - m  # devo avere g righe
tau3 <- mean(x[group2=='apt'])  - m 

tau = c(beta1 = tau1, beta2 = tau2, beta = tau3)

# d)
group = group2
g = b
n       <- length(group)

ng      <- table(group)       # number of obs. in each group
treat   <- levels(group)      # levels of the treatment

k <- g*(g-1)/2
alpha= 0.05

Mediag  <- tapply(x, group, mean)
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)

# BonfCI for all the DIFFERENCES
ICrange=NULL
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    print(paste(treat[i],"-",treat[j]))        
    print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
    ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }}
















