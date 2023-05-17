setwd("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 7 - 12042021")
load("D:/RTDA/Didattica/Applied Statistics MATE 20-21/Lab 5 - 23042021/mcshapiro.test.RData")

#_______________________________________________________________________________
##### Problem 1 of 09/07/08
##### (questions b and c)
#####------------------------

# The client.txt dataset contains data on 150 customers of
# PoliBank. For each customer we are given age [years], money invested
# at low risk [thousands of euros] (safemoney variable) and money invested
# at high risk [thousands of euros] (riskymoney variable).
# [a) Using only the variable age, cluster customers in three groups and
#     describe them in terms of age. Use a hierarchical agglomerative 
#     algorithm based on Euclidean distance and single linkage. Report
#     cophenetic coefficient and the size of the clusters.
#  b) Introducing the appropriate assumptions about distributions of the 
#     variables safemoney and riskymoney within the three groups, perform
#     a MANOVA to see if there is statistical evidence of a difference 
#     in the joint distributions of safemoney and riskymoney variables 
#     in the three groups.
#  c) Comment the result of MANOVA by means of suitable Bonferroni 
#     intervals with global confidence 90%.

client <- read.table('client_class.txt', header=T)  # Clustering already performed at 
                                                    # point a), labels inserted as 1st 
                                                    # column
head(client)
dim(client)

# prepare data
var.risp <- client[,2:3]
group.names <- client[,1]
dim(var.risp)[2]
levels(group.names)

# variables:
p <- 2
g <- 3
i1 <- which(group.names=='young')
i2 <- which(group.names=='adult')
i3 <- which(group.names=='old')
ng <- c(length(i1),length(i2),length(i3)) 
ng

N <- sum(ng)

### question b)

### One-way MANOVA:
### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), [p=2]
###       X.ij, mu, tau.i in R^2, i=1,2,3 

# verify assumptions
# 1) normality
Ps <- c(mcshapiro.test(var.risp[ i1, ])$p,
        mcshapiro.test(var.risp[ i2, ])$p,
        mcshapiro.test(var.risp[ i3, ])$p)
Ps

# 2) homogeneity in variance
S1 <-  cov(var.risp[ i1, ])
S2 <-  cov(var.risp[ i2, ])
S3 <-  cov(var.risp[ i3, ])

x11()
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
dev.off()

# Fit the model:
fit <- manova(as.matrix(var.risp) ~ group.names)
summary.manova(fit,test="Wilks")

# who's the responsible?
summary.aov(fit,test="Wilks")
# all the variables

# let's see if there is difference in the levels of the treatment -> question c)

### question c)
alpha <- 0.10
k <- p*g*(g-1)/2
k

qT <- qt(1-alpha/(2*k), N-g)

# I need the diagonal of W
fit$res   # residuals of the estimated model
W <- diag(t(fit$res) %*% fit$res)/(N-g)   
W
# mean within the groups
m1 <- colMeans(var.risp[i1,])
m2 <- colMeans(var.risp[i2,])
m3 <- colMeans(var.risp[i3,])
m1
m2
m3

Bf12 <- cbind(m1-m2 - qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[1]+1/ng[2])*W), m1-m2, m1-m2 + qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[1]+1/ng[2])*W))
Bf23 <- cbind(m2-m3 - qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[2]+1/ng[3])*W), m2-m3, m2-m3 + qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[2]+1/ng[3])*W))
Bf31 <- cbind(m3-m1 - qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[3]+1/ng[1])*W), m3-m1, m3-m1 + qt(1 -alpha/(2*k), N-g) * sqrt((1/ng[3]+1/ng[1])*W))

IC <- list(young_adult=Bf12, adult_old=Bf23, old_young=Bf31)
IC

#_______________________________________________________________________________
##### Problem 3 of 29/06/11
#####------------------------
# Juan de los Euros, a known driver of Santander, has measured the duration of
# his recent trips to / from the airport (file time.txt).
# a) Having fitted a two-factor ANOVA additive model (center-aero / aero-center,
#    Weekday / weekend) provide point estimated of the means and the variances  
#    of the four possible types of travel.
# b) On the basis of the model (a) perform a 90% test to test the significance
#    of the factor-center aero / aero-center.
# c) On the basis of the model (a) perform a 90% test to test the significance
#    of the factor weekday / weekend.
# d) Based on the tests (b) and (c) propose a possible reduced model and re-
#    estimate point-wise - coherently with the reduced model - the means and
#    variances of the four possible types of travel.

euros <- read.table('time.txt', header=T)
euros
attach(euros)

# question a)

### Two-ways ANOVA
### Model without interaction (additive model): 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
###     i=1,2 (effect direction centre-aero/aero-centre), 
###     j=1,2 (effect day weekday/weekend)

g <- 2
b <- 2
p <- 1
n <- 5
N <- n*g*b

# Verify the assumptions
# 1) normality (univariate) in each group
Ps <- c(shapiro.test(durata[ AR=='aero_centro' & FF=='festivo' ])$p,
        shapiro.test(durata[ AR=='aero_centro' & FF=='feriale' ])$p,
        shapiro.test(durata[ AR=='centro_aero' & FF=='festivo' ])$p,
        shapiro.test(durata[ AR=='centro_aero' & FF=='feriale' ])$p)
Ps

# 2) homogeneity of variances
bartlett.test(list(durata[ AR=='aero_centro' & FF=='festivo' ],
                   durata[ AR=='aero_centro' & FF=='feriale' ],
                   durata[ AR=='centro_aero' & FF=='festivo' ],
                   durata[ AR=='centro_aero' & FF=='feriale' ]))

# Fit the model:
fit <- aov(durata ~ AR + FF)
summary(fit)

names(fit)

# Estimate variances
W <- sum(fit$residuals^2)  # SS_res
var <- W/(g*b*n-g-b+1)     # SS_res/gdl(res)
var

# Estimate the great mean mu:
m <- mean(euros[,1])

# Estimate tau.i, beta.j:
tauAC  <- mean(euros[euros$AR=='aero_centro',1]) - m  # tau.1
tauCA  <- mean(euros[euros$AR=='centro_aero',1]) - m  # tau.2

betaFest <- mean(euros[euros$FF=='festivo',1]) - m  # beta.1
betaFer  <- mean(euros[euros$FF=='feriale',1]) - m  # beta.2

# Point-wise estimates of mean duration of travels
# (model without interaction!)
mAC_Fest <- m + tauAC + betaFest
mAC_Fer  <- m + tauAC + betaFer
mCA_Fest <- m + tauCA + betaFest
mCA_Fer  <- m + tauCA + betaFer

# questions b)/c)

summary(fit)

# question d)

### Reduced model (one-way ANOVA): 
### X.jk = mu + beta.j + eps.jk; eps.jk~N(0,sigma^2), 
###     j=1,2 (effect day weekday/weekend)

fit.red <- aov(durata ~ FF)
summary(fit.red)

W <- sum(fit.red$residuals^2)
var <- W/(g*b*n-b)
var

# point-wise estimate of the mean duration of the trips
# (one-way model!)
mAC_Fest <- m + betaFest
mAC_Fer  <- m + betaFer
mCA_Fest <- m + betaFest
mCA_Fer  <- m + betaFer

#_______________________________________________________________________________
##### Problem 1 of 12/02/08
#####------------------------
# In a medical research center, concentrations of interferon gamma-6 and 
# Interferon gamma-7 were measured in the blood of some patients who had
# been infected by Human Papilloma Virus (PV.txt file). The experiment
# included patients that have different medical profiles in terms
# Relapse and HPV-Clearance:
# Relapse: extinct infection (A) or ongoing infection (B);
# HPV-Clearance: extinct we aim to find out if Relapse and Clearance factors
# have an effect on interferon gamma-6 and interferon gamma-7 distributions.
# a) Introduce an appropriate statistical model that justifies the use of
#    MANOVA as a tool for the analysis of these data.
# b) Perform a test to verify the interaction between the factors.
# c) Identify the factor or factors that generate effects statistically
#    significant.
# c) Construct Bonferroni confidence intervals with global coverage 90% 
#    that clarify the conclusions drawn at point (c).

PV <- read.table('PV.txt', header = T)
PV

N <- dim(PV)[1]
p <- 2
g <- b <- 2
n <- N/(g*b)
n

var.risp <- PV[,1:2]

### question a)

### Two-ways MANOVA
### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N_p(0,Sigma), [p=2]
###     i=1,2 (effect REL), j=1,2 (effect HPV),
###     X.ijk, mu, tau.i, beta.j, gamma.ij in R^2

### We don't verify the assumptions (very few data)

### question b)

man.int <- manova(as.matrix(var.risp) ~  HPV + REL + HPV * REL, data = PV)
summary(man.int, test = 'Wilks')

### question c)

### Model without interaction (additive model): 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N_p(0,Sigma), [p=2]
###     i=1,2 (effect REL), j=1,2 (effect HPV),
###     X.ijk, mu, tau.i, beta.j in R^2

man     <- manova(as.matrix(var.risp) ~ HPV + REL, data = PV)
summary(man, test = 'Wilks')

summary.aov(man) # interesting!
# factor relapse seems to have effect on interferon gamma-6,
# factor HPV (virus) seems to have effect on interferon gamma-7
# but we don't know in which sense!

# let's clarify with Bonferroni:
SSres <- t(man$residuals) %*% man$residuals / (N-g-b+1) # model without interaction

k <- g*(g-1)/2*p + b*(b-1)/2*p
k

qT <- qt(1 - 0.1/(2*k), N-g-b+1)

attach(PV)

m6_rel <- tapply(gamma6, REL, mean)
m7_rel <- tapply(gamma7, REL, mean)
m6_HPV <- tapply(gamma6, HPV, mean)
m7_HPV <- tapply(gamma7, HPV, mean)

m6_rel
m7_rel
m6_HPV
m7_HPV

REL6   <- c(diff(m6_rel) - qT * sqrt( SSres[1,1] * (1/6+1/6) ),
            diff(m6_rel) + qT * sqrt( SSres[1,1] * (1/6+1/6) ))
REL7   <- c(diff(m7_rel) - qT * sqrt( SSres[2,2] * (1/6+1/6) ),
            diff(m7_rel) + qT * sqrt( SSres[2,2] * (1/6+1/6) ))

HPV6   <- c(diff(m6_HPV) - qT * sqrt( SSres[1,1] * (1/6+1/6) ),
            diff(m6_HPV) + qT * sqrt( SSres[1,1] * (1/6+1/6) ))
HPV7   <- c(diff(m7_HPV) - qT * sqrt( SSres[2,2] * (1/6+1/6) ),
            diff(m7_HPV) + qT * sqrt( SSres[2,2] * (1/6+1/6) ))

detach(PV)

Bf <- list(B_A_6 = REL6, B_A_7 = REL7, HPVest_HPVpres_6 = HPV6, HPVest_HPVpres_7 = HPV7)
Bf

# interf gamma 6 higher when the infection is extint (REL=A)
# interf gamma 7 higher when the virus is extint (HPVest = +)

#_______________________________________________________________________________
##### Problem 3 of 28/02/13
#####------------------------
# For security reasons, the direction of the Paris Louvre museum has undertaken
# a campaign of control of the tourists flow in the museum. During the first phase
# of monitoring ("louvre.txt" file), the durations [minutes] of the visits were
# measured in the Museum for 360 tourists from Europe, USA and Japan, also recording
# the type of visit (guided, with audio guide or without guide).
# a) Formulate a suitable (complete) model for the duration of a visit of the museum
#    with respect to the two factors nationality and type of visit; in particular,
#    introduce and verify the assumptions of the model.
# b) Through a suitable test, discuss the possibility of removing the interaction 
#    term and possibly reduce the model.
# c) On the basis of the model at step b), test the effect of the factors nationality
#    and type of visit on the average time of the visit and, if appropriate, propose
#    a reduced model.
# d) Provide the security managers of the four intervals museum with Bonferroni confidence
#    (globally 90%) for the mean and variance of the visit to the museum time for 
#    homogeneous groups of identified visitors at step c).

museo <- read.table('louvre.txt',header=TRUE)
attach(museo)

# question a)
### Two-ways ANOVA
### Model with interaction (complete model): 
### X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,Sigma), 
###     i=1,2,3 (effect TYPE OF VISIT), j=1,2,3 (effect NAZIONALITY),
###     X.ijk, mu, tau.i, beta.j, gamma.ij in R

fit <- aov(tempo ~ tipo + nazione + nazione:tipo, data=museo)
summary(fit)

# Verify assumptions
t <- tipo:nazione
st <- NULL
for(i in 1:9)
  st<-c(st, 
        shapiro.test(tempo[which(t==levels(t)[i])])$p )
st

bartlett.test(tempo,tipo:nazione)

# question b)
### Two-ways ANOVA
### Model without interaction (additive model): 
### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,Sigma), 
###     i=1,2,3 (effect TYPE OF VISIT), j=1,2,3 (effect NAZIONALITY),
###     X.ijs, mu, tau.i, beta.j in R

fit2 <- aov(tempo ~ tipo + nazione, data=museo)
summary(fit2)

# question c)
### Reduced model: One-way ANOVA
### X.ik = mu + tau.i + eps.ik; eps.ik~N(0,Sigma), 
###     i=1,2,3 (effect TYPE OF VISIT)
###     X.ijs, mu, tau.i, beta.j in R

fit3 <- aov(tempo ~ tipo, data=museo)
summary(fit3)

# question e)
# g=3 Groups identified by the levels of the factor TYPE
# k=4: g*(g-1)/2 comparisons between the means + 1 conf int on the variance 
n <- dim(museo)[1]
g <- 3
k <- g*(g-1)/2+1
S <- sum(residuals(fit3)^2)/(n-g)

alpha<- .1

Mg  <- tapply(museo[,1], tipo, mean) 

label <- levels(factor(tipo))
n1 <- length(museo[tipo==label[1],1])
n2 <- length(museo[tipo==label[2],1])
n3 <- length(museo[tipo==label[3],1])
t <- qt(1-alpha/(2*k),n-g)

# Conf int for the means
ICB1<-data.frame(L=Mg[1]-sqrt(S*(1/n1))*t,C=Mg[1],U=Mg[1]+sqrt(S/n1)*t)
ICB2<-data.frame(L=Mg[2]-sqrt(S*(1/n2))*t,C=Mg[2],U=Mg[2]+sqrt(S/n2)*t)
ICB3<-data.frame(L=Mg[3]-sqrt(S*(1/n3))*t,C=Mg[3],U=Mg[3]+sqrt(S/n3)*t)
ICB<-data.frame(rbind(ICB1,ICB2,ICB3))
ICB

# Conf int for variances
chi_u <- qchisq(alpha/(2*k),n-g)
chi_l <- qchisq(1-alpha/(2*k),n-g)
ICBV <- data.frame(L=(n-g)*S/chi_l,C=S,U=(n-g)*S/chi_u)
ICBV

detach(museo)

#_______________________________________________________________________________
##### Pb 3 of 05/09/08
#####-------------------

# The West Sussex Bread Association has randomly selected 60 business trade
# in which doughnuts are commonly sold. 30 activities are based
# in the city of Brighton and 30 in the town of Worthing. For each of the 
# two cities, in 10 activity the price of a plain bagel was recorded, in 10
# the price of a doughnut filled with cream and in 10 the price of a doughnut
# filled with jam. The data are reported in doughnut.txt dataset.
# a) Describe the ANOVA model you deem appropriate for the analysis of these data.
# b) Identifying the factors that significantly influence the distribution
#    of the price of doughnuts, identify a possible reduced model.
# c) using Bonferroni's inequality estimate through bilateral confidence
#    intervals (with global confidence 95%) the means and variances of the
#    subpopulations associated with the reduced model identified at step (b).

ciambelle <- read.table('doughnut.txt', header=TRUE)
ciambelle

attach(ciambelle)

# question a)
# ANOVA two-ways
# Model with interaction (complete model): 
# X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk; eps.ijk~N(0,sigma^2), 
#     i=1,2 (effect city), j=1,2,3 (effect type)

fit.c <- aov(prezzo ~ citta + tipo + citta:tipo)
summary(fit.c)

p.val <- c(shapiro.test(prezzo[which(citta==levels(citta)[1] & tipo==levels(tipo)[1])])$p,
           shapiro.test(prezzo[which(citta==levels(citta)[1] & tipo==levels(tipo)[2])])$p,
           shapiro.test(prezzo[which(citta==levels(citta)[1] & tipo==levels(tipo)[3])])$p,
           shapiro.test(prezzo[which(citta==levels(citta)[2] & tipo==levels(tipo)[1])])$p,
           shapiro.test(prezzo[which(citta==levels(citta)[2] & tipo==levels(tipo)[2])])$p,
           shapiro.test(prezzo[which(citta==levels(citta)[2] & tipo==levels(tipo)[3])])$p)
p.val

bartlett.test(prezzo, citta:tipo)

# question b)
# Model without interaction (additive model): 
# X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
#     i=1,2 (effect city), j=1,2,3 (effect type)
fit.c2 <- aov(prezzo ~ citta + tipo)
summary(fit.c2)

# one-way ANOVA 
# X.jk = mu + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
#     j=1,2,3 (effect type)
fit.c3 <- aov(prezzo ~ tipo)
summary(fit.c3)

# question c)
N <- dim(ciambelle)[1]
g <- length(levels(tipo))
DF <- N-g

alpha <- .05
k <- g+1

qT <- qt(1-alpha/(2*k), DF)
qCinf <- qchisq(1 - alpha / (2*k), DF)
qCsup <- qchisq(alpha / (2*k), DF)

Spooled <- (t(fit.c3$res) %*% fit.c3$res)/DF   
Spooled

m1 <- mean(ciambelle[which(tipo==levels(tipo)[1]),1])
m2 <- mean(ciambelle[which(tipo==levels(tipo)[2]),1])
m3 <- mean(ciambelle[which(tipo==levels(tipo)[3]),1])
medie <- c(m1,m2,m3)

ng <- c(length(which(tipo==levels(tipo)[1])),length(which(tipo==levels(tipo)[2])),length(which(tipo==levels(tipo)[3])))

BF    <- rbind(cbind(inf=medie - sqrt(as.vector(Spooled) / ng) * qT,
                     sup=medie + sqrt(as.vector(Spooled) / ng) * qT),
               c(inf=Spooled * DF / qCinf,
                 sup=Spooled * DF / qCsup))
BF

detach(ciambelle)

