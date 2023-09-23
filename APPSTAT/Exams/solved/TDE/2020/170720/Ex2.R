### Ex2

dataset = read.table("purity.txt", header = T)

pop1 = dataset[1:8, ]
pop2 = dataset[9:16, ]


n1 <- dim(pop1)[1] 
n2 <- dim(pop2)[1] 
p  <- dim(pop1)[2] 

pop1.mean <- sapply(pop1,mean)
pop2.mean <- sapply(pop2,mean)
pop1.cov  <-  cov(pop1)
pop2.cov  <-  cov(pop2)
Sp      <- ((n1-1)*pop1.cov + (n2-1)*pop2.cov)/(n1+n2-2)

# gaussianity
load("mcshapiro.test.RData") 
mcshapiro.test(pop1)$pvalue
mcshapiro.test(pop2)$pvalue


# H0: pop.mean1 - pop.mean2 >= 0  vs H1: pop.mean1 - pop.mean2 < 0
# rifiuto se z.i è troppo piccola -> z.i < - z(1-alpha)

z.i <- (pop1.mean-pop2.mean)/sqrt(diag(Sp)*(1/n1+1/n2))
p.i <- pnorm(z.i)  #p value


alpha = 0.1
p.BH <- p.adjust(p.i, method='BH')

which(p.BH < alpha)    # 2,3,5,7,9

# checkcon Bonf
alpha = 0.1
p.Bf <- p.adjust(p.i, method='bonferroni')
which(p.Bf < alpha)    # 2,3,5,9 

k = p
IC <- cbind(pop1.mean-pop2.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k), n1+n2-2),
            pop1.mean-pop2.mean,
            pop1.mean-pop2.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k), n1+n2-2))
IC
# visto che è unilaterale devo fare 1-alpha e non 1-alpha/2 ??

k = dim(pop1)[2]     # number of features = number of tests performed 

which(p.i*k < alpha) 

# perchè estrae il 5????



















