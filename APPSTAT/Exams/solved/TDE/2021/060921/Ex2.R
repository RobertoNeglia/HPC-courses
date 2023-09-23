### Ex2

dataset = read.table("waterquality.txt", header = T)


pop1 = dataset[1:8,]
pop2 = dataset[9:18,]
n1 <- dim(pop1)[1] 
n2 <- dim(pop2)[1] 
p  <- dim(pop1)[2] 


# we compute the sample mean, covariance matrices and the matrix Spooled
# we are assuming that the true covariance matrix is the same 

pop1.mean <- sapply(pop1,mean)
pop2.mean <- sapply(pop2,mean)
pop1.cov  <-  cov(pop1)
pop2.cov  <-  cov(pop2)
Sp      <- ((n1-1)*pop1.cov + (n2-1)*pop2.cov)/(n1+n2-2)

load("mcshapiro.test.RData") 
mcshapiro.test(pop1)$pvalue
mcshapiro.test(pop2)$pvalue

# b) BH

# caso non bernuolli metto Sp
z.i <- (pop1.mean-pop2.mean)/sqrt(diag(Sp)*(1/n1+1/n2))

# noi vogliamo che pop1.mean-pop2.mean < 0
# -> H0: pop1.mean-pop2.mean >= 0   vs H1: pop1.mean-pop2.mean < 0
# rifiuto se z.i è troppo piccola -> se z.i < - z(1-alpha)
# pvalue = pnorm(z.i)

# caso bilaterale 
# p.i <- ifelse(z.i<0, 2*pnorm(z.i),2*(1-pnorm(z.i)))  #p value

p.i <- pnorm(z.i)

# all of this is done with vector (but it's like we have done a for loop and one at the time)
alpha = 0.05
which(p.i< alpha)# we see one at the time which i are significant


k <- 7

which(p.i*k<alpha)    # -> 2,3,5  con il bilaterale
                      # -> 3,5 con unilaterale

# or we can use this function and choose the time of correction

p.Bf <- p.adjust(p.i, method='bonferroni')

which(p.Bf< alpha)     # -> 2,3,5 con bilaterale
                       # -> 3,5 con unilaterale


# Benjamini-Hockberg (control the false discovery rate)  
p.BH <- p.adjust(p.i, method='BH')

which(p.BH < alpha)    # -> 2,3,5,7  caso bilaterale 
                       # 3,5 caso unilaterale  


# c) Bonf
alpha <- 0.05
k = p
IC <- cbind(pop1.mean-pop2.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k), n1+n2-2),
            pop1.mean-pop2.mean,
            pop1.mean-pop2.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k), n1+n2-2))
IC
# Per unilaterale devo fare 1-alpha e non 1-alpha/2 ??

# 2,3,5 per il biaterale -> ok
# 3,5 per l'unilaterale -> ok
























