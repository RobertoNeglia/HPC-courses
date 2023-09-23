X <- read.table("~/GitHub/Applied-Statistics-Exam/Exams of previous years/2018/2018-09-13/Waiting.txt") # example, substitute
attach(X)
x11()
boxplot(waiting ~ course + city, las = 2)

# a) Two-ways anova and hypothesis checking
an <- aov(waiting ~ course + city + course:city)
summary(an)

# assumptions are gaussianity of each group and same variance
# verification of assumptions:
shapiro.test(an$residuals)
x11()
qqnorm(an$residuals)
qqline(an$residuals)
# we can assume gaussianity IF we can also assure homoschedasticity amongst groups
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Iasi" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Iasi" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Iasi" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Bucarest" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Bucarest" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Bucarest" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Main")])
# or more simply:
fact <- with(X, interaction(course, city))
bartlett.test(waiting, fact) # H0: same variance

# b) Model reduction
an_new <- aov(waiting ~ course + course:city)
summary(an_new)

# c) Confidence intervals (balanced) for the difference of the means of the groups identified in previous step
n <- dim(X)[1]
g <- 3
k <- g * (g - 1) / 2 # num of comparisons, for bonferroni correction
smean <- tapply(waiting, course, mean) # in alphabetical order of the factors
Spooled <- sum((an_new$residuals)^2) / (n - g)
n1 <- length(which(course == "Starter"))
n2 <- length(which(course == "Main"))
n3 <- length(which(course == "Dessert")) # balanced
n_balanced <- n3
alphaB <- 0.05 / k
# smean_a - smean_b ~ N(mu_a - mu_b, simga^2*(1/n_a + 1/n_b)) = N(mu_a - mu_b, 2*simga^2/n)
# (smean_a - smean_b)/sqrt(2*sigma^2/n) ~ N(0,1)
# s = SS_res/(n-g) ~ sigma^2*Chi_square(n-g)
# (smean_a - smean_b)/sqrt(2*s/n) ~ t(n-g)
band <- qt(1 - alphaB / 2, an_new$df) * sqrt(Spooled * 2 / n_balanced)
sdeltamean <- c(smean[1] - smean[2], smean[2] - smean[3], smean[3] - smean[1])
CI1 <- cbind(
    inf = sdeltamean - band,
    center = sdeltamean,
    sup = sdeltamean + band
)
CI1
SS <- c(sum((an_new$residuals[which(course == "Main")])^2), sum((an_new$residuals[which(course == "Starter")])^2), sum((an_new$residuals[which(course == "Dessert")])^2))
chi_up <- qchisq(1 - alphaB / 2, n_balanced - 1)
chi_down <- qchisq(alphaB / 2, n_balanced - 1)
CI2 <- cbind(
    inf = SS / chi_up,
    center = SS / (n_balanced - 1),
    sup = SS / chi_down
)
rownames(CI2) <- c("Main", "Starter", "Dessert")
CI2
