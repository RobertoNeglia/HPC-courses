setwd("~/shared-folder/HPC/APPSTAT/Exams/2023-07-07/E1")
rm(list = ls())

products <- read.table("products.txt", header = TRUE)
head(products)
dim(products)

attach(products)
boxplot(
    products,
    col = rainbow(12)
)

par(mfrow = c(1, 2))
svg("boxplot.svg", width = 6, height = 5)
boxplot(
    scale(products, center = TRUE, scale = FALSE),
    col = rainbow(12)
)
dev.off()

svg("boxplot_scaled.svg", width = 6, height = 5)
boxplot(
    scale(products, center = TRUE, scale = TRUE),
    col = rainbow(12)
)
dev.off()

pc.products <- princomp(tourists, scores = TRUE, )
summary(pc.products)

plot(cumsum(pc.products$sd^2) / sum(pc.products$sd^2), type = "l")

svg("screeplot.svg", width = 6, height = 5)
layout(matrix(c(2, 3, 1, 3), 2, byrow = TRUE))
plot(pc.products,
    las = 2,
    main = "Principal components",
)
barplot(
    sapply(products, sd)^2,
    las = 2,
    main = "Original Variables",
    ylab = "Variances"
)
plot(
    cumsum(pc.products$sd^2) / sum(pc.products$sd^2),
    type = "b",
    axes = FALSE,
    xlab = "number of components",
    ylab = "contribution to the total variance",
    ylim = c(0, 1)
)
abline(h = 1, col = "blue")
abline(h = 0.8, lty = 2, col = "blue")
box()
axis(2, at = 0:10 / 10, labels = 0:10 / 10)
axis(
    1,
    at = 1:ncol(products),
    labels = 1:ncol(products),
    las = 2
)
dev.off()

svg("plot.svg", width = 10, height = 10)
products.loads <- pc.products$loadings
par(mfrow = c(3, 1))
for (i in 1:3) {
    barplot(products.loads[, i], ylim = c(-1, 1))
}
dev.off()

plot(pc.products$scores[, 1:2], type = "n", xlab = "PC1", ylab = "PC2")

svg("biplot.svg", width = 10, height = 10)
biplot(pc.products, cex = 0.7)
dev.off()
