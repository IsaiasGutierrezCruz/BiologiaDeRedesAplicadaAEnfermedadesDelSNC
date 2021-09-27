library(minet)
library(Rgraphviz)

# early Onset 
data$countData[, 1:47]
write.table(data$countData[, 1:47], file = "earlyOnset_counts.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)


write.csv(mim, file = "mutual_information_earlyOnset.csv", row.names = TRUE, col.names = TRUE)
write.table(mim, file = "m.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)


earlyOnset <- as.data.frame(t(data$countData[, 1:47]))

mim <- build.mim(dataset = earlyOnset, estimator = "spearman")
red_aracne <- aracne(mim = mim, eps = 0)

plot( as(red_aracne,"graphNEL") )


#late onset 
data$countData[, 48:221]
write.table(data$countData[, 48:221], file = "lateOnset_counts.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)

mim <- build.mim(dataset = data$countData[, 48:221], estimator = "spearman")

red_aracne <- aracne(mim = mim, eps = 0)

plot( as(red_aracne,"graphNEL") )


prueba <- matrix(rnorm(10000), nrow = 100)
mi_quantil <- quantile(prueba, 0.9)
mi_quantil
prueba2 <- ifelse(prueba >= mi_quantil, 1, 0)
prueba[1:10, 1:10]
prueba2[1:10, 1:10]


# ----- normalized ------------

# early Onset 
y$counts[, 1:47]
write.table(y$counts[, 1:47], file = "earlyOnsetNormalized_counts.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)

mim <- build.mim(dataset = y$counts[, 1:47], estimator = "spearman")
red_aracne <- aracne(mim = mim, eps = 0)

plot( as(red_aracne,"graphNEL") )


#late onset 
y$counts[, 48:221]
write.table(y$counts[, 48:221], file = "lateOnsetNormalized_counts.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)

mim <- build.mim(dataset = y$counts[, 48:221], estimator = "spearman")

red_aracne <- aracne(mim = mim, eps = 0)

plot( as(red_aracne,"graphNEL") )

