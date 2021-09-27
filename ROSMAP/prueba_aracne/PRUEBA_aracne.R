library(minet)
library(Rgraphviz)

# early Onset 
y$counts[, 1:47]

saveRDS(y$counts[, 1:47], file = "earlyOnset_counts.rds")

saveRDS(mim, file = "mutual_information_earlyOnset.rds")
readRDS(file = "mutual_information_earlyOnset.rds")

earlyOnset <- as.data.frame(t(y$counts[, 1:47]))

mim <- build.mim(dataset = earlyOnset, estimator = "spearman")
relnet_earlyOnset <- clr(mim = mim, skipDiagonal = 1)

plot( as(red_aracne,"graphNEL") )


#late onset 
y$counts[, 48:221]

saveRDS(y$counts[, 48:221], file = "lateOnset_counts.rds")

saveRDS(mim, file = "mutual_information_lateOnset.rds")
readRDS(file = "mutual_information_lateOnset.rds")

lateOnset <- as.data.frame(t(y$counts[, 48:221]))

mim_lateOnset <- build.mim(dataset = lateOnset, estimator = "spearman")
relnet_lateOnset <- clr(mim = mim_lateOnset, skipDiagonal = 1)

plot( as(red_aracne,"graphNEL") )
