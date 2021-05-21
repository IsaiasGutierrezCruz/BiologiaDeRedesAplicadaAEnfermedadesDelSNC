# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts()

# get the metadata of the necessary samples for the study
source("getSamplesToStudy.R")
samplesToStudy <- getSamplesToStudy(RNAseqCounts = RNAseqCounts, 
                                    FullMeta = FullMeta)

# data analysis with DESeq2

# necesita RNAseqCounts y samplesToStudy

# load DESeq2
library(DESeq2)

# prepare countData 
samplesLess <- RNAseqCounts[, samplesToStudy[[1]]]
samplesGreater <- RNAseqCounts[, samplesToStudy[[2]]]
countData <- cbind(samplesLess, samplesGreater)
row.names(countData) <- RNAseqCounts$gene_id
countData <- as.matrix(countData)
countData <- apply(countData, 2, as.integer)
rownames(countData) <- RNAseqCounts$gene_id

countData <- countData[rowSums(countData)>1, ]

# prepare the metadata
condition <- as.factor(c(rep("EarlyOnset", times = 47), 
                         rep("LateOnset", times = 174)))
colData <- data.frame(condition)
row.names(colData) <- c(samplesToStudy[[1]], samplesToStudy[[2]])
colData <- as.matrix(colData)

# use of the package DESeq2
dds = DESeqDataSetFromMatrix(countData = countData, 
                             colData = colData, 
                             design =~condition)
dds = DESeq(dds)
dds

# results 
res = results(dds, contrast = c("condition", "EarlyOnset", "LateOnset"))
res = res[order(res$pvalue), ]
summary(res)

res$pvalue
res$padj
res$log2FoldChange

metadata(res)$filterThreshold

plotMA(res)


# ---- prueba edgeR ----
library(edgeR)
Label <- c(samplesToStudy[[1]], samplesToStudy[[2]])
colData <- cbind(colData, Label)

y <- DGEList(counts = countData[, 1:221], group = colData$condition)
colnames(y) <- colData$Label
dim(y)
head(y)

keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, ]
dim(y)

y$samples$lib.size
y$samples$lib.size <- colSums(y$counts)

# normalizar 
y <- calcNormFactors(y)
y$samples

# exploración de datos 
plotMDS(y)

# estimación de la dispesion
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)

plotBCV(y)

# pruebas 
et <- exactTest(y)
et

top <- topTags(et)
top
top$adjust.method
top2 <- topTags(et, n = 18591)
top2

# histograma 
table(top2$table$FDR < 0.05)

table(top2$table$FDR <0.05)/nrow(top2$table)

hist(top2$table$FDR, breaks = 100, main = "Histograma de FDR")
abline(v = 0.05, col = "red", lwd = 3)

# plotSmear 
de <- decideTestsDGE(et)
summary(de)

detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags = detags, main = "plotSmear")
abline(h=c(-1, 1), col = "blue")

