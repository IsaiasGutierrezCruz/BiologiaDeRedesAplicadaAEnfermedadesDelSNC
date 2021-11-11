analysisDifferenceExpression <- function(countData, colData, samplesToStudy){
  # Analysis of the data 
  library(edgeR)
  Label <- c(samplesToStudy[[1]], samplesToStudy[[2]])
  colData <- cbind(colData, Label)
  
  y <- DGEList(counts = countData[, 1:221], group = colData$condition)
  colnames(y) <- colData$Label
  
  keep <- rowSums(cpm(y)>1) >= 3
  y <- y[keep, ]
  
  y$samples$lib.size <- colSums(y$counts)
  
  # normalizar 
  y <- calcNormFactors(y)
  
  # estimación de la dispesion
  y <- estimateCommonDisp(y, verbose = TRUE)
  y <- estimateTagwiseDisp(y)
  
  plotBCV(y)
  
  # pruebas 
  et <- exactTest(y)
  
  top2 <- topTags(et, n = 18591)
  
  # plotSmear 
  de <- decideTestsDGE(et, p.value = 0.1)
  summary(de)
  
  detags <- rownames(y)[as.logical(de)]
  
  pdf("Plots/plotSmearAnalysis_edgeR.pdf", 
      width = 8, height = 7, 
      bg = "white",
      colormodel = "cmyk",
      paper = "A4")
  
  plotSmear(et, de.tags = detags, main = "plotSmear")
  abline(h=c(-1, 1), col = "blue")
  
  dev.off()
  
  top2
}