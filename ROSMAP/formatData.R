formatData <- function(RNAseqCounts, samplesToStudy){
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
  
  data <- list(countData = countData, colData = colData)
}