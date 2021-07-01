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
  
  # remove the version numbers of the gene ID's
  rowWithLength18 <- which((nchar(rownames(countData)) == 18) == TRUE)
  rowWithLength17 <- which((nchar(rownames(countData)) == 17) == TRUE)
  
  rownames(countData)[rowWithLength17] <- substr(rownames(countData)[rowWithLength17], 1, 
                                                 nchar(rownames(countData)[rowWithLength17])-2)
  
  rownames(countData)[rowWithLength18] <- substr(rownames(countData)[rowWithLength18], 1, 
                                                 nchar(rownames(countData)[rowWithLength18])-3)
  
  # prepare the metadata
  condition <- as.factor(c(rep("EarlyOnset", times = 47), 
                           rep("LateOnset", times = 174)))
  colData <- data.frame(condition)
  row.names(colData) <- c(samplesToStudy[[1]], samplesToStudy[[2]])
  
  data <- list(countData = countData, colData = colData)
}