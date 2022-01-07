formatData <- function(RNAseqCounts, samplesToStudy, group_names = c("EarlyOnset", "LateOnset")){
  # ---- Description ----
  # It change the format of the data to do a Differential Gene Expression Analysis and 
  # create factor for each group (control and study)
  #
  # ---- Parameters ----
  # RNASeqCounts: data frame
  #     Data frame with the counts from the RNASeq assay. In the columns are represented 
  #     the samples' names and in the rows are represented the genes' names 
  # SamplesToStudy: list
  #     A list that contains the names of samples in control and study group
  # group_names: character
  #     Vector of strings with the names of the control group and study group
  # 
  # 
  # ---- Returns ----
  # data: list
  #     A list that contains the count data with the correct format to be analyzed (index = 1) 
  #     and the group labels (index = 2)
  
  
  # prepare countData 
  countData <- cbind(RNAseqCounts[, samplesToStudy[[1]]], RNAseqCounts[, samplesToStudy[[2]]])
  countData <- apply(countData, 2, as.integer)
  rownames(countData) <- rownames(RNAseqCounts)
  countData <- countData[rowSums(countData)>1, ]
  
  # section to remove characters of the data set ROSMAP
  # remove the version numbers of the gene ID's
  rowWithLength18 <- which((nchar(rownames(countData)) == 18) == TRUE)
  rowWithLength17 <- which((nchar(rownames(countData)) == 17) == TRUE)
  
  rownames(countData)[rowWithLength17] <- substr(rownames(countData)[rowWithLength17], 1, 
                                                 nchar(rownames(countData)[rowWithLength17])-2)
  
  rownames(countData)[rowWithLength18] <- substr(rownames(countData)[rowWithLength18], 1, 
                                                 nchar(rownames(countData)[rowWithLength18])-3)
  
  # prepare the metadata
  condition <- as.factor(c(rep(group_names[1], times = length(samplesToStudy[[1]])), 
                           rep(group_names[2], times = length(samplesToStudy[[2]]))))
  colData <- data.frame(condition)
  row.names(colData) <- c(samplesToStudy[[1]], samplesToStudy[[2]])
  
  data <- list(countData = countData, colData = colData)
}