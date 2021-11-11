mergeCounts <- function(counts = "normalized"){
  # get the data of interest 
  if (counts == "normalized"){
    # RNAseq plates 1 - 6
    RNAseq16 <- read.table(file = 'Data/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv') 
    # RNAseq plates 7 - 8
    RNAseq78 <- read.table(file = 'Data/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv') 
    
    # merge the files 
    RNAseqCounts <- cbind(RNAseq16, RNAseq78[, c(-1, -2)])  
  } else if(counts == "un-normalized"){
    RNAseqCounts <- read.table(file = 'Data/ROSMAP_RNAseq_FPKM_gene.tsv')
  }
  
  
  # Add the names of the columns and remove the first row 
  colnames(RNAseqCounts) <- RNAseqCounts[1, ]
  RowsNewRNAseqCounts <- c(FALSE, rep(TRUE, times = nrow(RNAseqCounts) - 1))
  RNAseqCountsReady <- RNAseqCounts[RowsNewRNAseqCounts, ]
  
  # Remove Trialing undersocre and [0-8]
  colnamesToChange <- (colnames(RNAseqCountsReady) %in% c("tracking_id","gene_id") )==F
  colnames(RNAseqCountsReady)[  colnamesToChange ] <- substr(colnames(RNAseqCountsReady)[  colnamesToChange ], 1, 
                                                         nchar(colnames(RNAseqCountsReady)[  colnamesToChange ])-2)
  RNAseqCountsReady
}