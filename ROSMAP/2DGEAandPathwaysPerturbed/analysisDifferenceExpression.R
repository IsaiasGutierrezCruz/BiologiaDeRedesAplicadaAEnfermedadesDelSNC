analysisDifferenceExpression <- function(countData, colData, samplesToStudy, 
                                         makePlotBCV = FALSE, makePlotSmear = TRUE,
                                         output_path="Plots"){
  # ---- Description ----
  # It normalize the count data and compare the expression of the genes in each group. 
  # It also make the plots BCV and Smear 
  #
  # ---- Parameters ----
  # countData: data frame
  #     Data frame with the counts from the RNASeq assay. In the columns are represented 
  #     the samples' names and in the rows are represented the genes' names 
  # SamplesToStudy: list
  #     A list that contains the names of samples in control and study group
  # makePlotBCV: logical 
  #     Logical value to indicate if a Plot BCV has to be done
  # makePlotSmear: logical
  #     Logical value to indicate if a Plot Smear has to be done
  # 
  # 
  # ---- Returns ----
  # list: list
  #     A list that contains the data of the gene expression comparison and the 
  #     counts of RNAseq normalized 
  
  
  # Analysis of the data 
  library(edgeR)
  Label <- c(samplesToStudy[[1]], samplesToStudy[[2]])
  colData <- cbind(colData, Label)
  
  y <- DGEList(counts = countData[, 1:ncol(countData)], group = colData$condition)
  colnames(y) <- colData$Label
  
  keep <- rowSums(cpm(y)>1) >= 3
  y <- y[keep, ]
  
  y$samples$lib.size <- colSums(y$counts)
  
  # normalize
  y <- calcNormFactors(y)
  # dispersion estimation
  y <- estimateCommonDisp(y, verbose = TRUE)
  y <- estimateTagwiseDisp(y)
  
  
  # plot BCV
  if (makePlotBCV){
    pdf(paste0(output_path, "/plotBCVAnalysis_edgeR.pdf"), 
        width = 8, height = 7, 
        bg = "white",
        colormodel = "cmyk",
        paper = "A4")
    
    plotBCV(y)
    
    dev.off()
  }
  
  # proofs
  et <- exactTest(y)
  
  top2 <- topTags(et, n = nrow(y$counts))
  
  if (makePlotSmear){
    # plot Smear 
    de <- decideTestsDGE(et, p.value = 0.1)
    print(summary(de))
    
    detags <- rownames(y)[as.logical(de)]
    
    pdf(paste0(output_path, "/plotSmearAnalysis_edgeR.pdf"), 
        width = 8, height = 7, 
        bg = "white",
        colormodel = "cmyk",
        paper = "A4")
    
    plotSmear(et, de.tags = detags, main = "plotSmear")
    abline(h=c(-1, 1), col = "blue")
    
    dev.off()
  }
  
  
  list(top2 = top2, countDataNormalized = y$counts)
}