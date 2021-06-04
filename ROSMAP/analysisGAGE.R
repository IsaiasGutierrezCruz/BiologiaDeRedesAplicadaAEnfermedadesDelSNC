analysisGAGE <- function(countData, samplesToStudy){
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library(gage)
  library(gageData)
  
  # remove the version numbers of the gene ID's
  rowWithLength18 <- which((nchar(rownames(countData)) == 18) == TRUE)
  rowWithLength17 <- which((nchar(rownames(countData)) == 17) == TRUE)
  
  rownames(countData)[rowWithLength17] <- substr(rownames(countData)[rowWithLength17], 1, 
                                                 nchar(rownames(countData)[rowWithLength17])-2)
  
  rownames(countData)[rowWithLength18] <- substr(rownames(countData)[rowWithLength18], 1, 
                                                 nchar(rownames(countData)[rowWithLength18])-3)
  
  
  # change the format of the ID's 
  nombres <- data.frame(id = rownames(countData)) 
  nombres$entrez <- mapIds(org.Hs.eg.db,
                           keys=rownames(countData), 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
  
  rownames(countData) <- nombres$entrez
  
  
  # remove NA values 
  good <- complete.cases(rownames(countData))
  countData <- countData[good, ]
  
  # load the pathways
  data(kegg.sets.hs)
  data(sigmet.idx.hs)
  kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
  kegg.sets.hs <- kegg.sets.hs[1:131]
  
  # groups 
  earlyOnset <- which(colnames(countData)%in%samplesToStudy$LessThan70Years)
  lateOnset <- which(colnames(countData)%in%samplesToStudy$GreaterThan70Years)
  
  earlyOnset_v_LateOnsetSAMEDIR <- gage(exprs = countData, 
                                        gsets = kegg.sets.hs, 
                                        ref = earlyOnset,
                                        samp = lateOnset, 
                                        compare = "unpaired",
                                        same.dir = TRUE)
  
  earlyOnset_v_LateOnset.SigSAMEDIR <- sigGeneSet(earlyOnset_v_LateOnsetSAMEDIR)
  
  earlyOnset_v_LateOnsetBOTHDIR <- gage(exprs = countData, 
                                        gsets = kegg.sets.hs, 
                                        ref = earlyOnset,
                                        samp = lateOnset, 
                                        compare = "unpaired",
                                        same.dir = FALSE)
  
  earlyOnset_v_LateOnset.SigBOTHDIR <- sigGeneSet(earlyOnset_v_LateOnsetBOTHDIR, 
                                                  outname = "earlyOnset_v_lateOnset")
  earlyOnset_v_LateOnset.SigBOTHDIR
}