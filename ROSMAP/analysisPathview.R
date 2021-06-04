analysisPathview <- function(top2, earlyOnset_v_LateOnset.SigBOTHDIR){
  library(dplyr)
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library(pathview)
  
  rowWithLength18 <- which((nchar(rownames(top2)) == 18) == TRUE)
  rowWithLength17 <- which((nchar(rownames(top2)) == 17) == TRUE)
  
  rownames(top2$table)[rowWithLength17] <- substr(rownames(top2$table)[rowWithLength17], 1, 
                                                  nchar(rownames(top2$table)[rowWithLength17])-2)
  
  rownames(top2$table)[rowWithLength18] <- substr(rownames(top2$table)[rowWithLength18], 1, 
                                                  nchar(rownames(top2$table)[rowWithLength18])-3)
  
  
  top2$table$symbol <- mapIds(org.Hs.eg.db, 
                              keys = row.names(top2),
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")
  
  top2$table$entrez = mapIds(org.Hs.eg.db,
                             keys=row.names(top2), 
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
  
  top2$table$name =   mapIds(org.Hs.eg.db,
                             keys=row.names(top2), 
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first")
  
  
  top2$table
  
  # remove NA values 
  good <- complete.cases(top2$table$entrez)
  top2$table <- top2$table[good, ]
  
  
  library(gage)
  library(gageData)
  library(pathview)
  
  foldChanges <- top2$table$logFC
  names(foldChanges) <- top2$table$entrez
  
  # Get the pathways
  keggrespathways <- data.frame(id=rownames(earlyOnset_v_LateOnset.SigBOTHDIR$greater), 
                                earlyOnset_v_LateOnset.SigBOTHDIR$greater) %>% 
    tibble::as_tibble() %>% 
    filter(row_number()<=5) %>% 
    .$id %>% 
    as.character()
  
  # get the ID's
  keggresids <- substr(keggrespathways, start=1, stop=8)
  
  
  setwd("Graficos")
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldChanges, pathway.id=pid, species="hsa"))
  tmp
}