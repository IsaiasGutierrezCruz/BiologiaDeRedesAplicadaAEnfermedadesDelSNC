addInformationTo_top2 <- function(top2){
  library("AnnotationDbi")
  library("org.Hs.eg.db")
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
  
  # remove NA values 
  good <- complete.cases(top2$table$entrez)
  top2$table <- top2$table[good, ]
  
  top2
}