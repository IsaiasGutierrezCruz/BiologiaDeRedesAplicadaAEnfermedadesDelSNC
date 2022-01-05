addInformationTo_top2 <- function(table_to_change, info_to_add = c("SYMBOL", "ENTREZID", "GENENAME"), name_of_cols = c("symbol", "entrez", "name"), delete_NA = TRUE){
  # ---- Description ----
  # It add new information to a table in which each row is a gene and the row 
  # names are the code ENSEMBL
  #
  # ---- Parameters ----
  # table_to_change: data frame
  #     Data frame to add the new information. The row names are the names of genes 
  #     in ENSEMBL
  # info_to_add: character
  #     Character vector with information's names to add 
  # name_of_cols: character
  #     Character vector with column names to the columns added 
  # delete_NA: logical 
  #     Logical value to select delete or not the genes without the information added
  # 
  # ---- Returns ----
  # table_to_change: data frame
  #     Data frame with the new columns
  
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  
  for (i in seq_along(info_to_add)){
    table_to_change[[ name_of_cols[i] ]] <- mapIds(org.Hs.eg.db, 
                                                   keys = row.names(table_to_change),
                                                   column = info_to_add[i],
                                                   keytype = "ENSEMBL",
                                                   multiVals = "first")
  }
  
  # remove NA values 
  if (delete_NA){
    good <- complete.cases(table_to_change$entrez)
    table_to_change <- table_to_change[good, ]
  }
  
  table_to_change
}