relevanceNetworks <- function(dataCount){
  # ---- Description ----
  # Build the co-expression networks from counts RNAseq with mutual information 
  # relevance networks
  #
  # ---- Parameters ----
  # dataCount: data frame 
  #     Data frame with the counts from the RNASeq assay. In the columns are represented 
  #     the samples' names and in the rows are represented the genes' names 
  # 
  # ---- Returns ----
  # mimgraph: list 
  #     A list with information of the structure of the output network
  
  library(minet)
  library(Rgraphviz)
  library(igraph)
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  
  # add key symbol as row names
  nombresDataCount <- mapIds(org.Hs.eg.db, 
                            keys = row.names(dataCount),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
  nombresDataCount <- as.vector(nombresDataCount)
  row.names(dataCount) <- nombresDataCount
  # remove genes without entrez ID
  good <- complete.cases(row.names(dataCount))
  dataCount <- dataCount[good, ]
  dataCount <- dataCount[unique(row.names(dataCount)), ]
  
  # prepare the format of data to calculate mutual information
  dataCount <- as.data.frame(t(dataCount))
  # calculate mutual information
  mim <- build.mim(dataset = dataCount, estimator = "spearman")

  # relevant network 
  umbral <- quantile(mim, 0.99, na.rm = TRUE)
  mim <- ifelse(mim < umbral, 0, 1)
  #rm(mim)
  mim[is.na(mim)] <- 0
  mimgraph <- graph_from_adjacency_matrix(mim, mode = "undirected")
  
  mimgraph
}