moduleEnrichment <- function(graphs, modules, names = c("earlyOnset", "lateOnset"), 
                             output_path = "Results/ModuleEnrichmentAndProjection",
                             pvalueThreshold = 0.05){
  # ---- Description ----
  # It calculates the modules significant enriched considering a set of GO terms.
  # It also build and project the bipartite graph. Finally it keep them in .rds, .gml
  # and .txt files
  #
  # ---- Parameters ----
  # graphs: list
  #     A list of graphs 
  # modules: list 
  #     A list with information of the modules detected
  # names: character
  #     Character vector with the names of each graphs 
  # output_path: character
  #     Character with the path where the output will be keep 
  # pvalueThreshold: numeric
  #     A numeric value to use as a threshold in the selection of the GO terms
  
  #3) Module enrichment and bipartite network
  
  source("5ModuleEnrichmentAndProjection/HTSanalyzeR")
  library(tidyverse)
  library(plyr)
  library(data.table)
  library(igraph)
  load("Data/GOs_and_pathways.RData")
  #set variables
  pvalue_threshold = 0.05 
  
  if (!length(graphs) == length(modules)){
    stop("The number of networks and modules are different")
  }
  
  
  for (i in seq_along(graphs)) {
    #get list of communities 
    l_comm = communities( modules[[i]] )
    #in case we are provided with communities already identified for each vertex
    ## tmp_df = get.data.frame(g, "vertices")
    ## tmp_l  = lapply(X = unique(tmp_df$infomap), FUN = function(i){
    ##   filter(tmp_df, infomap==i)$name
    ## })
    
    #Enrich
    
    #enrichment_list_cases = parallel::mclapply(X = seq_along(l_comm_cases), mc.cores = 1, FUN = function(i){
    enrichment_list = lapply(X = seq_along(l_comm), FUN = function(i){
      print(i)
      nomen = names(l_comm)[i]
      
      my_enrichment = multiHyperGeoTest(collectionOfGeneSets = LIST_GO, 
                                        universe = V( graphs[[i]] )$name, 
                                        hits = l_comm[[i]], 
                                        minGeneSetSize = 1, 
                                        pAdjustMethod = "BH", 
                                        verbose = TRUE
      )
      print(nrow(my_enrichment))
      my_2 = tibble::rownames_to_column(as.data.frame(my_enrichment))#%>%filter(Adjusted.Pvalue<10)
      my_3 = cbind(comm = nomen, my_2)
      return(my_3)
      
    })
    
    enrichment_df = ldply(enrichment_list, data.frame)
    
    
    #calculate adjusted pvalue for all communities
    enrichment_df$Adjusted.Pvalue2 = p.adjust(enrichment_df$Pvalue, method = "BH")
    
    write.table(x = enrichment_df, 
                file = paste0(output_path, "/enrichment_df_", names[i] , ".txt"), 
                row.names = FALSE, 
                col.names = TRUE, 
                sep = "\t", 
                quote = FALSE)
    
    #which(enrichment_df_cases$Adjusted.Pvalue<0.05)
    #which(enrichment_df_cases$Adjusted.Pvalue2<0.05)
    #bipartite graph
    enrichment_df <- enrichment_df %>% mutate(globalPvalue = p.adjust(Adjusted.Pvalue)) %>% filter(globalPvalue < 0.05)
    
    saveRDS(enrichment_df, file=paste0(output_path, "/enrichment_df_", names[i], "_filtered.RDS"))
    
    b = graph_from_data_frame(enrichment_df, directed = FALSE)
    
    V(b)$type = TRUE
    V(b)$type[grep(pattern = "GO_", x = V(b)$name)] = FALSE
    
    #add a value that tells us if the comm shows enrichment or not
    #if it has neighbors, it enriches. Else, no
    V(b)$enrich = degree(b)
    V(b)$enrich = ifelse(V(b)$enrich==0, "NO", "YES")
    
    #write out bipartite
    write.graph(b, file = paste0(output_path, "/bipartite_", cases, ".gml"), "gml")
    
    #remove big enrichment list from memory
    rm(enrichment_df)
    
    #4) Enrichment projection
    
    #projection
    bp = bipartite_projection(graph = b, which = TRUE)
    
    write.graph(graph = bp, file = paste0(output_path, "/bp_", cases, ".gml"), "gml")
  }
}