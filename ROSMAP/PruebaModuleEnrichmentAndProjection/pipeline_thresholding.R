#libraries
library(stringi)
library(tidyverse)
library(data.table)
library(plyr)
library(HTSanalyzeR)
library(igraph)
source("MapFlow_2.R")
load("GOs_and_pathways.RData")

#set variables
path = "data/Basal_PBCMC_ARSYM.sif"
my_basename = "cases"
pvalue_threshold = 0.05

#my_i = c(0.9999, 0.99999)
my_i = c(0.992, 0.994, 0.996)
# my_i =        c(0.99,
#                 0.998,
#                 0.9982,#percolation value
#                 0.9983,
#                 0.9985,
#                 0.9988,
#                 0.999#,
# #                0.9999,
# #                0.99999
# )
#read the big sif
x = fread(input = path)
x = x[order(-V2)]
x = x[,c(1,3,2)]

for(i in my_i){
  print(i)
  proctor <- proc.time()

  #make paths to write outs
  print("making paths")
  
  my_pattern = paste0(my_basename, "_", i)
  my_pattern = gsub(pattern = "\\.", replacement = "", x = my_pattern)
  my_pattern = stringi::stri_pad_right(str = my_pattern, width = 12, pad = "0")
  
  my_dir_path = paste0("results/", my_pattern)
  dir.create(path = my_dir_path)
  
  my_graph_path      = paste0(my_dir_path, "/", my_pattern, ".gml")
  my_dict_path       = paste0(my_dir_path, "/", my_pattern, "_dict.txt")
  my_gp_path         = paste0(my_dir_path, "/", my_pattern, "_gp.gml")
  my_rds_path        = paste0(my_dir_path, "/", my_pattern, "modules.RDS")
  my_enrichtbl_path  = paste0(my_dir_path, "/", my_pattern, "_EnrichmentAll.txt")
  my_hitNumber_path  = paste0(my_dir_path, "/", my_pattern, "_mappingPathways.txt")
  my_bipartite_path  = paste0(my_dir_path, "/", my_pattern, "_bipartite.gml")
  my_bp_path         = paste0(my_dir_path, "/", my_pattern, "_bp.gml")
  my_shareLinks_path = paste0(my_dir_path, "/", my_pattern, "_sharedLinks.txt")
  
  #get the network from the big sif 
  print("getting the network")
  ## threshold value
  a = quantile(x = x$V2, i)
  
  ##filter sif
  y = x[V2>=a]
  
  #my network
  g = graph_from_data_frame(y, directed = FALSE)
  E(g)$weight = E(g)$V2

  #detect modules 
  print("detecting modules")
  modules = infomap.community(graph = g, nb.trials = 1000)
  
  ##write out network with modules
  V(g)$infomap = membership(modules)
  write.graph(graph = g, 
              file  = my_graph_path, 
              format = "gml")
  
  ##Write out dictionary
  write.table(x = get.data.frame(x = g, "vertices"), 
              file = my_dict_path, 
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
  
  ##writeout modules
  saveRDS(object = modules, file = my_rds_path)
  
  ##module projection
  gp = mapflow(g, modules)
  
  #export 
  write.graph(graph = gp, 
              file = my_gp_path, 
              format = "gml"
              )

  #module enrichment
  print("enrichment start")
    #get list of communities 
    l_comm = communities(modules)
    
    #check how many pathways from my list have genes mapped in the network
    my_check = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = LIST_GO, 
                                             universe = V(g)$name, 
                                             hits = c("A", "B", "C"), 
                                             minGeneSetSize = 1, 
                                             pAdjustMethod = "BH", 
                                             verbose = TRUE
    )
    my_check = nrow(my_check)
    
    writeLines(text = paste0(my_check, " pathways out of ", length(LIST_GO)), 
               con = my_hitNumber_path)
    
    #begin enrichment
    enrichment_list = parallel::mclapply(X = seq_along(l_comm), 
                                               mc.cores = 4, 
                                               FUN = function(i){
      
      
      nomen = names(l_comm)[i]
      
      my_enrichment = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = LIST_GO, 
                                                     universe = V(g)$name, 
                                                     hits = l_comm[[i]], 
                                                     minGeneSetSize = 1, 
                                                     pAdjustMethod = "BH", 
                                                     verbose = TRUE
      )
      print(nrow(my_enrichment))
      my_2 = tibble::rownames_to_column(as.data.frame(my_enrichment))
      my_3 = cbind(comm = nomen, my_2)
      return(my_3)
      
    })
    
    enrichment_df = ldply(enrichment_list, data.frame)
    enrichment_df$Adjusted.Pvalue2 = p.adjust(enrichment_df$Pvalue, method = "BH")
    
    ##write out big enrichment table
    print("writeout enrichment table")
    fwrite(x = enrichment_df, 
                file = my_enrichtbl_path, 
                row.names = FALSE, 
                col.names = TRUE, 
                sep = "\t", 
                quote = FALSE)
    
    
    #bipartite graph
    print("making bipartite")
    b = graph_from_data_frame(enrichment_df, directed = FALSE)
    ##free up memory, remove enrichment df
    rm(enrichment_df)
    
    #filter by pvalue threshold
    b = delete.edges(graph = b, 
                           edges = E(b)[Adjusted.Pvalue2>pvalue_threshold]
    )
    
    
    #add types to nodes
    V(b)$type = TRUE
    V(b)$type[grep(pattern = "GO_", x = V(b)$name)] = FALSE
    
    #add a value that tells us if the comm shows enrichment or not
    #if it has neighbors, it enriches. Else, no
    V(b)$enrich = degree(b)
    V(b)$enrich = ifelse(V(b)$enrich==0, "NO", "YES")
    
    ##writeout bipartite
    print("writing bipartite")
    write.graph(graph = b, 
                file = my_bipartite_path,
                format = "gml"
                )

  #module projection
    print("projecting bipartite")
    bp = bipartite_projection(graph = b, which = TRUE)
  ##writeout projection
    write.graph(graph = bp, 
                file = my_bp_path, 
                format = "gml")

#module GP GB comparison
    print("begin GP BP comparison")
    xxx = as.matrix(get.adjacency(gp))
    yyy = as.matrix(get.adjacency(bp))
    xxx = xxx[order(as.numeric(rownames(xxx))),
                          order(as.numeric(colnames(xxx)))
                          ]
    yyy = yyy[order(as.numeric(rownames(yyy))),
                          order(as.numeric(colnames(yyy)))
                          ]
    zzz = xxx*yyy
    my_replicated_links = which(zzz!=0, arr.ind = TRUE)
#writeout comparison
    write.table(x = my_replicated_links, 
                file = my_shareLinks_path,
                row.names = FALSE, 
                col.names = TRUE, 
                sep = "\t", 
                quote = FALSE)
    
    print(proc.time()-proctor)
}
