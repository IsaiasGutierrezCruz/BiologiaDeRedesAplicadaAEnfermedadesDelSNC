moduleDetection <- function(graphs, names = c("earlyOnset", "lateOnset"), output_path = "Results/ModuleEnrichmentAndProjection"){
  # ---- Description ----
  # It calculate the modules into co-expression networks and generates a table
  # to represent the graphs with the modules and a .rds file to keep the modules
  #
  # ---- Parameters ----
  # graphs: list
  #     A list of graphs 
  # names: character
  #     Character vector with the names of each graphs 
  # output_path: character
  #     Character with the path where the output will be keep 
  # 
  # ---- Returns ----
  # modules_graphsCoexp: list 
  #     A list of lists. The first list has information of the graphs' modules 
  #     and the second list has the graphs with the modules 
  
  #libraries
  library(igraph)
  library(data.table)
  library(tidyverse)
  
  if (!is.list(graphs)){
    stop("The parameter graphs should be a list")
  }
  
  
  if (!length(graphs) == length(names)){
    stop("The number of netwroks and names are different")
  }
  
  modules_graphsCoexp <- list()
  
  for (i in seq_along(graphs)) {
    #detect modules
    modules = infomap.community(graph = graphs[[i]], nb.trials = 1)
    modules_graphsCoexp[[ names[i] ]] <- modules
    #assign module membership to vertices
    V( graphs[[i]] )$infomap = membership(modules)
    
    write.table(x = get.data.frame(x = graphs[[i]] , "vertices"), 
                file = paste0(output_path, "/dict_commInfomap_", names[i],".txt"), 
                row.names = FALSE, 
                col.names = TRUE, 
                sep = "\t", 
                quote = FALSE)
    
    #write out community structures
    saveRDS(object = modules, file = paste0(output_path, "/modules_", names[i], ".RDS"))
  }
  
  result <- list(modules_graphsCoexp = modules_graphsCoexp, graphs_with_modules = graphs)
  
  result
}