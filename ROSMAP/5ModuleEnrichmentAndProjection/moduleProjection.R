moduleProjection <- function(graphs, modules, names = c("earlyOnset", "lateOnset"), 
                             output_path = "Results/ModuleEnrichmentAndProjection"){
  # ---- Description ----
  # It calculates the flow of the graphs' modules and makes gml files of them
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
  
  # Module flow projection
  
  #libraries
  library(igraph)
  library(data.table)
  library(tidyverse)
  library(dplyr)
  source("5ModuleEnrichmentAndProjection/MapFlow_2.R")
  
  
  if (!length(graphs) == length(modules)){
    stop("The number of netwroks and modules are different")
  }
  
  for (i in seq_along(graphs)) {
    #Make projection to the community space
    flow_cases = mapflow(graphs[[i]], modules[[i]])
    #export 
    write.graph(flow_cases, paste0(output_path, "/flow_", names[i], ".gml"), "gml")
  }
}