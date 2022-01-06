analysisKeggGraphs <- function(pathways_names, objects_names, species_name, database, etiquetas, output_path = "Plots/"){
  # ---- Description ----
  # It analyze specific Kegg networks to calculate their network properties and make plots considering the nodes' degree
  #
  # ---- Parameters ----
  # pathways_names: character
  #     Character vector that keep the pathways names 
  # objects_names: character
  #     Character vector with abbreviations of the pathways names
  # species_name: character
  #     Character with the name of the specie of interest 
  # database: character
  #     Character with the database name of interest
  # etiquetas: character 
  #     Character with the labels to the plots
  # output_path: character
  #     Character with the path where the output will be keep 
  # 
  # ---- Returns ----
  # prop_graphs: list
  #     A list with the properties of the network (index = 1) and a list with the 
  #     networks with the properties (index = 2)
  
  
  library(graphite)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(ggpubr)
  source("supportFunctions/calculateGraphProperties.R")
  humanKegg <- pathways(species_name, database) 
  
  prop_graphs <- list()
  for(i in seq_along(pathways_names)){
    
    graph <- humanKegg[[pathways_names[i]]]
    graphNEL <- pathwayGraph(graph)
    graphIgraph <- graph_from_graphnel(graphNEL)
    
    prop <- calculateGraphProperties(networks = list(graphIgraph), names = c(objects_names[i]), 
                             calculate_communities = FALSE, keep_network_info = TRUE)
    
    prop_graphs[[i]] <- prop[[1]]
    graphIgraph <- prop[[2]][[1]]
    
    # ---- PLOTS ----
    gt <- as_tbl_graph(graphIgraph)
    
    Plotg <- gt %>% 
      ggraph(layout = 'kk') + 
      geom_edge_link(aes(), show.legend = FALSE) + 
      geom_node_point(aes(size = degree), colour = "lightskyblue1") +
      theme_graph() + ggtitle(etiquetas[i])
    
    ggsave(filename = paste0(output_path, "Graph", pathways_names[i], ".png"), plot = Plotg, units = c("cm"),
           width = 30, height = 20, dpi = 200, limitsize = FALSE)
  }
  
  prop_graphs
}