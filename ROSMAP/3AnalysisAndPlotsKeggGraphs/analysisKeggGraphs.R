analysisKeggGraphs <- function(pathways_names, objects_names, species_name, database, etiquetas, output_path = "Plots/"){
  # ---- Description ----
  # It analyze specific Kegg networks to calculate their network properties and make plots of the networks
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
  #     A list with the properties of networks and plots of the networks. It keeps
  #     the plots in the path set up
  
  
  library(graphite)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(ggpubr)
  
  source("supportFunctions/calculateGraphProperties.R")
  humanKegg <- pathways(species_name, database) 
  
  # get the networks
  graphs <- list()
  for(i in seq_along(pathways_names)){
    graph <- humanKegg[[ pathways_names[i] ]]
    graphNEL <- pathwayGraph(graph)
    graphs[[i]] <- graph_from_graphnel(graphNEL)
  }  
    
  # calculate the network properties  
  prop <- calculateGraphProperties(networks = graphs, names = objects_names, 
                             calculate_communities = FALSE, keep_network_info = TRUE)
    
  prop_graphs <- prop$prop_graphs
  graphs <- prop$graph_with_properties
  
  # plot the networks 
  for (i in seq_along(pathways_names)){
    # ---- PLOTS ----
    gt <- as_tbl_graph(graphs[[i]])
    
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