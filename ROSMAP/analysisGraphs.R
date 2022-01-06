analysisGraphs <- function(pathways_names, objects_names, species_name, database, etiquetas){
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
    
    ggsave(filename = paste0("Plots/Graph", pathways_names[i], ".png"), plot = Plotg, units = c("cm"),
           width = 30, height = 20, dpi = 200, limitsize = FALSE)
  }
  #names(prop_graphs) <- objects_names
  
  prop_graphs
}