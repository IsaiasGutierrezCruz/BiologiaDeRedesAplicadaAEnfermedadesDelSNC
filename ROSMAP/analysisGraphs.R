analysisGraphs <- function(pathways_names, objects_names, species_name, database, etiquetas){
  library(graphite)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(ggpubr)
  humanKegg <- pathways(species_name, database) 
  prop_graphs <- list()
  for(i in seq_along(pathways_names)){
    prop_graphs[[i]] <- list()
    graph <- humanKegg[[pathways_names[i]]]
    graphNEL <- pathwayGraph(graph)
    graphIgraph <- graph_from_graphnel(graphNEL)
    # components
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(components = components(graphIgraph)))
    # average path length 
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(average.path.length = average.path.length(graphIgraph))) 
    # Clustering coefficient global  
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(transitivity = transitivity(graphIgraph, type = "global")))
    # Diameter 
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(diameter = diameter(graphIgraph)))
    # ------------- nodes properties analysis  ------------
    # degree
    V(graphIgraph)$degree <- degree(graphIgraph)
    # betweenness
    V(graphIgraph)$betweenness <- betweenness(graphIgraph)
    # get the local clustering coefficient 
    V(graphIgraph)$transitivity <- transitivity(graphIgraph, type = "local", isolates = "zero")
    # get the full nodes information
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(nodes = get.data.frame(x = graphIgraph, what = "vertices")))
    # ---------- analysis of the properties of the edges ------------
    E(graphIgraph)$intermediacion <- edge.betweenness(graphIgraph)
    prop_graphs[[i]] <- append(prop_graphs[[i]], list(edges = get.data.frame(x = graphIgraph, what = "edges")))
    
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
  names(prop_graphs) <- objects_names
  
  prop_graphs
}