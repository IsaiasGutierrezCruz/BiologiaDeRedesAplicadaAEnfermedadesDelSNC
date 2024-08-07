graphProps <- function(mimgraph, mimgraph_lo){
  library(igraph)
  library(dplyr)
  
  
  prop_graphs <- list()
  
  # ---- earlyOnset ----
  # components
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(components = components(mimgraph)))
  # average path length
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(average.path.length = average.path.length(mimgraph)))
  # Clustering coefficient global 
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(transitivity = transitivity(mimgraph)))
  # Diameter 
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(diameter = diameter(mimgraph)))
  
  # nodes properties analysis
  # degree
  V(mimgraph)$degree <- degree(mimgraph)
  # betweenness
  V(mimgraph)$betweenness <- betweenness(mimgraph)
  # get the local clustering coefficient 
  V(mimgraph)$transitivity <- transitivity(mimgraph, type = "local", isolates = "zero")
  
  # comunidades 
  comm.louvain <- cluster_louvain(mimgraph)
  V(mimgraph)$comm.louvain <- membership(comm.louvain)
  
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(nodes = get.data.frame(x = mimgraph, what = "vertices")))
  
  
  # edges 
  E(mimgraph)$betweenness <- edge.betweenness(mimgraph)
  
  prop_graphs[["earlyOnset"]] <- append(prop_graphs[["earlyOnset"]], list(edges = get.data.frame(x = mimgraph, what = "edges")))
  
  # ------ lateonset ------
  
  # components
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(components = components(mimgraph_lo)))
  # average path length
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(average.path.length = average.path.length(mimgraph_lo)))
  # Clustering coefficient global 
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(transitivity = transitivity(mimgraph_lo)))
  # Diameter 
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(diameter = diameter(mimgraph_lo)))
  
  # nodes properties analysis
  # degree
  V(mimgraph_lo)$degree <- degree(mimgraph_lo)
  # betweenness
  V(mimgraph_lo)$betweenness <- betweenness(mimgraph_lo)
  # get the local clustering coefficient 
  V(mimgraph_lo)$transitivity <- transitivity(mimgraph_lo, type = "local", isolates = "zero")
  
  # comunidades 
  comm.louvain <- cluster_louvain(mimgraph_lo)
  V(mimgraph_lo)$comm.louvain <- membership(comm.louvain)
  
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(nodes = get.data.frame(x = mimgraph_lo, what = "vertices")))
  
  # edges 
  E(mimgraph_lo)$betweenness <- edge.betweenness(mimgraph_lo)
  
  prop_graphs[["lateOnset"]] <- append(prop_graphs[["lateOnset"]], list(edges = get.data.frame(x = mimgraph_lo, what = "edges")))
  
  prop_graphs
}