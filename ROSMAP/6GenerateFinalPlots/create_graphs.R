create_graphs <- function(plot_directory, directory){
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(ggpubr)
  library(grDevices)
  g <- igraph::read_graph(file=file.path(directory, 'Results/CoexpressionGraphs/graphcoexp_earlyonset.gml'), format = "gml")
  subg_213_earlyOnset <- induced_subgraph(g, which(V(g)$infomap == 213))
  
  modules_cases = infomap.community(graph =  subg_213_earlyOnset, nb.trials = 1)
  #assign module membership to vertices
  V( subg_213_earlyOnset )$infomap = membership(modules_cases)
  
  V(subg_213_earlyOnset)$infomap
  V(subg_213_earlyOnset)$clu <- as.character(V(subg_213_earlyOnset)$infomap)
  
  E(subg_213_earlyOnset)$weight <- seq(ecount(subg_213_earlyOnset))
  
  gt <- as_tbl_graph(subg_213_earlyOnset)
  
  pal <- colorRampPalette(c("#ff0000", "#ff00ee", "#7700ff"))
  net_pal<- pal(6212)
  
  
  gplot<-ggraph(gt,layout = "graphopt")+
    geom_edge_link0(aes(edge_width=weight),edge_colour = "grey66")+
    geom_node_point(aes(fill = clu,size = degree),shape = 21)+
    #geom_node_text(aes(filter = degree >= 26, label = name),family="serif")+
    geom_node_text(aes(label = name),family="serif")+
    #scale_fill_manual(values = net_pal)+
    #scale_fill_brewer(palette = "Dark2")+
    scale_color_viridis_d()+
    scale_edge_width(range = c(0.2,2))+
    scale_size(range = c(8,25))+
    theme_graph()+
    theme(legend.position = "none")
  
  ggsave(filename = file.path(directory, plot_directory, "GraphCluster213_earlyOnset.png"), plot = gplot, units = c("cm"),
         width = 30, height = 20, dpi = 200, limitsize = FALSE)
  
  
  #late Onset
  g2 <- igraph::read.graph(file=file.path(directory, "Results/CoexpressionGraphs/graphcoexp_lateOnset.gml"), format = "gml")
  
  #subg <- subgraph(g, which(V(g)$infomap == 213))
  subg_24and211_lateOnset <- induced_subgraph(g2, which(V(g2)$infomap == 24 | V(g2)$infomap == 211))
  
  modules_cases = infomap.community(graph =  subg_24and211_lateOnset, nb.trials = 1)
  #assign module membership to vertices
  V( subg_24and211_lateOnset )$infomap = membership(modules_cases)
  
  V(subg_24and211_lateOnset)$infomap
  V(subg_24and211_lateOnset)$clu <- as.character(V(subg_24and211_lateOnset)$infomap)
  
  E(subg_24and211_lateOnset)$weight <- seq(ecount(subg_24and211_lateOnset))
  
  gt2 <- as_tbl_graph(subg_24and211_lateOnset)
  
  gplot2<-ggraph(gt2, layout = 'lgl', maxiter = 500, maxdelta = 0.1, area = 5000, coolexp = 1.5, repulserad = 10)+
    geom_edge_link0(aes(edge_width=weight),edge_colour = "grey66")+
    geom_node_point(aes(fill = clu,size = degree),shape = 21)+
    #geom_node_text(aes(filter = degree >= 26, label = name),family="serif")+
    geom_node_text(aes(label = name),family="serif")+
    #scale_fill_manual(values = net_pal)+
    #scale_fill_brewer(palette = "Dark2")+
    scale_color_viridis_d()+
    scale_edge_width(range = c(0.2,2))+
    scale_size(range = c(8,25))+
    theme_graph()+
    theme(legend.position = "none")
  
  ggsave(filename = file.path(directory, plot_directory, "GraphCluster24and211_lateOnset.png"), plot = gplot2, units = c("cm"),
         width = 40, height = 30, dpi = 200, limitsize = FALSE)
}