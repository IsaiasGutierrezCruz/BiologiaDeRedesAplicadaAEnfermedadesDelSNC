#' Create Venn Diagrams
#'
#' This function generates Venn diagrams based on the coexpression graph data.
#'
#' @param plot_directory The directory where the generated plots will be saved.
#' @param directory The directory where the coexpression graph data is located.
#'
#' @return None
#'
#' @import igraph
#' @import tidygraph
#' @import ggraph
#' @import dplyr
#' @import ggpubr
#' @import grDevices
#' @import VennDiagram
#'
#' @examples
#' create_venn_diagrams(plot_directory = "Plots", directory = "ROSMAP")
#'
#' @export
create_venn_diagrams <- function(plot_directory, directory){

  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(ggpubr)
  library(grDevices)
  library(VennDiagram)
  
  g <- igraph::read_graph(file=file.path(directory, 'Results/CoexpressionGraphs/graphcoexp_earlyonset.gml'), format = "gml")
  
  #late Onset
  g2 <- igraph::read_graph(file=file.path(directory, "Results/CoexpressionGraphs/graphcoexp_lateOnset.gml"), format = "gml")
  
  subg_5_earlyOnset <- induced_subgraph(g, which(V(g)$infomap == 5))
  subg_2_earlyOnset <- induced_subgraph(g, which(V(g)$infomap == 2))
  
  subg_1_lateOnset <- induced_subgraph(g2, which(V(g2)$infomap == 1))
  subg_2_lateOnset <- induced_subgraph(g2, which(V(g2)$infomap == 2))
  
  names_spanish <- c("Módulo 5 Inicio temprano", "Módulo 2 Inicio temprano" ,  "Módulo 1 Inicio tardío", "Módulo 2 Inicio tardío")
  venn.diagram(
    x = list(names(V(subg_5_earlyOnset)), names(V(subg_2_earlyOnset)), names(V(subg_1_lateOnset)), names(V(subg_2_lateOnset))),
    category.names = c("Module 5 Early Onset", "Module 2 Early Onset" ,  "Module 1 Late Onset", "Module 2 Late Onset"),
    #filename = 'Plots/allclusters_venn_diagramEng.png',
    filename = file.path(directory,plot_directory,  'allclusters_venn_diagramEng.png'),
    output=TRUE,
    # Output features
    imagetype="png" ,
    height = 1500 , 
    width = 1500 , 
    resolution = 200,
    compression = "lzw",
    # Circles
    lwd = 1,
    col=c("#440154ff", '#21908dff', '#fde725ff', '#3f9dfcff'),
    #lty = 'blank',
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha('#3f9dfcff', 0.3)),
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(-27, 27, 150),
    #cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    #rotation = 1
  )
  
  subg_213_earlyOnset <- induced_subgraph(g, which(V(g)$infomap == 213))
  
  modules_cases = infomap.community(graph =  subg_213_earlyOnset, nb.trials = 1)
  #assign module membership to vertices
  V( subg_213_earlyOnset )$infomap = membership(modules_cases)
  
  V(subg_213_earlyOnset)$infomap
  V(subg_213_earlyOnset)$clu <- as.character(V(subg_213_earlyOnset)$infomap)
  
  E(subg_213_earlyOnset)$weight <- seq(ecount(subg_213_earlyOnset))

  subg_211_lateOnset <- induced_subgraph(g2, which(V(g2)$infomap == 211))
  
  # cluster 24
  subg_24_lateOnset <- induced_subgraph(g2, which(V(g2)$infomap == 24))
  
  modules_cases = infomap.community(graph =  subg_24_lateOnset, nb.trials = 1)
  #assign module membership to vertices
  V( subg_24_lateOnset )$infomap = membership(modules_cases)
  
  V(subg_24_lateOnset)$infomap
  V(subg_24_lateOnset)$clu <- as.character(V(subg_24_lateOnset)$infomap)
  
  E(subg_24_lateOnset)$weight <- seq(ecount(subg_24_lateOnset))
  
  # comparison between the 213 module early onset and the modules 211 and 24 late onset
  venn.diagram(
    x = list(names(V(subg_213_earlyOnset)), names(V(subg_211_lateOnset)), names(V(subg_24_lateOnset))),
    category.names = c("Module 213 Early Onset" , "Module 211 Late Onset" , "Module 24 Late Onset"),
    #filename = 'Plots/inmunesys_allclusters_venn_diagramEng.png',
    filename = file.path(directory, plot_directory, 'inmunesys_allclusters_venn_diagramEng.png'),
    output=TRUE,
    # Output features
    imagetype="png" ,
    height = 1500 , 
    width = 1500 , 
    resolution = 200,
    compression = "lzw",
    # Circles
    lwd = 1,
    col=c("#440154ff", '#21908dff', '#fde725ff'),
    #lty = 'blank',
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 150),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  

}



