create_degree_plots <- function(plot_directory, directory){
  library(ggplot2)
  library(ggthemes)
  # graph flow
  prop_graphs_flow <- readRDS(file.path(directory, "Results/ModuleEnrichmentAndProjection/prop_graphs_flow.RDS"))
  dataHistProj <- table(prop_graphs_flow[[1]]$nodes$degree)
  dataProj <- as.data.frame(dataHistProj)
  
  dataHistProj2 <- table(prop_graphs_flow[[2]]$nodes$degree)
  dataProj2 <- as.data.frame(dataHistProj2)
  
  plot_degree_proj <- data.frame(num_nodes = c(dataProj$Freq, dataProj2$Freq), 
                                 Degree = c(as.numeric(dataProj$Var1), 
                                            as.numeric(dataProj2$Var1)), 
                                 type = c(rep("Early Onset", each = nrow(dataProj)), 
                                          rep("Late Onset", each = nrow(dataProj2))))
  p3 <- ggplot(data = plot_degree_proj, mapping = aes(x= Degree, y = num_nodes, col = type)) + 
    geom_point(alpha = 0.6, size = 2) + scale_x_log10() + scale_y_log10() + 
    ggthemes::scale_color_ptol() + ylab("Number of nodes") + xlab("Degree")+
    theme(legend.position = "top", 
          panel.background = element_rect(fill = "white", colour="black"),
          panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
          panel.grid.minor = element_blank())
  ggsave(filename = file.path(directory, plot_directory, 'Degree_graph_flow.png'), plot = p3, units = c("cm"),
         width = 30, height = 25, dpi = 200, limitsize = FALSE)
  
  # relevance networks
  prop_graphs_relevance_networks <- readRDS(file.path(directory, "Results/CoexpressionGraphs/prop_graphsCoexp.rds"))
  dataHist <- table(prop_graphs_relevance_networks[[1]]$nodes$degree)
  prueba <- as.data.frame(dataHist)
  dataHist2 <- table(prop_graphs_relevance_networks[[2]]$nodes$degree)
  prueba2 <- as.data.frame(dataHist2)
  plot_degree <- data.frame(num_nodes = c(prueba$Freq, prueba2$Freq), Degree = c(as.numeric(prueba$Var1),
                                                                                 as.numeric(prueba2$Var1)), type = c(rep("Early Onset", each = nrow(prueba)), rep("Late Onset", each = nrow(prueba2))))
  p1 <- ggplot(data = plot_degree, mapping = aes(x= Degree, y = num_nodes, col = type)) + 
    geom_point(alpha = 0.3, size = 2) + scale_x_log10() + scale_y_log10() + 
    ggthemes::scale_color_ptol() + ylab("Number of nodes") + 
    ggpubr::theme_pubclean() + theme(legend.title = element_blank()) #+ geom_smooth()
  
  ggsave(filename = file.path(directory, plot_directory, 'Degree_graph_relevance.png'), plot = p1, units = c("cm"),
         width = 30, height = 25, dpi = 200, limitsize = FALSE)
  
  # modules
  table_graph1 <- read.table(file.path(directory, "Results/ModuleEnrichmentAndProjection/dict_commInfomap_earlyOnset.txt"), header = TRUE, sep = "\t")
  dataHistModules1 <- table(table_graph1$infomap)
  dataModules1 <- as.data.frame(dataHistModules1)
  table_graph2 <- read.table(file.path(directory, "Results/ModuleEnrichmentAndProjection/dict_commInfomap_lateOnset.txt"), header = TRUE, sep = "\t")
  dataHistModules2 <- table(table_graph2$infomap)
  dataModules2 <- as.data.frame(dataHistModules2)
  # Prepare the data to be ploted 
  plot_modules <- data.frame(size_module = c(dataModules1$Freq, dataModules2$Freq), rank_module = c(as.numeric(dataModules1$Var1),
                                                                                                    as.numeric(dataModules2$Var1)), type = c(rep("Early Onset", each = nrow(dataModules1)), rep("Late Onset", each = nrow(dataModules2))))
  
  p2 <- ggplot(data = plot_modules, mapping = aes(x= rank_module, y = size_module, col = type)) + 
    geom_point(alpha = 0.3, size = 2)  + scale_y_log10() +
    ggthemes::scale_color_ptol() + ylab("Module size") + xlab("Rank module") +
    theme(legend.position = "top", 
          panel.background = element_rect(fill = "white", colour="black"),
          panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
          panel.grid.minor = element_blank())
  ggsave(filename = file.path(directory, plot_directory, 'modulesize_graph_relevance.png'), plot = p2, units = c("cm"),
         width = 30, height = 25, dpi = 200, limitsize = FALSE)

}