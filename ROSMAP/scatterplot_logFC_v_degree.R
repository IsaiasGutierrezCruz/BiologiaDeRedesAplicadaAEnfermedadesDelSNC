scatterplot_logFC_v_degree <- function(prop_graphs, objects_names, top2){
  scatterplots_dataframes <- list()
  
  for (i in seq_along(prop_graphs)){
    properties <- prop_graphs[[i]]
    # Remove description of the row names 
    rownames(properties$nodes) <- substr(rownames(properties$nodes), 10, 
                                          nchar(rownames(properties$nodes)))
    # make the data frame
    genes <- properties$nodes[which(rownames(properties$nodes) %in% top2$table$entrez), ]
    # add degree to the data frame
    scatterplot_dataframe <- data.frame(degree = genes[, "degree"])
    rownames(scatterplot_dataframe) <- rownames(genes)
    
    # add logFC to the data frame
    table_logFC <- top2$table[which(top2$table$entrez %in% rownames(genes)), ]
    scatterplot_dataframe <- cbind(scatterplot_dataframe, data.frame(logFC = table_logFC$logFC))
    
    scatterplot_dataframe <- mutate(scatterplot_dataframe, names = table_logFC$name)
    
    
    scatterplots_dataframes[[i]] <- scatterplot_dataframe
  }
  names(scatterplots_dataframes) <- objects_names
  
  
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggpubr)
  
  F1 <- ggplot(scatterplots_dataframes$OxiPho, aes(x = degree, y = logFC, color = degree)) + 
    geom_point(size = 3) + 
    geom_label_repel(
      data = scatterplots_dataframes$OxiPho %>% filter(names == "ATPase H+ transporting V0 subunit a1"), 
      aes(label = names),
      nudge_x = -20, 
      nudge_y = -0.08,
      color = "black", 
      fill = "lightskyblue1", 
      force = 600
    ) + ggtitle("Via de la Fosforilacion Oxidativa") +labs(x = "Degree") + 
    ggplot2::theme_minimal()
  
  F2 <- ggplot(scatterplots_dataframes$CardiacM, aes(x = degree, y = logFC, color = degree)) + 
    geom_point(size = 4) + 
    geom_label_repel(
      data = scatterplots_dataframes$CardiacM %>% filter(logFC > 0.5 | logFC < -0.25), 
      aes(label = names),
      nudge_x = 0, 
      nudge_y = -0.08,
      color = "black", 
      fill = "lightskyblue1", 
      force = 100
    ) + ggtitle("Via de la Contraccion del Musculo Cardiaco") +labs(x = "Degree") + 
    ggplot2::theme_minimal()
  
  F1 <- F1 + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  F2 <- F2 + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  save <- ggarrange(F1, F2, ncol = 1, nrow = 2)
  
  ggsave(filename = "Plots/Scatterplots_logFC_v_degree.png", plot = save, 
         units = c("cm"), width = 20, height = 28, dpi = 200, limitsize = FALSE,
         bg = "white")
}