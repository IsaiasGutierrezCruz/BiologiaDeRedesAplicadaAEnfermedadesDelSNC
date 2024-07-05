create_volcano_plot <- function(plot_directory, directory){
  volcanodata <- read.csv(file.path(directory, 'Results/DGE/top2_table.csv'))
  
  volcanodata <- volcanodata[!duplicated(volcanodata$symbol), ]
  row.names(volcanodata) <- make.names(volcanodata[, "symbol"])
  
  pval_threshold <- 0.05
  logfc_threshold <- 1
  
  GED <- as.factor(ifelse(abs(volcanodata$logFC)>= logfc_threshold &
                            volcanodata$FDR < pval_threshold, 
                          "GED", "No GED"))
  
  xi <- which(GED == "GED")
  
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  g <- ggplot(data = volcanodata, aes(x=logFC, y = -log10(FDR), colour = GED))+
    geom_point(alpha = 0.4, size = 3) + theme(legend.position = "none", 
                                              panel.background = element_rect(fill = "white", colour="black"),
                                              panel.grid.major = element_line(colour = "gray", linetype = "dotted"),
                                              panel.grid.minor = element_blank()) + 
    geom_vline(xintercept = logfc_threshold) + geom_vline(xintercept = -logfc_threshold) +
    geom_hline(yintercept = -log10(pval_threshold)) + xlab("log2 fold change") + ylab("-log10 FDR") +
    geom_text_repel(aes(label = ifelse(FDR < 0.1e-6 & abs(logFC) >= logfc_threshold,
                                       row.names(volcanodata), ''))) + scale_color_manual(values = c("#6064e0", "#52854C")) 
  
  ggsave(filename = file.path(directory, plot_directory, 'volcano_plt.png'), plot = g, units = c("cm"),
         width = 30, height = 25, dpi = 200, limitsize = FALSE)
  
}

