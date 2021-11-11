library(minet)
library(Rgraphviz)
library(igraph)
library("AnnotationDbi")
library("org.Hs.eg.db")

# early Onset 
y$counts[, 1:47]

saveRDS(y$counts[, 1:47], file = "earlyOnset_counts.rds")
earlyOnset <- readRDS(file = "E:/DataROSMAPNetwork/Data/earlyOnset_counts.rds")

# add key symbol
nombresEarlyOnset <- mapIds(org.Hs.eg.db, 
                            keys = row.names(earlyOnset),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
nombresEarlyOnset <- as.vector(nombresEarlyOnset)

row.names(earlyOnset) <- nombresEarlyOnset

good <- complete.cases(row.names(earlyOnset))
earlyOnset <- earlyOnset[good, ]

earlyOnset <- earlyOnset[unique(row.names(earlyOnset)), ]


earlyOnset <- as.data.frame(t(earlyOnset))

mim <- build.mim(dataset = earlyOnset, estimator = "spearman")


saveRDS(mim, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/mutual_information_earlyOnset.rds")
mim <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/mutual_information_earlyOnset.rds")

# relevant network 

umbral <- quantile(mim, 0.99, na.rm = TRUE)
mim2 <- ifelse(mim < umbral, 0, 1)
mim2[is.na(mim2)] <- 0
mimgraph <- graph_from_adjacency_matrix(mim2, mode = "undirected")

saveRDS(mimgraph, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_earlyOnset.rds")
mimgraph <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_earlyOnset.rds")

plot(mimgraph, vertex.label.cex = 0.5, layout = layout_nicely)



#late onset 
y$counts[, 48:221]

saveRDS(y$counts[, 48:221], file = "lateOnset_counts.rds")
lateOnset <- readRDS(file = "E:/DataROSMAPNetwork/Data/lateOnset_counts.rds")

# add key symbol
nombresLateOnset <- mapIds(org.Hs.eg.db, 
                            keys = row.names(lateOnset),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
nombresLateOnset <- as.vector(nombresLateOnset)

row.names(lateOnset) <- nombresLateOnset

good <- complete.cases(row.names(lateOnset))
lateOnset <- lateOnset[good, ]

lateOnset <- lateOnset[unique(row.names(lateOnset)), ]



lateOnset <- as.data.frame(t(lateOnset))

mimLO <- build.mim(dataset = lateOnset, estimator = "spearman")

saveRDS(mimLO, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/mutual_information_lateOnset.rds")
mimLO <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/mutual_information_lateOnset.rds")

# relevance network
umbral_lo <- quantile(mimLO, 0.99, na.rm = TRUE)
mimLO2 <- ifelse(mimLO < umbral_lo, 0, 1)
mimLO2[is.na(mimLO2)] <- 0
mimgraph_lo <- graph_from_adjacency_matrix(mimLO2, mode = "undirected")

saveRDS(mimgraph_lo, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_lateOnset.rds")
mimgraph_lo <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_lateOnset.rds")

plot(mimgraph_lo, vertex.label.cex = 0.5, layout = layout_nicely)

# ---- calculate the description measurements ----

library(igraph)
library(dplyr)

mimgraph <- readRDS(file = "E:/DataROSMAPNetwork/Data/network_earlyOnset.rds")
mimgraph_lo <- readRDS(file = "E:/DataROSMAPNetwork/Data/network_lateOnset.rds")

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


# ---- lateOnset ----
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


saveRDS(prop_graphs, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/prop_graphs.rds")
prop_graphs <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/prop_graphs.rds")

# ---- keep the graphs with info

saveRDS(mimgraph, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_earlyOnset_with_info.rds")
mimgraph <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_earlyOnset_with_info.rds")
saveRDS(mimgraph_lo, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_lateOnset_with_info.rds")
mimgraph_lo <- readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_lateOnset_with_info.rds")


# ---- plot ----

library(tidygraph)
library(ggraph)

gt <- as_tbl_graph(mimgraph)
saveRDS(gt, file = "E:/DataROSMAPNetwork/Data/earlyOnsetTryToGraph/tbl_graph_earlyOnset.rds")

plot <- gt %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(), show.legend = FALSE) + 
  geom_node_point(aes(colour = as.factor(comm.louvain)), size = 7) + 
  guides(colour=guide_legend(title="Louvain community")) +
  theme_graph()

saveRDS(plot, file = "E:/DataROSMAPNetwork/Data/earlyOnsetTryToGraph/plot_graph_earlyOnset.rds")

ggsave(filename = "Network_earlyOnset.png", plot = plot, units = c("cm"),
       width = 30, height = 20, dpi = 200, limitsize = FALSE)




# ---- get the formats gml  ----
write_graph(mimgraph, file = "E:/DataROSMAPNetwork/Data/pruebas_graphs/earlyOnset_graph.gml", format = "gml")

write_graph(mimgraph_lo, file = "E:/DataROSMAPNetwork/Data/pruebas_graphs/lateOnset_graph.gml", format = "gml")
