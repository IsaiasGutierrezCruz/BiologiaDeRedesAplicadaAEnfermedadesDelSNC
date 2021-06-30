# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts(counts = "normalized")

# get the metadata of the necessary samples for the study
source("getSamplesToStudy.R")
samplesToStudy <- getSamplesToStudy(RNAseqCounts = RNAseqCounts, 
                                    FullMeta = FullMeta)

# format the data for the analysis
source("formatData.R")
data <- formatData(RNAseqCounts = RNAseqCounts, 
                   samplesToStudy = samplesToStudy)

# data analysis with edgeR
source("analysisDifferenceExpression.R")
top2 <- analysisDifferenceExpression(countData = data$countData, 
                             colData = data$colData)

# data analysis with GAGE
source("analysisGAGE.R")
earlyOnset_v_LateOnset.SigBOTHDIR <- analysisGAGE(countData = countData, 
                                                  samplesToStudy = samplesToStudy)

# representation of the pathways
source("analysisPathview.R")
tmp <- analysisPathview(top2 = top2, 
                        earlyOnset_v_LateOnset.SigBOTHDIR = earlyOnset_v_LateOnset.SigBOTHDIR)

# network analysis 
library(graphite)
library(igraph)

# importar las pathways
humanKegg <- pathways("hsapiens", "kegg")
OxiPho <- humanKegg[["Oxidative phosphorylation"]]
CardiacMuscle <- humanKegg[["Cardiac muscle contraction"]]

OxiPhoNEL <- pathwayGraph(OxiPho)
CardiacMuscleNEL <- pathwayGraph(CardiacMuscle)

OxiPhoIgraph <- graph_from_graphnel(OxiPhoNEL)
CardiacMuscleIgraph <- graph_from_graphnel(CardiacMuscleNEL)

# create a list to keep the properties of each pathway
prop_OxiPho <- list()
prop_CardiacM <- list()

# -------- Analizar propiedades globales --------------

# Componentes 
prop_OxiPho <- append(prop_OxiPho, list(components = components(OxiPhoIgraph)))
prop_CardiacM <- append(prop_CardiacM, list(components = components(CardiacMuscleIgraph)))

# Longitud promedio de caminos mas cortos 
prop_OxiPho <- append(prop_OxiPho, list(average.path.length = average.path.length(OxiPhoIgraph)))
prop_CardiacM <- append(prop_CardiacM, list(average.path.length = average.path.length(CardiacMuscleIgraph)))

# Clusterring coefficient global 
prop_OxiPho <- append(prop_OxiPho, list(transitivity = transitivity(OxiPhoIgraph, type = "global")))
prop_CardiacM <- append(prop_CardiacM, list(transitivity = transitivity(CardiacMuscleIgraph, type = "global")))

# Diametro 
prop_OxiPho <- append(prop_OxiPho, list(diameter = diameter(OxiPhoIgraph)))
prop_CardiacM <- append(prop_CardiacM, list(diameter = diameter(CardiacMuscleIgraph)))

# ------------- Analisis de propiedades de los nodos ------------

# degree
V(OxiPhoIgraph)$degree <- degree(OxiPhoIgraph)
V(CardiacMuscleIgraph)$degree <- degree(CardiacMuscleIgraph)

# Intermediacion (betweenness)
V(OxiPhoIgraph)$betweenness <- betweenness(OxiPhoIgraph)
V(CardiacMuscleIgraph)$betweenness <- betweenness(CardiacMuscleIgraph)

# Obtención del coeficiente de agrupamiento local
V(OxiPhoIgraph)$transitivity <- transitivity(OxiPhoIgraph, type = "local", isolates = "zero")
V(CardiacMuscleIgraph)$transitivity <- transitivity(CardiacMuscleIgraph, type = "local", isolates = "zero")

# Obtención completa de los nodos
prop_OxiPho <- append(prop_OxiPho, list(nodes = get.data.frame(x = OxiPhoIgraph, what = "vertices")))
prop_CardiacM <- append(prop_CardiacM, list(nodes = get.data.frame(x = CardiacMuscleIgraph, what = "vertices")))


# ---------- Analizar propieades de los enlaces ------------
E(OxiPhoIgraph)$intermediacion <- edge.betweenness(OxiPhoIgraph)
E(CardiacMuscleIgraph)$intermediacion <- edge.betweenness(CardiacMuscleIgraph)

prop_OxiPho <- append(prop_OxiPho, list(edges = get.data.frame(x = OxiPhoIgraph, what = "edges")))
prop_CardiacM <- append(prop_CardiacM, list(edges = get.data.frame(x = CardiacMuscleIgraph, what = "edges")))



# -------------------- Scatterplot -------------------------------------

library("AnnotationDbi")
library("org.Hs.eg.db")

rowWithLength18 <- which((nchar(rownames(top2)) == 18) == TRUE)
rowWithLength17 <- which((nchar(rownames(top2)) == 17) == TRUE)

rownames(top2$table)[rowWithLength17] <- substr(rownames(top2$table)[rowWithLength17], 1, 
                                                nchar(rownames(top2$table)[rowWithLength17])-2)

rownames(top2$table)[rowWithLength18] <- substr(rownames(top2$table)[rowWithLength18], 1, 
                                                nchar(rownames(top2$table)[rowWithLength18])-3)


top2$table$symbol <- mapIds(org.Hs.eg.db, 
                            keys = row.names(top2),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")

top2$table$entrez = mapIds(org.Hs.eg.db,
                           keys=row.names(top2), 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")

top2$table$name =   mapIds(org.Hs.eg.db,
                           keys=row.names(top2), 
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")

# remove NA values 
good <- complete.cases(top2$table$entrez)
top2$table <- top2$table[good, ]


prop_OxiPho$nodes

# Remove description of the row names 
# OxiPho
rownames(prop_OxiPho$nodes) <- substr(rownames(prop_OxiPho$nodes), 10, 
                                                       nchar(rownames(prop_OxiPho$nodes)))
rownames(prop_CardiacM$nodes) <- substr(rownames(prop_CardiacM$nodes), 10, 
                                                       nchar(rownames(prop_CardiacM$nodes)))


# hacer el data frame
genesOxiPho <- prop_OxiPho$nodes[which(rownames(prop_OxiPho$nodes) %in% top2$table$entrez), ]

OxiPho_scatterplot <- data.frame(degree = genesOxiPho[, "degree"])
rownames(OxiPho_scatterplot) <- rownames(genesOxiPho)

genesOxiPho <- which(top2$table$entrez %in% rownames(genesOxiPho))
OxiPho_logFC <- top2$table[genesOxiPho, ]
OxiPho_scatterplot <- cbind(OxiPho_scatterplot, data.frame(logFC = OxiPho_logFC$logFC))
rownames(OxiPho_scatterplot) <- OxiPho_logFC$name
OxiPho_scatterplot <- mutate(OxiPho_scatterplot, names = OxiPho_logFC$name)

# CardiacM
# hacer el data frame
genesCardiacM <- prop_CardiacM$nodes[which(rownames(prop_CardiacM$nodes) %in% top2$table$entrez), ]

CardiacM_scatterplot <- data.frame(degree = genesCardiacM[, "degree"])
rownames(CardiacM_scatterplot) <- rownames(genesCardiacM)

genesCardiacM <- which(top2$table$entrez %in% rownames(genesCardiacM))
CardiacM_logFC <- top2$table[genesCardiacM, ]
CardiacM_scatterplot <- cbind(CardiacM_scatterplot, data.frame(logFC = CardiacM_logFC$logFC))
rownames(CardiacM_scatterplot) <- CardiacM_logFC$name
CardiacM_scatterplot <- mutate(CardiacM_scatterplot, names = CardiacM_logFC$name)

library(ggplot2)
library(dplyr)
library(ggrepel)

F1 <- ggplot(OxiPho_scatterplot, aes(x = degree, y = logFC, color = degree)) + 
  geom_point(size = 3) + 
  geom_label_repel(
    data = OxiPho_scatterplot %>% filter(names == "ATPase H+ transporting V0 subunit a1"), 
    aes(label = names),
    nudge_x = -20, 
    nudge_y = -0.08,
    color = "black", 
    fill = "lightskyblue1", 
    force = 600
  ) + ggtitle("Vía de la Fosforilación Oxidativa") +labs(x = "Degree") + 
  ggplot2::theme_minimal()
  


F2 <- ggplot(CardiacM_scatterplot, aes(x = degree, y = logFC, color = degree)) + 
  geom_point(size = 4) + 
  geom_label_repel(
    data = CardiacM_scatterplot %>% filter(logFC > 0.5 | logFC < -0.25), 
    aes(label = names),
    nudge_x = 0, 
    nudge_y = -0.08,
    color = "black", 
    fill = "lightskyblue1", 
    force = 100
  ) + ggtitle("Vía de la Contracción del Musculo Cardiaco") +labs(x = "Degree") + 
  ggplot2::theme_minimal()


F1 + theme(plot.title = element_text(face = "bold", hjust = 0.5))
F2 + theme(plot.title = element_text(face = "bold", hjust = 0.5))

# graficar networks
plot(OxiPhoIgraph, vertex.label.cex = 0.75, vertex.size = degree(OxiPhoIgraph), layout = layout_nicely)

plot(OxiPhoIgraph, 
     vertex.label.cex = 0.75, 
     edge.curved = TRUE, 
     layout = layout_in_circle)

library(tidygraph)
library(ggraph)

gt <- as_tbl_graph(OxiPhoIgraph)

gt %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(), show.legend = FALSE) + 
  geom_node_point(aes(size = as.factor(degree)), colour = "lightskyblue1") +
  theme_graph()


plot(CardiacMuscleIgraph, vertex.label.cex = 0.75, vertex.size = degree(OxiPhoIgraph), layout = layout_nicely)

plot(CardiacMuscleIgraph, 
     vertex.label.cex = 0.75, 
     edge.curved = TRUE, 
     layout = layout_in_circle)

gt <- as_tbl_graph(CardiacMuscleIgraph)

gt %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(), show.legend = FALSE) + 
  geom_node_point(aes(filter = degree == 4 | degree == 5, size = as.factor(degree)), colour = "lightskyblue1") +
  theme_graph()


