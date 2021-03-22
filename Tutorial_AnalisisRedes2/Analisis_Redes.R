library(igraph)
library(dplyr)
# Leer red 
g <- read.graph(file = "data/les_miserables.graphml", format = "graphml")
g

plot(g)

# Escribir una red 
write.graph(graph = g, file = "red_mis.graphml", format = "graphml")

write.graph(graph = g, file = "red_mis.txt", format = "edgelist")

# --------- Analizar propiedades globales -------

# Componentes 
mis_componentes <- components(g)
mis_componentes

# Longitud promedio de caminos más cortos 
mis_caminos <- average.path.length(g)
mis_caminos

# Clustering coefficient global 
mi_clustering <- transitivity(g, type = "global")
mi_clustering

# Diametro 
mi_diametro <- diameter(g)
mi_diametro

# ------------- Analisis de propiedades de los nodos ------------

# Acceder a los nodos
V(g)

# Obtención completa de los nodos
mi_df_nodos <- get.data.frame(x = g, what = "vertices")
mi_df_nodos

# Añadir una nueva propiedad (en este caso una copia de label nombrada name)
V(g)$name <- V(g)$label 
mi_df_nodos <- get.data.frame(x = g, what = "vertices")
mi_df_nodos


# ---------- Calculo de algunas medidas de centralidad -------
# grado 
V(g)$grado <- degree(g)
V(g)$grado

# Interrmadiación (betweenness)
V(g)$intermediacion <- betweenness(g)
V(g)

# Obtención del coeficiente de agrupamiento local 
V(g)$agrupamiento <- transitivity(g, type = "local", isolates = "zero")

mi_df_nodos <- get.data.frame(x = g, what = "vertices")
mi_df_nodos


# ---------- Analizar propieades de los enlaces ------------
E(g)

# Propiedad edge betweenness
E(g)$intermediacion <- edge.betweenness(g)

mi_edges <- get.data.frame(x = g, what = "edges")
mi_edges




# ---------------- Clustering (comunidades) ----------------

# deteccion de las comunidades 
comm.louvain <- cluster_louvain(g)
comm.louvain

# agregar la comunidad a la que pertenece a cada nodo 
V(g)$comm.louvain <- membership(comm.louvain)


# mas algoritmos para clustering 
comm.infomap <- cluster_infomap(g)
comm.walktrp <- cluster_walktrap(g)
comm.fstgred <- cluster_fast_greedy(g)

V(g)$comm.infomap <- membership(comm.infomap)
V(g)$comm.walktrp <- membership(comm.walktrp)
V(g)$comm.fstgred <- membership(comm.fstgred)

mi_df_nodos <- get.data.frame(x = g, what = "vertices")
mi_df_nodos %>% head


# --------------- Generación de plots -------------------

# Genera una visualización conel tamaño de los nodos proporcional al grado
# escogiendo heurísticamente el "mejor acomodo"
# y con el tamaño de las etiquetas más chico 
plot(g, vertex.label.cex = 0.75, vertex.size = degree(g), layout = layout_nicely)

# Visualización con los nodos acomodados en un círculo y pone los enlaces curveados

plot(g, 
     vertex.label.cex = 0.75, 
     edge.curved = TRUE, 
     layout = layout_in_circle)




library(tidygraph)
library(ggraph)

gt <- as_tbl_graph(g)

gt %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(), show.legend = FALSE) + 
  geom_node_point(aes(colour = as.factor(comm.louvain)), size = 7) + 
  guides(colour=guide_legend(title="Louvain community")) +
  theme_graph()


