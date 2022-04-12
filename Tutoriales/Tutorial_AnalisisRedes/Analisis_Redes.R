library(igraph)

# -------- Leer redes ----------
# graphml
g <- igraph::read.graph(file = "data/ZacharyKarateNetwork.graphml", format = "graphml")
g

# edge list
g2 <- igraph::read.graph(file = "data/ZacharyKarateNetwork_edgelist.txt",
                         format = "edgelist", directed = F)
g2

# comprobar si se trata del mismo grafo 
igraph::is_isomorphic_to(graph1 = g, graph2 = g2)

# sin embargo no son identicos porque g tiene diferentes atributos de nodo
igraph::identical_graphs(g, g2)

# leer grafo desde un frame
df <- read.table(file = "data/ZacharyKarateNetwork_edgelist.txt")
print(head(df))

g3 <- graph_from_data_frame(d = df, directed = F)
g3

# es posible convertir un grafo a una matriz de adyacencia 
mat <- get.adjacency(graph = g, sparse = F)
mat[1:5, 1:5]

g4 <- graph_from_adjacency_matrix(adjmatrix = mat, mode = "undirected")
g4

# escribir una red 
write.graph(graph = g, file = "ejemplo_escribir.graphml", format = "graphml")


# ------ Analizar redes ------------
set.seed(725)
plot(g)

# Medir caminos 
matriz_de_caminos <- shortest.paths(g)
matriz_de_caminos[1:5, 1:5]

avg_shortest_path <- average.path.length(g)
avg_shortest_path

my_diameter <- diameter(g)
my_diameter

# Revisar los componentes conectados 
components_g <- components(g)
components_g$membership

components_g$csize

components_g$no

# Determinar las centralidades (se accede a los nodos)
V(g)

# Determinar el Grado 
g <- set.vertex.attribute(graph = g, name = "degree", value = degree(g))
# ó: 
# V(g)$degree_alt <- degree(g)

head(get.data.frame(g, what = "vertices"))

# Intermediación 
V(g)$betweenness_centrality <- betweenness(g)
head(get.data.frame(g, what = "vertices"))

E(g)$edge_betweenness <- edge.betweenness(g)
head(get.data.frame(g, what = "edges"))

# Coeficiente de agrupamiento/Transitividad 
V(g)$cc <- transitivity(graph = g, type = "local", isolates = "zero")
cc_global <- transitivity(graph = g, type = "global")
head(get.data.frame(g, what = "vertices"))
