library(igraph)
library(plyr)

mapflow = function(g, modulos){
  #copy graph
  d = g
  
  #detect modules
  #modulos = infomap.community(d, nb.trials = 1000)
  
  modulos = modulos
  print(length(modulos))
  #assign membership to graph
  V(d)$infomap = membership(modulos)
  
  #get community of heads and tails of crossover edges
  crossers = igraph::crossing(modulos,d)
  #if communities are not from a community igraph structure this workaround
  ##E(d)$head_comm = head_of(d, E(d))$infomap
  ##E(d)$tail_comm = tail_of(d, E(d))$infomap
  ##crossers = E(d)[tail_comm!=head_comm]
  ##
  
  x = head_of(graph = d, es = E(d)[crossers])$infomap
  y = tail_of(graph = d, es = E(d)[crossers])$infomap
  
  #make a dataframe with these communities
  z = as.data.frame(cbind(x,y))
  
  #arrange data frame so that head is always the lesser numbered community
  
  z1 = apply(X = z, MARGIN = 1, FUN = min)
  z2 = apply(X = z, MARGIN = 1, FUN = max)
  z3 = as.data.frame(cbind(z1, z2))
  
  #reduce data frame and add the count of edges between communities
  w =  plyr::ddply(z3, .(z1,z2), nrow)
  
  #get graph
  r = graph_from_data_frame(d = w, directed = FALSE)
  E(r)$weight = w$V1
  
  #number of intra community links
  intra_links = E(g)[!igraph::crossing(modulos,d)]
  #if communities are not from a community igraph structure this workaround
  ##intra_links = E(d)[tail_comm==head_comm]
  
  intras = head_of(graph = d, es = E(d)[intra_links])$infomap
  intras = table(intras)
  
  #print(intras)
  #print(V(r))
  my_intras = intras[match(names(V(r)), names(intras))]
  V(r)$intra = as.numeric(my_intras)
  
  #community size
  my_commsize   = sizes(modulos)
  
  my_commsize   = my_commsize[match(names(V(r)), names(my_commsize))]
  #print(my_commsize)
  #print(names(V(r)))
  V(r)$commsize = as.numeric(my_commsize)
  
  #add communities that have no links 
  my_missing = names(sizes(modulos))[which(!(names(sizes(modulos))%in%V(r)$name))]
  r = add_vertices(r, 
               length(my_missing), 
               attr = list(name=names(sizes(modulos))[which(!(names(sizes(modulos))%in%V(r)$name))],
                           intra=0,
                           commsize=0)
  )
  
  #return graph
  return(r)
}
