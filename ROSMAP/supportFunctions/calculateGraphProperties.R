calculateGraphProperties <- function(networks, names, calculate_communities = FALSE, keep_network_info = FALSE){
  # ---- Description ----
  # It calculate graph properties of several network and it keep it in a list
  #
  # ---- Parameters ----
  # networks: list
  #     A list of networks
  # names: character
  #     Character vector with the names of each network
  # calculate_communities: logical 
  #     Logical value to determine calculate communities or not 
  # keep_network_info: logical 
  #     Logical value to determine keep the network or not 
  # 
  # 
  # ---- Returns ----
  # result: list
  #     A list with the properties of each network (index = 1) and a list of
  #     the network with the information (index = 2). 
  #     IMPORTANT: the networks in the list should be obtained with two indexes 
  #     example: [[1]][[1]]
  
  if (!length(networks) == length(names)){
    stop("The number of netwroks and names are different")
  }
  
  prop_graphs <- list()
  
  for (i in seq_along(networks)){
    prop_graphs[[ names[i] ]] <- list()
    # components
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(components = components( networks[[i]] )))
    # average path length
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(average.path.length = average.path.length( networks[[i]] )))
    # Clustering coefficient global 
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(transitivity = transitivity( networks[[i]] )))
    # Diameter 
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(diameter = diameter( networks[[i]] )))
    
    # nodes properties analysis
    # degree
    V( networks[[i]] )$degree <- degree( networks[[i]] )
    # betweenness
    V( networks[[i]] )$betweenness <- betweenness( networks[[i]] )
    # get the local clustering coefficient 
    V( networks[[i]] )$transitivity <- transitivity( networks[[i]] , type = "local", isolates = "zero")
    
    if (calculate_communities){
      # comunidades 
      comm.louvain <- cluster_louvain( networks[[i]] )
      V( networks[[i]] )$comm.louvain <- membership(comm.louvain)
      
      modules_cases = infomap.community(graph =  networks[[i]], nb.trials = 1)
      #assign module membership to vertices
      V( networks[[i]] )$infomap = membership(modules_cases)
    }
    
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(nodes = get.data.frame(x =  networks[[i]], what = "vertices")))
    
    
    # edges 
    E( networks[[i]] )$betweenness <- edge.betweenness( networks[[i]] )
    
    prop_graphs[[ names[i] ]] <- append(prop_graphs[[ names[i] ]], list(edges = get.data.frame(x =  networks[[i]], what = "edges")))
  }
  
  # determine keep networks and properties or only properties
  if(keep_network_info){
    result <- list(prop_graphs = prop_graphs, graph_with_properties = networks)
  } else if (keep_network_info == FALSE){
    result <- prop_graphs
  }
  
  result
}


