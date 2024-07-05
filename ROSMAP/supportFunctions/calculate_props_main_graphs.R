calculate_props_main_graphs <- function(directory){
  library(igraph)
  source(file.path(directory, 'supportFunctions/calculate_prop_graphs_to_plot.R'))
  g_flow_early_onset <- read.graph(file = file.path(directory, 'Results/ModuleEnrichmentAndProjection/flow_earlyOnset.gml'), format = "gml")
  g_flow_late_onset <- read.graph(file = file.path(directory, 'Results/ModuleEnrichmentAndProjection/flow_lateOnset.gml'), format = "gml")
  prop_graphs_flow <- graphProps(g_flow_early_onset, g_flow_late_onset)
  
  saveRDS(prop_graphs_flow, file = file.path(directory, "Results/ModuleEnrichmentAndProjection/prop_graphs_flow.RDS"))
  
  g_bipartite_early_onset <- read.graph(file = file.path(directory, 'Results/ModuleEnrichmentAndProjection/bipartite_earlyOnset.gml'), format = "gml")
  g_bipartite_late_onset <-read.graph(file = file.path(directory, 'Results/ModuleEnrichmentAndProjection/bipartite_lateOnset.gml'), format = "gml") 
  prop_graphs_bipartite <- graphProps(g_bipartite_early_onset, g_bipartite_late_onset)
  
  saveRDS(prop_graphs_bipartite, file = file.path(directory, "Results/ModuleEnrichmentAndProjection/prop_graphs_bipartite.RDS"))
}