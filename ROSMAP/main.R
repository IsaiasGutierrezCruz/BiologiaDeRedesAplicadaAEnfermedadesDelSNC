main <- function(directory = "~/"){
  setwd(directory)
  dir.create("Plots")
  dir.create("Results")
  # -------------------------------------------------------------------------------------------
  # --------------------------- Pre-processing of the data ------------------------------------
  # -------------------------------------------------------------------------------------------
  
  source("DataPreprocessingFunction.R")
  dataPreprocessed <- DataPreprocessingFunction()
  
  samplesToStudy <- dataPreprocessed$ListSamplesToStudy
  data <- dataPreprocessed$ListData
  
  # -------------------------------------------------------------------------------------------
  # ----------- Differential gene expression analysis and pathways perturbed ------------------
  # -------------------------------------------------------------------------------------------
  
  # data analysis with edgeR
  # It normalize the count data and compare the expression of the genes in each group. 
  source("2DGEAandPathwaysPerturbed/analysisDifferenceExpression.R")
  DGEA <- analysisDifferenceExpression(countData = data$countData, 
                                       colData = data$colData, samplesToStudy = samplesToStudy, 
                                       makePlotBCV = TRUE, makePlotSmear = TRUE,
                                       output_path="Plots")
  top2 <- DGEA$top2
  countDataNormalized <- DGEA$countDataNormalized
  
  # data analysis with GAGE
  source("2DGEAandPathwaysPerturbed/analysisGAGE.R")
  earlyOnset_v_LateOnset.SigBOTHDIR <- analysisGAGE(countData = data$countData, 
                                                    samplesToStudy = samplesToStudy, 
                                                    range_kegg_pathways = c(1, 131),
                                                    same_dir = FALSE, 
                                                    output_path = "Plots/earlyOnset_v_lateOnsetGAGE")
  
  # add information about the genes to top2
  source("2DGEAandPathwaysPerturbed/addInformationTo_top2.R")
  top2$table <- addInformationTo_top2(table_to_change = top2$table,
                                      info_to_add = c("SYMBOL", "ENTREZID", "GENENAME"), 
                                      name_of_cols = c("symbol", "entrez", "name"), 
                                      delete_NA = TRUE)
  
  # representation of the pathways
  source("2DGEAandPathwaysPerturbed/analysisPathview.R")
  tmp <- analysisPathview(top2 = top2, 
                          earlyOnset_v_LateOnset.SigBOTHDIR = earlyOnset_v_LateOnset.SigBOTHDIR,
                          start_id = 1L, stop_id = 8L)
  
  # -------------------------------------------------------------------------------------------
  # --------------------------- Analysis and plots of the graphs ------------------------------
  # -------------------------------------------------------------------------------------------
  
  # analysis and plots of the graphs
  pathways_names <- c("Oxidative phosphorylation", "Cardiac muscle contraction") 
  objects_names <- c("prop_OxiPho", "prop_CardiacM") 
  etiquetas <- c("Red de la Fosforilacion Oxidativa", 
                 "Red de la Contraccion del Musculo Cardiaco")
  
  source("3AnalysisAndPlotsKeggGraphs/analysisKeggGraphs.R")
  prop_graphs <- analysisKeggGraphs(pathways_names = pathways_names, 
                                objects_names = objects_names,
                                species_name = "hsapiens", database = "kegg",
                                etiquetas = etiquetas, output_path = "Plots/")
  
  
  # Scatter plots 
  objects_names <- c("OxiPho", "CardiacM") 
  source("3AnalysisAndPlotsKeggGraphs/scatterplot_logFC_v_degree.R")
  scatterplot_logFC_v_degree(prop_graphs = prop_graphs, objects_names = objects_names, 
                             top2 = top2)
  
  # -------------------------------------------------------------------------------------------
  # --------------------------------- Relevance Networks --------------------------------------
  # -------------------------------------------------------------------------------------------
  
  endControlGroup <- length(samplesToStudy[[1]])
  endStudyGroup <- length(samplesToStudy[[2]])
  source("4RelevanceNetworks/relevanceNetworks.R")
  # get the co-expression networks 
  graphEarlyOnset <- relevanceNetworks(dataCount= countDataNormalized[, 1:endControlGroup])
  graphLateOnset <- relevanceNetworks(dataCount = countDataNormalized[, (endControlGroup + 1):(endControlGroup + endStudyGroup)])
  
  # get the networks' properties 
  source("supportFunctions/calculateGraphProperties.R")
  
  propAndgraphs_coexp <- calculateGraphProperties(networks = list(graphEarlyOnset, graphLateOnset), 
                                                names = c("graphEarlyOnset", "graphLateOnset"),
                                                calculate_communities = TRUE,
                                                keep_network_info = TRUE)
  
  # keep the properties of the co-expression networks 
  library(igraph)
  
  saveRDS(propAndgraphs_coexp[[1]], file = "Results/CoexpressionGraphs/prop_graphsCoexp.rds")
  
  # keep the co-expression networks
  saveRDS(propAndgraphs_coexp[[2]][[1]], file = "Results/CoexpressionGraphs/graphCoexp_earlyOnset.rds")
  write_graph(propAndgraphs_coexp[[2]][[1]], file = "Results/CoexpressionGraphs/graphcoexp_earlyonset.gml", format = "gml")
  
  saveRDS(propAndgraphs_coexp[[2]][[2]], file = "Results/CoexpressionGraphs/graphCoexp_lateOnset.rds")
  write_graph(propAndgraphs_coexp[[2]][[2]], file = "Results/CoexpressionGraphs/graphcoexp_lateonset.gml", format = "gml")
  
  
  
  # -------------------------------------------------------------------------------------------
  # -------------------------- Module Enrichment and Projection -------------------------------
  # -------------------------------------------------------------------------------------------
  
  # module detection
  
  g_earlyOnset = readRDS(file = "Results/CoexpressionGraphs/graphCoexp_earlyOnset.rds")
  g_lateOnset = readRDS(file = "Results/CoexpressionGraphs/graphCoexp_lateOnset.rds")
  
  source("5ModuleEnrichmentAndProjection/moduleDetection.R")
  modulesAndGraphs <- moduleDetection(graphs = list(g_earlyOnset, g_lateOnset), 
                                         names = c("earlyOnset", "lateOnset"), 
                                         output_path = "Results/ModuleEnrichmentAndProjection")
  
  modules_graphsCoexp <- modulesAndGraphs[[1]]
  graphs_with_modules <- modulesAndGraphs[[2]]
  
  # modules projection
  source("5ModuleEnrichmentAndProjection/moduleProjection.R")
  moduleProjection(graphs = graphs_with_modules, modules = modules_graphsCoexp, 
                   names = c("earlyOnset", "lateOnset"), 
                   output_path = "Results/ModuleEnrichmentAndProjection")
  
  # module enrichment
  # NOTE: expensive computational cost
  source("5ModuleEnrichmentAndProjection/moduleEnrichment.R")
  moduleEnrichment(graphs = graphs_with_modules, modules = modules_graphsCoexp, 
                   names = c("earlyOnset", "lateOnset"), 
                   output_path = "Results/ModuleEnrichmentAndProjection",
                   pvalueThreshold = 0.05)
}