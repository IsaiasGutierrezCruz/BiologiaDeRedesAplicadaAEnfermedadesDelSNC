main <- function(directory = "~/"){
  setwd(directory)
  dir.create("Plots")
  
  # ----------- Pre-processing of the data -----------------
  
  source("DataPreprocessingFunction.R")
  dataPreprocessed <- DataPreprocessingFunction()
  
  samplesToStudy <- dataPreprocessed$ListSamplesToStudy
  data <- dataPreprocessed$ListData
  
  
  # ----------- Differential gene expression analysis and pathways perturbed ------------------
  
  # data analysis with edgeR
  # It normalize the count data and compare the expression of the genes in each group. 
  source("DGEAandPathwaysPerturbed/analysisDifferenceExpression.R")
  DGEA <- analysisDifferenceExpression(countData = data$countData, 
                                       colData = data$colData, samplesToStudy = samplesToStudy, 
                                       makePlotBCV = FALSE, makePlotSmear = TRUE)
  top2 <- DGEA$top2
  countDataNormalized <- DGEA$countDataNormalized
  
  # data analysis with GAGE
  source("DGEAandPathwaysPerturbed/analysisGAGE.R")
  earlyOnset_v_LateOnset.SigBOTHDIR <- analysisGAGE(countData = data$countData, 
                                                    samplesToStudy = samplesToStudy, 
                                                    range_kegg_pathways = c(1, 131),
                                                    same_dir = FALSE, 
                                                    output_path = "Plots/earlyOnset_v_lateOnsetGAGE")
  
  # add information about the genes to top2
  source("DGEAandPathwaysPerturbed/addInformationTo_top2.R")
  top2$table <- addInformationTo_top2(table_to_change = top2$table,
                                      info_to_add = c("SYMBOL", "ENTREZID", "GENENAME"), 
                                      name_of_cols = c("symbol", "entrez", "name"), 
                                      delete_NA = TRUE)
  
  # representation of the pathways
  source("DGEAandPathwaysPerturbed/analysisPathview.R")
  tmp <- analysisPathview(top2 = top2, 
                          earlyOnset_v_LateOnset.SigBOTHDIR = earlyOnset_v_LateOnset.SigBOTHDIR,
                          start_id = 1L, stop_id = 8L)
  
  # ------------------ Analysis and plots of the graphs -----------------------------
  
  # analysis and plots of the graphs
  pathways_names <- c("Oxidative phosphorylation", "Cardiac muscle contraction") 
  objects_names <- c("prop_OxiPho", "prop_CardiacM") 
  etiquetas <- c("Red de la Fosforilacion Oxidativa", 
                 "Red de la Contraccion del Musculo Cardiaco")
  
  source("analysisGraphs.R")
  prop_graphs <- analysisGraphs(pathways_names = pathways_names, 
                                objects_names = objects_names,
                                species_name = "hsapiens", database = "kegg",
                                etiquetas = etiquetas)
  
  
  # Scatter plots 
  objects_names <- c("OxiPho", "CardiacM") 
  source("scatterplot_logFC_v_degree.R")
  scatterplot_logFC_v_degree(prop_graphs = prop_graphs, objects_names = objects_names, 
                             top2 = top2)
}
