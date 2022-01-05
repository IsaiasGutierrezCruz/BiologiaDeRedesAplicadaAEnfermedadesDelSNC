main <- function(directory = "~/"){
  setwd(directory)
  dir.create("Plots")
  
  # ----------- Pre-processing of the data -----------------
  source("DataPreprocessingFunction.R")
  dataPreprocessed <- DataPreprocessingFunction()
  
  samplesToStudy <- dataPreprocessed$ListSamplesToStudy
  data <- dataPreprocessed$ListData
  
  
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
                               colData = data$colData, samplesToStudy = samplesToStudy)
  
  # data analysis with GAGE
  source("analysisGAGE.R")
  earlyOnset_v_LateOnset.SigBOTHDIR <- analysisGAGE(countData = data$countData, 
                                                    samplesToStudy = samplesToStudy)
  
  # add information about the genes to top2
  source("addInformationTo_top2.R")
  top2 <- addInformationTo_top2(top2 = top2)
  
  # representation of the pathways
  source("analysisPathview.R")
  tmp <- analysisPathview(top2 = top2, 
                          earlyOnset_v_LateOnset.SigBOTHDIR = earlyOnset_v_LateOnset.SigBOTHDIR)
  
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
  
  
  # -------------------- Scatter plots -------------------------------------
  objects_names <- c("OxiPho", "CardiacM") 
  source("scatterplot_logFC_v_degree.R")
  scatterplot_logFC_v_degree(prop_graphs = prop_graphs, objects_names = objects_names, 
                             top2 = top2)
}
