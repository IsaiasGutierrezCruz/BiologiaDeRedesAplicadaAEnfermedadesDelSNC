analysisGAGE <- function(countData, samplesToStudy, range_kegg_pathways = c(1, 131), same_dir = FALSE, output_path = "Plots/earlyOnset_v_lateOnsetGAGE"){
  # ---- Description ----
  # It identify the pathways with perturbations given the values of RNAseq
  #
  # ---- Parameters ----
  # countData: data frame
  #     Data frame with the counts from the RNASeq assay. In the columns are represented 
  #     the samples' names and in the rows are represented the genes' names 
  # SamplesToStudy: list
  #     A list that contains the names of samples in control and study group
  # range_kegg_pathways: numeric
  #     A numeric vector with the first and last pathways to be analyzed from kegg
  # same_dir: logical 
  #     A logical value to consider perturbations in same direction or in both directions
  # output_path: character 
  #     Path where the plot will be keep 
  # 
  # ---- Returns ----
  # controlGroup_v_studyGroup: list
  #     A list that contains the data frame genes perturbed in each pathway (index = 1)
  #     and un data frame with the statistic of each pathway (index = 2)
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library(gage)
  library(gageData)
  
  # change the format of the ID's 
  nombres <- data.frame(id = rownames(countData)) 
  nombres$entrez <- mapIds(org.Hs.eg.db,
                           keys=rownames(countData), 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
  
  rownames(countData) <- nombres$entrez
  
  
  # remove NA values (genes without the entrez ID)
  good <- complete.cases(rownames(countData))
  countData <- countData[good, ]
  
  # load the pathways
  data(kegg.sets.hs)
  data(sigmet.idx.hs)
  kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
  kegg.sets.hs <- kegg.sets.hs[ range_kegg_pathways[1] : range_kegg_pathways[2] ]
  
  # groups 
  control_group <- which(colnames(countData)%in%samplesToStudy[[1]])
  study_group <- which(colnames(countData)%in%samplesToStudy[[2]])
  
  # gage analysis
  controlGroup_v_studyGroupSAMEDIR <- gage(exprs = countData, 
                                           gsets = kegg.sets.hs, 
                                           ref = control_group,
                                           samp = study_group, 
                                           compare = "unpaired",
                                           same.dir = same_dir)
  
  controlGroup_v_studyGroup <- sigGeneSet(controlGroup_v_studyGroupSAMEDIR, 
                                          outname = output_path)
  
  controlGroup_v_studyGroup
}