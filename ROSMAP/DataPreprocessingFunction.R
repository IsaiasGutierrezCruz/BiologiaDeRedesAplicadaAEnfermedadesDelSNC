DataPreprocessingFunction <- function(){
  # it merges the metadata  
  source("1DataPreprocessing/mergeMetadata.R")
  FullMeta <- mergeMetadata(meta_directory = 'Data/ROSMAP_biospecimen_metadata.csv', 
                            meta2_directory = 'Data/ROSMAP_clinical.csv',
                            variable_name = "individualID") 
  
  # it merges the files of the counts, change the column names and remove the useless columns
  source("1DataPreprocessing/mergeCounts.R")
  RNAseqCounts <- mergeCounts(file1 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv", 
                              file2 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv", 
                              num_to_remove_beginning = 0L, num_to_remove_end = 2L, 
                              num_col_with_name = c(1, 2))
  
  # get the metadata of the necessary samples to the study
  source("1DataPreprocessing/getSamplesToStudy.R")
  samplesToStudy <- getSamplesToStudy(RNAseqCounts = RNAseqCounts, 
                                      FullMeta = FullMeta, assay = 'rnaSeq', 
                                      cogdx = 4)
  # format the data for the analysis
  source("1DataPreprocessing/formatData.R")
  data <- formatData(RNAseqCounts = RNAseqCounts, 
                     samplesToStudy = samplesToStudy, 
                     group_names = c("EarlyOnset", "LateOnset"))
  
  list(ListSamplesToStudy = samplesToStudy, ListData = data)
}