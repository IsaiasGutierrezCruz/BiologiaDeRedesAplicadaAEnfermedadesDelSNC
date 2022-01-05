DataPreprocessingFunction <- function(){
  # it merges the metadata  
  source("DataPreprocessing/mergeMetadata.R")
  FullMeta <- mergeMetadata(meta_directory = 'Data/ROSMAP_biospecimen_metadata.csv', 
                            meta2_directory = 'Data/ROSMAP_clinical.csv',
                            variable_name = "individualID") 
  
  # it merges the files of the counts, change the column names and remove the useless columns
  source("DataPreprocessing/mergeCounts.R")
  RNAseqCounts <- mergeCounts(file1 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv", 
                              file2 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv", 
                              num_to_remove_beginning = 0L, num_to_remove_end = 2L, 
                              num_col_with_name = c(1, 2))
  
}