DataPreprocessingFunction <- function(){
  # it merges the metadata  
  source("DataPreprocessing/mergeMetadata.R")
  FullMeta <- mergeMetadata(meta_directory = 'Data/ROSMAP_biospecimen_metadata.csv', 
                            meta2_directory = 'Data/ROSMAP_clinical.csv',
                            variable_name = "individualID") 
  
}