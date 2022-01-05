mergeMetadata <- function(meta_directory = 'Data/ROSMAP_biospecimen_metadata.csv', 
                          meta2_directory = 'Data/ROSMAP_clinical.csv',
                          variable_name = "individualID"){
  # ---- Description ----
  # It merge two files of metadata. Both of them have the same individual ID but one of them have multiples samples for each patient
  #
  # ---- Parameters ----
  # meta_directory: character 
  #     Directory of the Metadata file 1. This file has different samples from each patient
  # meta2_directory: character
  #     Directory of the Metadata file 2. This file has a unique value for each patient 
  # variable_name: character
  #     Name of the variable shared by both data frames 
  # 
  # ---- Returns ----
  # FullMeta: data.frame
  #     Data frame with metadata from each sample 
  
  
  # Biospecimen Metadata 
  Meta <- read.csv(file = meta_directory)   
  # Clinical Metadata
  Meta2 <- read.csv(file = meta2_directory) 
  
  # merge the two files of the metadata 
  FullMeta <- merge(Meta, Meta2, by = c(variable_name, variable_name),  all.x = TRUE)
  FullMeta
}

