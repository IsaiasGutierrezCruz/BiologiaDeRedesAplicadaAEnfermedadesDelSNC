mergeMetadata <- function(){
  # Clinical Metadata
  Meta <- read.csv(file = 'Data/ROSMAP_clinical.csv') 
  # Biospecimen Metadata 
  Meta2 <- read.csv(file = 'Data/ROSMAP_biospecimen_metadata.csv')   
  
  # merge the two files of the metadata 
  row.names(Meta) <- Meta$individualID
  FullMeta <- cbind(Meta2, Meta[Meta2$individualID,])
  FullMeta <- FullMeta[, c(-37)]
  FullMeta
}
