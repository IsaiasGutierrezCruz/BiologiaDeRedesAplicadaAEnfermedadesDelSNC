# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts()

# get the metadata of the necessary samples for the study
source("getSamplesToStudy.R")

