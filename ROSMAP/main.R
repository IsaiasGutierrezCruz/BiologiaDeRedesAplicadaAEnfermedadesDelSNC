# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts(counts = "normalized")

# get the metadata of the necessary samples for the study
source("getSamplesToStudy.R")
samplesToStudy <- getSamplesToStudy(RNAseqCounts = RNAseqCounts, 
                                    FullMeta = FullMeta)

# data analysis with edgeR
source("analysisDifferenceExpression.R")
analysisDifferenceExpression(RNAseqCounts = RNAseqCounts, 
                             samplesToStudy = samplesToStudy)

