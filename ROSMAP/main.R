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

# format the data for the analysis
source("formatData.R")
data <- formatData(RNAseqCounts = RNAseqCounts, 
                   samplesToStudy = samplesToStudy)

countData <- data$countData
colData <- data$colData


# data analysis with edgeR
source("analysisDifferenceExpression.R")
top2 <- analysisDifferenceExpression(countData = countData, 
                             colData = colData)

# data analysis with GAGE
source("analysisGAGE.R")
earlyOnset_v_LateOnset.SigBOTHDIR <- analysisGAGE(countData = countData, 
                                                  samplesToStudy = samplesToStudy)

# representation of the pathways
source("analysisPathview.R")
tmp <- analysisPathview(top2 = top2, 
                        earlyOnset_v_LateOnset.SigBOTHDIR = earlyOnset_v_LateOnset.SigBOTHDIR)
