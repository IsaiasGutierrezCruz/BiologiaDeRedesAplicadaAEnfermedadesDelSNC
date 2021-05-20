# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts()

# get the metadata of the necessary samples for the study
samples <- names(RNAseqCounts[, 3:642])
samples <- as.data.frame(samples)
MetaSamples <- FullMeta[FullMeta$specimenID %in% samples$samples, ]
MetaSamplesOnlyRNAseq <- MetaSamples[MetaSamples$assay == 'rnaSeq', ]
MetaSamplesOnlyRNAseq <- MetaSamplesOnlyRNAseq[MetaSamplesOnlyRNAseq$cogdx != 1, ]
MetaSamplesOnlyRNAseq <- MetaSamplesOnlyRNAseq[MetaSamplesOnlyRNAseq$cogdx != 6, ]


table(MetaSamplesOnlyRNAseq$age_first_ad_dx)
typeof(MetaSamplesOnlyRNAseq$age_first_ad_dx)
