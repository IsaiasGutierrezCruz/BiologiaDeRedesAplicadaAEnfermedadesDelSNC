getSamplesToStudy <- function(RNAseqCounts, FullMeta){
  # get the metadata of the samples of interest
  columns <- ncol(RNAseqCounts)
  samples <- names(RNAseqCounts[, 3:columns])
  samples <- as.data.frame(samples)
  MetaSamples <- FullMeta[FullMeta$specimenID %in% samples$samples, ]
  
  # filter with the conditions of interest
  MetaSamplesOnlyRNAseq <- MetaSamples[MetaSamples$assay == 'rnaSeq', ]
  MetaSamplesOnlyRNAseq <- MetaSamplesOnlyRNAseq[MetaSamplesOnlyRNAseq$cogdx == 4, ]
  
  LessThan70Years <- MetaSamplesOnlyRNAseq[MetaSamplesOnlyRNAseq$age_first_ad_dx == "", ]
  LessThan70Years <- LessThan70Years[1:47, ]
  GreaterThan70Years <- MetaSamplesOnlyRNAseq[MetaSamplesOnlyRNAseq$age_first_ad_dx != "", ]
  GreaterThan70Years <- GreaterThan70Years[1:174, ]
  
  # get the samples of interest 
  samplesLessThan70Years <- LessThan70Years$specimenID
  samplesGreaterThan70Years <- GreaterThan70Years$specimenID
  dividedSamples <- list(LessThan70Years = samplesLessThan70Years, 
                         GreaterThan70Years = samplesGreaterThan70Years)
  dividedSamples
}