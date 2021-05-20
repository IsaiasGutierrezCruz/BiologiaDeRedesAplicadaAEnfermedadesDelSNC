getMetaSamples <- function(RNAseqCounts, FullMeta){
  samples <- names(RNAseqCounts[, 3:642])
  samples <- as.data.frame(samples)
  MetaSamples <- FullMeta[FullMeta$specimenID %in% samples$samples, ]
  MetaSamplesOnlyRNAseq <- MetaSamples[MetaSamples$assay == 'rnaSeq', ]
  
}