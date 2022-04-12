getSamplesToStudy <- function(RNAseqCounts, FullMeta, assay = 'rnaSeq', cogdx = 4, age_cutoff=70){
  # ---- Description ----
  # Selection of the metadata of samples used in the RNAseq assay for later apply several 
  # filters to select specific variables of interest 
  #
  # ---- Parameters ----
  # RNASeqCounts: data frame
  #     Data frame with the counts from the RNASeq assay. In the columns are represented 
  #     the samples' names and in the rows are represented the genes' names 
  # FullMeta: data frame
  #     Data frame with the metadata of the samples
  # assay: character
  #     Assay of interest 
  # cogdx: integer
  #     Number of the category of cogdx
  # age_cutoff: integer 
  #     Integer value with the cutoff to use
  # 
  # ---- Returns ----
  # SamplesToStudy: list
  #     A list that contains the names of samples in control (index =1) and study group (index = 2)
  
  
  # get the metadata of the samples of interest
  samples <- as.data.frame(names(RNAseqCounts))
  MetaSamples <- FullMeta[FullMeta$specimenID %in% samples[[1]], ]
  
  # filter with the conditions of interest
  metaSamplesRNAseq <- MetaSamples[MetaSamples$assay == assay, ]
  metaSamplesRNAseq <- metaSamplesRNAseq[metaSamplesRNAseq$cogdx == cogdx, ]
  
  # remove NA values 
  metaSamplesRNAseq <- metaSamplesRNAseq[!is.na(metaSamplesRNAseq$age_first_ad_dx), ]
  metaSamplesRNAseq[metaSamplesRNAseq$age_first_ad_dx == "90+", ]$age_first_ad_dx <- "90"
  
  # divide the data metadata in groups
  LessThan70Years <- metaSamplesRNAseq[metaSamplesRNAseq$age_first_ad_dx < age_cutoff, ]
  GreaterThan70Years <- metaSamplesRNAseq[metaSamplesRNAseq$age_first_ad_dx > age_cutoff, ]
  
  # get the samples of interest 
  samplesLessThan70Years <- LessThan70Years$specimenID
  samplesGreaterThan70Years <- GreaterThan70Years$specimenID
  dividedSamples <- list(LessThan70Years = samplesLessThan70Years, 
                         GreaterThan70Years = samplesGreaterThan70Years)
  dividedSamples
}