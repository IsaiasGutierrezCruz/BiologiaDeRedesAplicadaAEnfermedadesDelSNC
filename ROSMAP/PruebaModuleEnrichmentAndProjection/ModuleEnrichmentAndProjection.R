library(readr)

load("E:/DataROSMAPNetwork/ModuleEnrichmentAndProjection/GOs_and_pathways.RData")

modules_cases <- readRDS("resultsWithOnlySymbl/modules_cases.RDS")
modules_cntrl <- readRDS("resultsWithOnlySymbl/modules_cntrl.RDS")

flow_cntrl <- readRDS("E:/DataROSMAPNetwork/DataOnlyWithSymbol/flow_cntrl.RDS")
flow_cases <- readRDS("E:/DataROSMAPNetwork/DataOnlyWithSymbol/flow_cases.RDS")


source("HTSanalyzeR.R")

# cases 
saveRDS(enrichment_list_cntrl, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/enrichment_list_cntrl_earlyOnset.RDS")

# enrichment_list_cases
enrichment_df_cases <- read_table("E:/DataROSMAPNetwork/DataOnlyWithSymbol/enrichment_df_cases.txt", col_names = TRUE)

saveRDS(bp_cases, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/bp_cases.RDS")
readRDS("E:/DataROSMAPNetwork/DataOnlyWithSymbol/bp_cases.RDS")


# cntrl 
saveRDS(enrichment_list_cases, file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/enrichment_list_cases_lateOnset.RDS")

