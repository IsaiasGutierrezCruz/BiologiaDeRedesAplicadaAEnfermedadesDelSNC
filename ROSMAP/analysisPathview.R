analysisPathview <- function(top2, earlyOnset_v_LateOnset.SigBOTHDIR){
  library(dplyr)
  library(gage)
  library(gageData)
  library(pathview)
 
  foldChanges <- top2$table$logFC
  names(foldChanges) <- top2$table$entrez
  
  # Get the pathways
  keggrespathways <- data.frame(id=rownames(earlyOnset_v_LateOnset.SigBOTHDIR$greater), 
                                earlyOnset_v_LateOnset.SigBOTHDIR$greater) %>% 
    tibble::as_tibble() %>% 
    filter(row_number()<=3) %>% 
    .$id %>% 
    as.character()
  
  # get the ID's
  keggresids <- substr(keggrespathways, start=1, stop=8)
  
  
  setwd("Plots")
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldChanges, pathway.id=pid, species="hsa"))
  
  setwd("../")
  
  tmp
}