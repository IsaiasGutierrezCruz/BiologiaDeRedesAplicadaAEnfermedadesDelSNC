analysisPathview <- function(top2, earlyOnset_v_LateOnset.SigBOTHDIR, start_id = 1L, stop_id = 8L){
  # ---- Description ----
  # Make representation of the pathways perturbed 
  #
  # ---- Parameters ----
  # top2: list
  #     A list that contains the data of the gene expression comparison and the 
  #     variables adjust.method, comparison and test. The columns logFc and entrez
  #     in the data frame is needed 
  # earlyOnset_v_LateOnset.SigBOTHDIR: list
  #     A list that contains the data frame genes perturbed in each pathway (index = 1)
  #     and un data frame with the statistic of each pathway (index = 2)
  # start_id: integer 
  #     Integer to locate the beginning of the pathway' id
  # stop_id: integer 
  #     Integer to locate the end of the pathway' id
  # 
  # ---- Returns ----
  # tmp: list
  #     List with metadata of the representations
  
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