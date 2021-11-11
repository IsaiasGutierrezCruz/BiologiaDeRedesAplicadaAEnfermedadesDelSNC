#1) Module detection

#libraries
library(igraph)
library(data.table)
library(tidyverse)
#read graph

g_cntrl = readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_earlyOnset_with_info.rds")
g_cases = readRDS(file = "E:/DataROSMAPNetwork/DataOnlyWithSymbol/network_lateOnset_with_info.rds")

#cases
#E(g_cases)$weight = E(g_cases)$V2

#detect modules
modules_cases = infomap.community(graph = g_cases, nb.trials = 1)
#assign module membership to vertices
V(g_cases)$infomap = membership(modules_cases)

write.table(x = get.data.frame(x = g_cases, "vertices"), 
            file = "resultsWithOnlySymbl/dict_commInfomap_cases.txt", 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

#cntrl
#E(g_cntrl)$weight = E(g_cntrl)$V2

#detect modules
modules_cntrl = infomap.community(graph = g_cntrl, nb.trials = 1)
#assign module membership to vertices
V(g_cntrl)$infomap = membership(modules_cntrl)

write.table(x = get.data.frame(x = g_cntrl, "vertices"), 
            file = "resultsWithOnlySymbl/dict_commInfomap_cntrl.txt", 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

#write out community structures
saveRDS(object = modules_cntrl, file = "resultsWithOnlySymbl/modules_cntrl.RDS")
saveRDS(object = modules_cases, file = "resultsWithOnlySymbl/modules_cases.RDS")
