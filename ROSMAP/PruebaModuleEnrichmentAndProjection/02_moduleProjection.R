#2) Module flow projection

#libraries
library(igraph)
library(data.table)
library(tidyverse)
library(dplyr)
source("MapFlow_2.R")

#cases
#Make projection to the community space
flow_cases = mapflow(g_cases, modules_cases)

#export 
write.graph(flow_cases, "resultsWithOnlySymbl/flow_cases.gml", "gml")

#cntrl
#Make projection to the community space
flow_cntrl = mapflow(g_cntrl, modules_cntrl)

#export 
write.graph(flow_cntrl, "resultsWithOnlySymbl/flow_cntrl.gml", "gml")
