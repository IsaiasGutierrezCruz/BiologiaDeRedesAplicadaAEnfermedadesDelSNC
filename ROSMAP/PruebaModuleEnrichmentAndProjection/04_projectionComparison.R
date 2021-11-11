#5) Network comparison

#cases
#get matrices from graphs
xxx_cases = as.matrix(get.adjacency(flow_cases))
xxx_cases <- xxx_cases[c("1", "2", "24", "71", "211"),c("1", "2", "24", "71", "211")]
xxx_cases[1:10,1:10]
yyy_cases = as.matrix(get.adjacency(bp_cases))
yyy_cases[1:10,1:10]

#matrices may be unordered. As commnames are numbers, just sort
xxx_cases = xxx_cases[order(as.numeric(rownames(xxx_cases))),
          order(as.numeric(colnames(xxx_cases)))
          ]
xxx_cases[1:10,1:10]
yyy_cases = yyy_cases[order(as.numeric(rownames(yyy_cases))),
          order(as.numeric(colnames(yyy_cases)))
          ]
yyy_cases[1:10,1:10]

#to find which links are in both networks, multiply both matrices
#if there is a link in ModuleProjection and EnrichmentProjection, 
#product will be != 0
zzz_cases = xxx_cases*yyy_cases

my_replicated_links_cases = which(zzz_cases!=0, arr.ind = TRUE)

write.table(x = my_replicated_links_cases, 
            file = "results/replicatedLinks_cases.txt",
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)
#gives us the row and column
#which will be the community ids for communities linked in both networks

##check them
# neighbors(graph = mappy, v = "293")
# neighbors(graph = bp, v = "293")
# neighbors(graph = mappy, v = "326")
# neighbors(graph = bp, v = "326")
# 
# l_comm["293"]
# l_comm["134"]
# l_comm["326"]
# l_comm["137"]
# 

#cntrl
#get matrices from graphs
xxx_cntrl = as.matrix(get.adjacency(flow_cntrl))
xxx_cntrl[1:10,1:10]
yyy_cntrl = as.matrix(get.adjacency(bp_cntrl))
yyy_cntrl[1:10,1:10]

#matrices may be unordered. As commnames are numbers, just sort
xxx_cntrl = xxx_cntrl[order(as.numeric(rownames(xxx_cntrl))),
                      order(as.numeric(colnames(xxx_cntrl)))
                      ]
xxx_cntrl[1:10,1:10]
yyy_cntrl = yyy_cntrl[order(as.numeric(rownames(yyy_cntrl))),
                      order(as.numeric(colnames(yyy_cntrl)))
                      ]
yyy_cntrl[1:10,1:10]

#to find which links are in both networks, multiply both matrices
#if there is a link in ModuleProjection and EnrichmentProjection, 
#product will be != 0
zzz_cntrl = xxx_cntrl*yyy_cntrl

my_replicated_links_cntrl = which(zzz_cntrl!=0, arr.ind = TRUE)

write.table(x = my_replicated_links_cntrl, 
            file = "results/replicatedLinks_cntrl.txt",
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)
#gives us the row and column