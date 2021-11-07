library(HTSanalyzeR)

test_commie = filter(commDict_cases, infomap==1)$name

my_enrichment = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = LIST_GO, 
                                               universe = V(g_cases)$name, 
                                               hits = test_commie, 
                                               minGeneSetSize = 1, 
                                               pAdjustMethod = "BH", 
                                               verbose = TRUE
)
print(nrow(my_enrichment))
