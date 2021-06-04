# merge the metadata  
source("mergeMetadata.R")
FullMeta <- mergeMetadata() 

# merge the files of the counts 
source("mergeCounts.R")
RNAseqCounts <- mergeCounts(counts = "normalized")

# get the metadata of the necessary samples for the study
source("getSamplesToStudy.R")
samplesToStudy <- getSamplesToStudy(RNAseqCounts = RNAseqCounts, 
                                    FullMeta = FullMeta)

# format the data for the analysis
source("formatData.R")
data <- formatData(RNAseqCounts = RNAseqCounts, 
                   samplesToStudy = samplesToStudy)

countData <- data$countData
colData <- data$colData


# data analysis with edgeR
source("analysisDifferenceExpression.R")
top2 <- analysisDifferenceExpression(countData = countData, 
                             colData = colData)


# -------- gage -----------
# necesita: 


# ---- all the data wit kegg----

library("AnnotationDbi")
library("org.Hs.eg.db")

library(gage)
library(gageData)
library(pathview)

# remove the version numbers of the gene ID's
rowWithLength18 <- which((nchar(rownames(countData)) == 18) == TRUE)
rowWithLength17 <- which((nchar(rownames(countData)) == 17) == TRUE)

rownames(countData)[rowWithLength17] <- substr(rownames(countData)[rowWithLength17], 1, 
                                               nchar(rownames(countData)[rowWithLength17])-2)

rownames(countData)[rowWithLength18] <- substr(rownames(countData)[rowWithLength18], 1, 
                                               nchar(rownames(countData)[rowWithLength18])-3)


# change the format of the ID's 
nombres <- data.frame(id = rownames(countData)) 
nombres$entrez <- mapIds(org.Hs.eg.db,
                           keys=rownames(countData), 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")

rownames(countData) <- nombres$entrez


# remove NA values 
good <- complete.cases(rownames(countData))
countData <- countData[good, ]


data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
kegg.sets.hs <- kegg.sets.hs[1:131]

# groups 
earlyOnset <- which(colnames(countData)%in%samplesToStudy$LessThan70Years)
lateOnset <- which(colnames(countData)%in%samplesToStudy$GreaterThan70Years)

earlyOnset_v_LateOnsetSAMEDIR <- gage(exprs = countData, 
                               gsets = kegg.sets.hs, 
                               ref = earlyOnset,
                               samp = lateOnset, 
                               compare = "unpaired",
                               same.dir = TRUE)

earlyOnset_v_LateOnset.SigSAMEDIR <- sigGeneSet(earlyOnset_v_LateOnsetSAMEDIR)

earlyOnset_v_LateOnsetBOTHDIR <- gage(exprs = countData, 
                                      gsets = kegg.sets.hs, 
                                      ref = earlyOnset,
                                      samp = lateOnset, 
                                      compare = "unpaired",
                                      same.dir = FALSE)

earlyOnset_v_LateOnset.SigBOTHDIR <- sigGeneSet(earlyOnset_v_LateOnsetBOTHDIR, 
                                                outname = "earlyOnset_v_lateOnset")

tmp = sapply(c("hsa03010", "hsa00190", "hsa04260", "hsa04971", "hsa04971"), function(pid) pathview(gene.data=countData, 
                                                pathway.id=pid, species="hsa"))

# ---- all the data with go -----
# biological process (BP)
# cellular component (cc)
# molecular function (MF)

# pathways
data(go.sets.hs)
data("go.subs.hs")

earlyOnset_v_LateOnsetGO <- gage(exprs = countData, 
                               gsets = go.sets.hs[go.subs.hs$BP], 
                               ref = earlyOnset,
                               samp = lateOnset, 
                               compare = "unpaired")


# ---- after edgeR ----
# necesita: top2, 

library(dplyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

rowWithLength18 <- which((nchar(rownames(top2)) == 18) == TRUE)
rowWithLength17 <- which((nchar(rownames(top2)) == 17) == TRUE)

rownames(top2$table)[rowWithLength17] <- substr(rownames(top2$table)[rowWithLength17], 1, 
                                               nchar(rownames(top2$table)[rowWithLength17])-2)

rownames(top2$table)[rowWithLength18] <- substr(rownames(top2$table)[rowWithLength18], 1, 
                                               nchar(rownames(top2$table)[rowWithLength18])-3)


top2$table$symbol <- mapIds(org.Hs.eg.db, 
                      keys = row.names(top2),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

top2$table$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(top2), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

top2$table$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(top2), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")


top2$table


library(gage)
library(gageData)
library(pathview)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)


foldChanges <- top2$table$logFC
names(foldChanges) <- top2$table$entrez
head(foldChanges)


keggres <- gage(foldChanges, gsets = kegg.sets.hs, same.dir = TRUE)
lapply(keggres, head)

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tibble::as_tibble() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()

keggrespathways

# get the ID's
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldChanges, pathway.id=pid, species="hsa"))
