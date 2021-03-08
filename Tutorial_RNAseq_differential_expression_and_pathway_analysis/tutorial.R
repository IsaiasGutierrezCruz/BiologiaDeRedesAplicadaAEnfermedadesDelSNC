library(dplyr)
library(DESeq2)

# load the file 
countDataFile = "./GSE37704_featurecounts.csv"

# Import countdata
countData = read.csv(countDataFile, row.names=1) %>% 
  dplyr::select(-length) %>% 
  as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)

# Import metadata 
colData = read.csv("d:/ProyectosProgramacion/BiologiadeRedesSNC/Tutorial_RNAseq_differential_expression_and_pathway_analysis/GSE37704_metadata.csv", row.names=1)
colData

# Set up the DESeqDataSet Objetct and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
dds

# Get results for the HoxA1 knockdown versus control siRNA, and reorder them by p-value
# summary shows how many genes are up or down-regulated at FDR 0.1
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

# get the Entrez IDs, gene symbols, and full gene names
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
head(res, 10)


# Pathway analysis 

# get the “cleaner” gene sets of sinaling and metabolic pathways only
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

# set the vector of fold changes (the names of the values are the Entrez gene IDs)
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# get the results 
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

# get the upregulated pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways

# get the IDs
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# ---------make the plots --------
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

