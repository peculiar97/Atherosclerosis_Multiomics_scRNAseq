## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
library(DESeq2)
library(ggplot2) 
library(knitr)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(textshaping)
library(tibble)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata = read.table("/new_workflow/results/5_final_counts/MQ_final_counts.txt", header = T, skip = 1, row.names = 1)
head(countdata)

colnames(countdata) = gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) = gsub("Aligned.sortedByCoord.out", "", colnames(countdata), fixed = T)
colnames(countdata) = gsub("..", "", colnames(countdata), fixed = T)

countdata = countdata[ , c(-1:-5)]
colnames(countdata)

countdata = countdata[, c(2, 8, 14, 1, 7, 13, 4, 10, 16, 5, 11, 17, 6, 12, 18, 3, 9,15)]


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Filtering out low counts

# Exclude rows with count number = 0
countdata_woZeros = countdata[rowSums(countdata)>0, ]

# Exclude rows that dont have at least a count number > 10
countdata_atLeast10 = countdata_woZeros %>% rownames_to_column("gene") %>% filter_if(is.numeric, any_vars(. > 10)) %>% column_to_rownames("gene")
# source: https://stackoverflow.com/questions/45427444/how-to-preserve-base-data-frame-rownames-upon-filtering-in-dplyr-chain


## --------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare metadata file
id = paste0(rep(paste0(rep("MQ_S", 3), 1:3), 6), rep(c("_0h", "_05h", "_2h", "_4h", "_8h", "_24h"), each=3))
sample_info = data.frame(id= id ,
                         Time=factor(c(rep("UN", 3), rep("0.5h", 3), rep("2h", 3), rep("4h",3), rep("8h", 3), rep("24h", 3)), 
                                     levels = c("UN", "0.5h", "2h", "4h", "8h", "24h")),
                         Replicate= factor(rep(c("1", "2" , "3"), 6), levels = c("1", "2", "3")), 
                         row.names = 1)
sample_info$sampleid = row.names(sample_info)

sample_info = sample_info %>% rownames_to_column() %>% arrange(Time) %>% column_to_rownames()

sample_info

# Reorder sampleID's to match featureCounts column order.
sample_info = sample_info[match(colnames(countdata), sample_info$sampleid),]

# Make sure ID's are correct
DT::datatable(sample_info)



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Make DESeq2 object from counts and metadata
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. See https://support.bioconductor.org/p/67600/#67612
#                                                      https://www.biostars.org/p/325009/    DESeq2 compare all levels
#   Use (~) before the name of the column variable to compare

ddsMat = DESeqDataSetFromMatrix(countData = countdata_atLeast10, colData = sample_info, design = ~ Replicate + Time)

ddsMat$Time = relevel(ddsMat$Time, ref = "UN")
ddsMat$Time

# Check the expression data
assay(ddsMat) %>% head(10)

# Check the sample info
colData(ddsMat) %>% head(20)

# Check the difference before/after DESeq2 running
setdiff(colSums(counts(ddsMat)), colSums(countdata_atLeast10))


# Find differential expressed genes
ddsMat = DESeq(ddsMat, test = "LRT", reduced = ~Replicate)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#Get basic statisics about the number of significant genes

results = results(ddsMat, alpha = 0.05)

# Generate summary of testing
summary(results)

# # Check directionality of the log2 fold changes
view(mcols(results, use.names = T))



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Results from different time points
resultsNames(ddsMat)

# prepare the tibble
for(i in resultsNames(ddsMat)){
  res = results(ddsMat, name=i, alpha = 0.05)
  assign(i, res %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble())
    }
res_0.5_2 = full_join(Time_0.5h_vs_UN, Time_2h_vs_UN, by="gene") %>% dplyr::select(gene, log2FoldChange.x, log2FoldChange.y) %>% 
  dplyr::rename(c(log2FC_0.5h=log2FoldChange.x, log2FC_2h=log2FoldChange.y ))

res_4_8_24 = full_join(Time_4h_vs_UN, Time_8h_vs_UN, by="gene") %>% dplyr::select(gene, log2FoldChange.x, log2FoldChange.y) %>%
  dplyr::rename(c(log2FC_4h = log2FoldChange.x, log2FC_8h = log2FoldChange.y)) %>% full_join(., Time_24h_vs_UN, by="gene") %>% 
  dplyr::rename(log2FC_24h = log2FoldChange)

res_all = full_join(res_0.5_2, res_4_8_24, by="gene")  
res_all = res_all[, c(1,2,3,4,5,7,6,8,9,10,11)]


# more on joining 
#https://stackoverflow.com/questions/32066402/how-to-perform-multiple-left-joins-using-dplyr-in-r


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Finding specific genes
v = c("Stat1", "Stat2", "Irf1", "Irf9", "Isg15")

Time_2h_vs_UN %>% filter(gene %in% v)
#OR
Time_4h_vs_UN %>% filter_all(any_vars(. %in% v))



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Gather gene annotation information
library(org.Mm.eg.db)

# Add gene full name
res_all$description = mapIds(x = org.Mm.eg.db, keys = res_all$gene, column = "GENENAME", keytype = "SYMBOL", multiVals = "first")

# Add ENSEMBL
res_all$ensembl = mapIds(x = org.Mm.eg.db, keys = res_all$gene, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

# Add ENTREZ ID
res_all$entrez = mapIds(x = org.Mm.eg.db, keys = res_all$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Subset for only significat genes
res_all_padj05 = subset(res_all, padj < 0.05)

# write significant annotated results to a file
write.table(x = as.data.frame(res_all_padj05), file= "MQ_results/MQ_results_annotated_significant.tsv",
            sep = '\t', quote = F, col.names = NA)



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Sample distances

#1 
ddsMat_rlog = rlog(ddsMat, blind = FALSE)

sampleDists <- dist( t( assay(ddsMat_rlog) ) )
sampleDistMatrix <- as.matrix( sampleDists )
#rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
#colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(file = "MQ_sample_distance.pdf")

pheatmap(sampleDistMatrix,
	  clustering_distance_rows=sampleDists,
	  clustering_distance_cols=sampleDists,
	  col=colors)

dev.off()

#2 using Poisson Distance

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsMat)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colData(ddsMat)$sampleid
colnames(samplePoisDistMatrix) <- colData(ddsMat)$sampleid

pdf(file = "MQ_Poisson_distance.pdf")
pheatmap(samplePoisDistMatrix,
	  clustering_distance_rows=poisd$dd,
	  clustering_distance_cols=poisd$dd,
	  col=colors)

dev.off()



## --------------------------------------------------------------------------------------------------------------------------------------------------
# PCA plot

ddsMat_rlog = rlog(ddsMat, blind = FALSE)

# par(mfrow= c(1,2))

pcaData = plotPCA(ddsMat_rlog, intgroup = c("Time", "Replicate"), ntop=500, returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = "MQ_PCA_plot.pdf")
ggplot(pcaData, aes(PC1, PC2, color=Time, shape=Replicate)) +   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance"))

dev.off()

plotPCA(ddsMat_rlog, intgroup = "Time", ntop=500) +
  theme_bw() +
  geom_point(size=5) + 
  ggtitle(label = "PCA", subtitle = "Top 500 most variable genes")
  
plotPCA(ddsMat_rlog, intgroup = "Replicate", ntop=500) +
  theme_bw() + 
  geom_point(size=5)


# dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Heatmap plot
ddsMat_rlog = rlog(ddsMat, blind = FALSE)

# Gather 40 most significatnt genes
mat = assay(ddsMat_rlog[order(res_all_padj05$padj, decreasing = F),])[1:50,]

annotation_col = data.frame(
  Time = factor(colData(ddsMat_rlog)$Time),
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

pdf(file = "MQ_most_sig_genes.pdf")

# Scale genes to Z-score (how many standard deviations)
pheatmap(mat=mat, scale = "row", annotation_col = annotation_col, fontsize = 6.5)

dev.off()

#########################
# the most highly variable genes.

library(ComplexHeatmap)

topVarGenes <- head(order(rowVars(assay(ddsMat_rlog)),decreasing=TRUE),50)
mat <- assay(ddsMat_rlog)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
#df = as.data.frame(colData(ddsMat_rlog)) %>% select(c(Time, Replicate)) %>% rownames_to_column() %>% arrange(Time) %>% column_to_rownames() # First Arrange rows
df = as.data.frame(colData(ddsMat_rlog)[, c("Time" , "Replicate")])

pdf(file = "MQ_most_variable_genes.pdf")
ComplexHeatmap::pheatmap(mat, annotation_col=df, cluster_cols = T, cluster_rows = T, fontsize_row = 8)

dev.off()



## --------------------------------------------------------------------------------------------------------------------------------------------------
#plotMA
 
# Remember that the res object is  " Time 24h vs UN""
plotMA(res, ylim = c(-5,5))



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Single gene Plot

ddsMat_rlog = rlog(ddsMat, blind = FALSE)

plotCounts(ddsMat, gene = "Cxcl10", intgroup = c("Time", "Replicate"), normalized = T, transform = T)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
#Preparing the gene list for clustering 
res_padjE_3 = res_all_padj05 %>% as.data.frame() %>% filter(padj < 0.001) %>% arrange(padj) %>% head(1000)



## --------------------------------------------------------------------------------------------------------------------------------------------------
library(DEGreport)
#First, the lasso package should be installed using install.packages("lasso2_1.2-22.tar.gz", type = "source", repos = NULL)

ddsMat_rlog = rlog(ddsMat, blind = FALSE)
ddsMat_rlog_counts = assay(ddsMat_rlog)

cluster_rlog_padjE_3 = ddsMat_rlog_counts[res_padjE_3$gene,]


cl_rlog_padjE_3 = degPatterns(ma = cluster_rlog_padjE_3, metadata = sample_info, minc = 40, time = "Time", reduce = F)


#Extracting & saving 
cl_rlog_padjE_3$df %>% filter(cluster == 2) %>% select(genes) %>% write.table("MQ_cl2.tsv", sep = "\t", quote = F, row.names = F, col.names = F )


# Extracting the genes in cluster 1
cluster_1 = cl_rlog_padjE_3$df %>% filter( cluster == 1)

# Extracting our interested genes
cl_rlog_padjE_3$df %>% filter( genes %in% v)




## --------------------------------------------------------------------------------------------------------------------------------------------------
# Finding pathways from differential expressed genes

#Remove any genes that do not have any entrez identifiers
results_sig_entrez = subset(res_all_padj05, is.na(entrez) == FALSE)

# Create a data frame for the average of log2FC
df_meanslog2FC = data.frame(gene=results_sig_entrez$gene, log2FC = rowMeans(results_sig_entrez[, 2:6]))

gene_matrix = df_meanslog2FC$log2FC

names(gene_matrix) = results_sig_entrez$entrez

###Enrichment analysis using KEGG
kegg_enrich = enrichKEGG(gene = names(gene_matrix), organism = "mouse", qvalueCutoff = 0.10)
barplot(kegg_enrich, drop = T, showCategory = 10, font.size = 8)

###Enrichment analysis using GO
go_enrich = enrichGO(gene = names(gene_matrix), OrgDb = "org.Mm.eg.db", readable = T, ont = "BP", qvalueCutoff = 0.10)
barplot(go_enrich, drop =T, showCategory = 10, font.size = 8)

###Plotting KEGG pathway
library(pathview)

pathview(gene.data = gene_matrix, pathway.id = "mmu05014", species = "mouse")


## --------------------------------------------------------------------------------------------------------------------------------------------------
for (i in 1:5){ cl = cl_rlog_padjE_3$df %>% filter(cluster == i) 
cl$entrezid = mapIds(x = org.Mm.eg.db, keys = cl$genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
assign(x=paste0("cl_", i), cl)
}

myls = list(cl_1$entrezid, cl_2$entrezid, cl_3$entrezid, cl_4$entrezid, cl_5$entrezid)
names(myls) = paste0("X", 1:5)

ck1 = compareCluster(geneClusters = myls, OrgDb = org.Mm.eg.db )
ck1 = setReadable(ck1, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

dotplot(ck1)


## --------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()


