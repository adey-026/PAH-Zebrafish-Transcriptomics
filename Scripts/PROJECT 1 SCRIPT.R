# Transcriptomic Analysis of PAH Exposure
# Zebrafish RNA-seq Reanalysis

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("EnhancedVolcano")
BiocManager::install("tximport")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("biomaRt")
library(clusterProfiler)
library(org.Dr.eg.db)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(tximport)
library(DESeq2)
library(biomaRt)
library(enrichplot)
setwd("D:/MISSION UK 2026/BIOINFORMATICS PROJECT/PROJECT 1 - Transcriptomic Fingerprint of PAH Pollution in Zebrafish Embryos/Data/ZF all data")
getwd()

# 1. Import Salmon Quantification Files

files <- list.files(pattern="quant.genes.sf", full.names=TRUE)
length(files)
sample_names <- gsub("_quant.genes.sf","", basename(files))
names(files) <- sample_names
files
txi <- tximport(files,type="salmon", txOut=TRUE)
dim(txi$counts)
head(txi$counts)

# 2. Prepare Count Matrix

counts <- txi$counts[,c("GSM6180939_Control_10hpf_10","GSM6180940_Control_10hpf_11","GSM6180941_Control_10hpf_12","GSM6180942_Control_10hpf_9","GSM6180943_BaP_10hpf_21","GSM6180944_BaP_10hpf_22","GSM6180945_BaP_10hpf_23","GSM6180946_BaP_10hpf_24")]
dim(counts)
condition <- c("control","control","control","control","BaP","BaP","BaP","BaP")
metadata <- data.frame(row.names = colnames(counts),condition)
metadata
counts <- round(counts)
head(counts)

# 3. Differential Expression Analysis

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)
res <- results(dds)
nrow(dds)

# 3. Differential Expression Analysis

vsd <- vst(dds)
plotPCA(vsd, intgroup="condition")

# 4. significant genes analysis

summary(res)
sig_genes <- subset(res, padj < 0.05)
nrow(sig_genes)

# 5. Volcano Plot

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
EnhancedVolcano(res,lab = rownames(res), x = 'log2FoldChange', y = 'padj')

# 6. Heatmap of Top Genes

library(pheatmap)
topgenes <- head(order(res$padj),50)
pheatmap(assay(vsd)[topgenes,], scale="row")


# 7. Heatmap of Top Genes

gene_list <- rownames(sig_genes)

gene_list <- gsub("\\..*", "", gene_list)
head(gene_list)


mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),filters = "ensembl_transcript_id", values = gene_list, mart = mart)
mart <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl", mirror = "useast")
head(gene_list)
gene_ids <- gsub("ENSDART", "ENSDARG", gene_list)
head(gene_ids)
gene_ids <- unique(gene_ids)

length(gene_ids)
sig_genes <- subset(res, padj < 0.05)
gene_list <- rownames(sig_genes)
gene_list
summary(res)
length(gene_ids)
nrow(sig_genes)

gene_df <- bitr(gene_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
nrow(gene_df)
entrez_genes <- unique(gene_df$ENTREZID)
length(entrez_genes)

ego <- enrichGO(gene = entrez_genes, OrgDb = org.Dr.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.2, qvalueCutoff = 0.2)
head(ego)
dotplot(ego, showCategory = 20)
dotplot(ego, showCategory = 10, font.size = 10)


# 8. KEGG analysis

dotplot(ego, showCategory = 10) +
  ggplot2::theme(axis.text.y = element_text(size = 8))
kegg <- enrichKEGG(gene = entrez_genes, organism = "dre",pvalueCutoff = 0.2)
as.data.frame(kegg)
dotplot(kegg, showCategory = 10)

# Save Results 
write.csv(sig_genes,"significant_genes.csv")









