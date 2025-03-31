# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)  # Mouse database
library(clusterProfiler)
library(biomaRt)

# Source the external file containing the function (if it's in a separate file)
source("C:/Users/ztba231/Documents/convert_gene_names.R")

############################################################
# Read count data for D14
#counts <- read.table("C:/Users/ztba231/Documents/genes.readcount.cntlD14.txt", 
#                     header = TRUE, 
#                     row.names = 1)
# Delete column, 'C_3_280N' as a control if you are in D14, because it is actually a D5 treatment
#counts <- subset(counts, select = -C_3_280N)


# Create metadata
#coldata <- data.frame(
#  row.names = colnames(counts),
#  condition = factor(c("control", "control",
#                       "treatment", "treatment", "treatment"))
#)
############################################################
# Read count data for D5
counts <- read.table("C:/Users/ztba231/Documents/genes.readcount.cntlD5.txt", 
                     header = TRUE, 
                     row.names = 1)

# Create metadata
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("control", "control", "treatment",
                       "treatment", "treatment", "treatment"))
)
############################################################

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Pre-filtering: keep only rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, alpha = 0.05, pAdjustMethod = "BH")
res_ordered <- res[order(res$padj),]
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)

# Convert Ensembl IDs to gene symbols and names
res_df <- convert_ensembl_to_symbol(res_df, "gene")

# Write results to file with gene names
write.csv(res_df, "deseq2_results_BH_with_names.csv")

# Get significant genes
sig_genes <- subset(res_df, padj < 0.05)
write.csv(sig_genes, "significant_genes_with_names.csv")

# GO Pathway Analysis
sig_genes_list <- sig_genes$gene[!is.na(sig_genes$gene)]

if (length(sig_genes_list) > 0) {
  go_enrichment <- enrichGO(
    gene = sig_genes_list,
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(go_enrichment) && nrow(as.data.frame(go_enrichment)) > 0) {
    # Get the genes in each GO term
    go_with_genes <- setReadable(go_enrichment, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
    
    # Extract gene clusters for each GO term
    gene_clusters <- data.frame(go_with_genes@result)
    
    # Write results with genes to CSV
    write.csv(gene_clusters, "go_enrichment_with_genes.csv")
    
    # Original dotplot
    png("go_dotplot.png", width = 800, height = 600)
    print(dotplot(go_enrichment, showCategory = 20))
    dev.off()
    
    # Create simplified GO network plot
    if (requireNamespace("enrichplot", quietly = TRUE)) {
      png("go_network.png", width = 1000, height = 800)
      print(enrichplot::emapplot(enrichplot::pairwise_termsim(go_enrichment)))
      dev.off()
      
      # Cnetplot shows genes connected to GO terms
      png("go_cnetplot.png", width = 1200, height = 900)
      print(enrichplot::cnetplot(go_with_genes, showCategory = 10, 
                                 foldChange = NULL, # Add fold change if available
                                 colorEdge = TRUE))
      dev.off()
    }
  }
}

# Volcano plot with gene symbols
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), size = 1) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("volcano_plot.png")

# Sample distance heatmap
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Matrix",
         file = "sample_distance_heatmap.png")

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()
ggsave("pca_plot.png")

# MA plot
png("ma_plot.png")
plotMA(res, ylim = c(-5,5))
dev.off()

# Get top differentially expressed genes
topGenes <- head(res_df[order(res_df$padj),], 20)
write.csv(topGenes, "top_20_DEGs.csv")

# Generate expression heatmap for top genes using gene names
valid_rows <- rownames(topGenes) %in% rownames(vsd)
if(sum(valid_rows) > 0) {
  mat <- assay(vsd)[rownames(topGenes)[valid_rows],]
  mat <- mat - rowMeans(mat)
  display_names <- topGenes$gene_symbol
  display_names[is.na(display_names)] <- rownames(topGenes)[is.na(display_names)]
  rownames(mat) <- display_names[valid_rows]
  pheatmap(mat,
           main = "Top 20 Differentially Expressed Genes",
           show_rownames = TRUE,
           scale = "row",
           file = "top_genes_heatmap.png")
}
