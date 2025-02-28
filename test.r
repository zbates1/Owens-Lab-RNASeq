# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)

# Read count data
counts <- read.table("C:/Users/ztba231/Documents/genes.readcount.cntlD14.txt", 
                     header = TRUE, 
                     row.names = 1)

# Create metadata
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("control", "control", "control",
                       "treatment", "treatment", "treatment"))
)

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

# Convert gene IDs using AnnotationDbi instead of biomaRt
res_df$gene_symbol <- mapIds(org.Hs.eg.db,
                             keys = res_df$gene,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")

res_df$gene_name <- mapIds(org.Hs.eg.db,
                           keys = res_df$gene,
                           column = "GENENAME",
                           keytype = "ENSEMBL",
                           multiVals = "first")

# Write results to file with gene names
write.csv(res_df, "deseq2_results_BH_with_names.csv")

# Get significant genes
sig_genes <- subset(res_df, padj < 0.05)
write.csv(sig_genes, "significant_genes_with_names.csv")

# GO Pathway Analysis
# Get significant gene list and remove NAs
sig_genes_list <- sig_genes$gene[!is.na(sig_genes$gene)]

# Check if we have any significant genes
if(length(sig_genes_list) > 0) {
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(
    gene = sig_genes_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # Save GO results if we got any
  if(!is.null(go_enrichment) && nrow(as.data.frame(go_enrichment)) > 0) {
    write.csv(as.data.frame(go_enrichment), "go_enrichment_results.csv")
    
    # Plot top GO terms
    png("go_dotplot.png", width = 800, height = 600)
    print(dotplot(go_enrichment, showCategory = 20))
    dev.off()
    
    # Create simplified GO network plot
    if(requireNamespace("enrichplot", quietly = TRUE)) {
      png("go_network.png", width = 1000, height = 800)
      print(enrichplot::emapplot(enrichplot::pairwise_termsim(go_enrichment)))
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
# Make sure we have valid row names
valid_rows <- rownames(topGenes) %in% rownames(vsd)
if(sum(valid_rows) > 0) {
  mat <- assay(vsd)[rownames(topGenes)[valid_rows],]
  mat <- mat - rowMeans(mat)
  # Use gene symbols where available, otherwise keep original IDs
  display_names <- topGenes$gene_symbol
  display_names[is.na(display_names)] <- rownames(topGenes)[is.na(display_names)]
  rownames(mat) <- display_names[valid_rows]
  pheatmap(mat,
           main = "Top 20 Differentially Expressed Genes",
           show_rownames = TRUE,
           scale = "row",
           file = "top_genes_heatmap.png")
}