# Load required libraries
library(DESeq2)
library(tidyverse)

# Read count data
# Replace 'genes.readcount' with your file path
counts <- read.table("C:/Users/ztba231/Documents/genes.readcount.xls", header=TRUE, row.names=1)

print(counts.shape)

counts <- counts[,-1]

# Create the coldata data frame
coldata <- data.frame(
  row.names = colnames(counts)[-1],  # Exclude the first column (geneID)
  condition = factor(c("control", "control", "control",
                       "treatment", "treatment", "treatment",
                       "treatment", "treatment"))
)

# Print coldata to verify
print(coldata)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filter low count genes (optional but recommended)
# Keep genes with at least 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results with BH adjusted p-values (padj)
res <- results(dds, alpha = 0.05, pAdjustMethod = "BH")

# Sort by adjusted p-value
res_ordered <- res[order(res$padj),]

# Convert to data frame for easier handling
res_df <- as.data.frame(res_ordered)

# Add gene names as column
res_df$gene <- rownames(res_df)

# Write results to file
write.csv(res_df, "deseq2_results_BH.csv")

# Basic summary of results
summary(res)

# Get significantly differentially expressed genes
sig_genes <- subset(res_df, padj < 0.05)

# Create volcano plot using res_df instead of res
plot(res_df$log2FoldChange, -log10(res_df$pvalue),
     pch=20, col="gray",
     xlab="Log2 Fold Change",
     ylab="-Log10 P-value")

# Add points for significant genes
points(res_df$log2FoldChange[res_df$padj < 0.05],
       -log10(res_df$pvalue[res_df$padj < 0.05]),
       pch=20, col="red")
