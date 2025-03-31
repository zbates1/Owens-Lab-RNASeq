# ===== Setting up environment =====
# Clean environment
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)

# Install missing packages if needed
required_packages <- c("tidyverse", "RColorBrewer", "pheatmap", "clusterProfiler", 
                       "org.Mm.eg.db", "DOSE", "enrichplot", "ggupset")
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages) > 0) {
  if("ggupset" %in% new_packages) {
    install.packages("ggupset", repos = "http://cran.us.r-project.org", dependencies = TRUE)
  }
  if(any(c("clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot") %in% new_packages)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_packages[new_packages %in% c("clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot")])
  }
  install.packages(new_packages[!new_packages %in% c("clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "ggupset")])
}

# Load libraries with error handling
for (pkg in required_packages) {
  try(library(pkg, character.only = TRUE), silent = TRUE)
}

# Set paths with robustness checks
base_dir <- getwd()
in_path <- file.path(base_dir, "Datasets")
out_path <- file.path(base_dir, "PEA/Results")
bg_path <- file.path(base_dir, "PEA/Background_genes")

# Create directories if they don't exist
dir.create(file.path(base_dir, "PEA"), showWarnings = FALSE)
dir.create(out_path, showWarnings = FALSE)
dir.create(file.path(out_path, "clusterProfiler"), showWarnings = FALSE)
dir.create(bg_path, showWarnings = FALSE)

# ===== Read and prepare data =====
# Define file path with user dialog option
deseq_file <- "C:/Users/ztba231/Documents/deseq2_results_BH_with_names.csv"
if (!file.exists(deseq_file)) {
  deseq_file <- file.choose() # Opens file dialog if path doesn't exist
}

# Read DESeq2 results
df <- read.csv(deseq_file, row.names = 1)
print(paste("Loaded", nrow(df), "genes from DESeq2 results"))

# Annotate differentially expressed genes
df <- df %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  TRUE ~ 'NO'
))

# Get summary of DE genes
de_summary <- table(df$diffexpressed)
print("Differential expression summary:")
print(de_summary)

# ===== Pathway database preparation =====
# Extract genes present in data
genes_in_data <- df$gene_symbol

# Process GMT file with user dialog option
gmt_file <- "C:/Users/ztba231/Documents/mh.all.v2024.1.Mm.symbols.gmt"
if (!file.exists(gmt_file)) {
  gmt_file <- file.choose() # Opens file dialog if path doesn't exist
}

# Process GMT file and extract database name
db_name <- basename(gmt_file) %>% 
  gsub("\\.gmt$", "", .) %>%
  strsplit("\\.")
db_name <- db_name[[1]][1]  # Extract first element for short name

# Read and filter GMT file
pwl2 <- read.gmt(gmt_file)
pwl2_filtered <- pwl2[pwl2$gene %in% genes_in_data,]

# Save filtered gene set
filtered_file <- file.path(bg_path, paste0(db_name, ".RDS"))
saveRDS(pwl2_filtered, filtered_file)
print(paste("Saved filtered", db_name, "database with", length(unique(pwl2_filtered$gene)), "genes"))

# ===== Run ClusterProfiler =====
# Remove non-significant genes
df_sig <- df[df$diffexpressed != 'NO', ]
print(paste("Performing enrichment analysis on", nrow(df_sig), "DE genes"))

# Settings
name_of_comparison <- 'treatment_vs_control'
padj_cutoff <- 0.05
genecount_cutoff <- 5

# Split by expression direction
deg_results_list <- split(df_sig, df_sig$diffexpressed)

# Run enrichment analysis
res <- lapply(names(deg_results_list), function(x) {
  print(paste("Running enrichment for", x, "genes"))
  enricher(gene = deg_results_list[[x]]$gene_symbol,
           TERM2GENE = pwl2_filtered,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.1,  # More lenient for initial analysis
           qvalueCutoff = 0.1)
})
names(res) <- names(deg_results_list)

# Check if any results were found
if (all(sapply(res, function(x) nrow(x@result) == 0))) {
  stop("No enriched pathways found. Consider relaxing thresholds or using a different database.")
}

# Process results
res_df <- lapply(names(res), function(x) {
  if (nrow(res[[x]]@result) > 0) {
    result <- res[[x]]@result
    result$direction <- x
    return(result)
  } else {
    return(NULL)
  }
})
res_df <- do.call(rbind, res_df[!sapply(res_df, is.null)])

# Add useful columns
res_df <- res_df %>% 
  mutate(minuslog10padj = -log10(p.adjust),
         GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))

# Filter significant pathways
sig_res_df <- res_df[res_df$p.adjust < padj_cutoff & as.numeric(sapply(strsplit(res_df$GeneRatio, "/"), `[`, 1)) > genecount_cutoff, ]
print(paste("Found", nrow(sig_res_df), "significantly enriched pathways"))

# ===== Export results =====
# Save complete results
complete_results_file <- file.path(out_path, "clusterProfiler", paste0(name_of_comparison, "_", db_name, "_complete.csv"))
write.csv(res_df, complete_results_file, row.names = FALSE)

# Save significant results
sig_results_file <- file.path(out_path, "clusterProfiler", paste0(name_of_comparison, "_", db_name, "_significant.csv"))
write.csv(sig_res_df, sig_results_file, row.names = FALSE)

# ===== Visualization =====
# Check if there are results to visualize
if (nrow(sig_res_df) > 0) {
  # Make sure the output directory exists
  plots_dir <- file.path(out_path, "clusterProfiler", "plots")
  dir.create(plots_dir, showWarnings = FALSE)
  
  # 1. Dotplot of top pathways
  top_n <- min(20, nrow(sig_res_df))
  top_pathways <- sig_res_df %>%
    group_by(direction) %>%
    top_n(top_n, -p.adjust)
  
  dot_plot <- ggplot(top_pathways, aes(x = direction, y = Description, size = GeneRatio_numeric, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(title = "Top Enriched Pathways", x = "Regulation Direction", y = "Pathway") +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(plots_dir, paste0(name_of_comparison, "_", db_name, "_dotplot.png")), 
         dot_plot, width = 10, height = 8)
  
  # 2. Network plot of pathways using enrichplot
  for (direction in unique(sig_res_df$direction)) {
    dir_results <- res[[direction]]
    if (nrow(dir_results@result) > 0) {
      # Create and save emapplot if more than one pathway
      if (nrow(dir_results@result) > 1) {
        try({
          emap <- emapplot(pairwise_termsim(dir_results), showCategory = min(30, nrow(dir_results@result)))
          ggsave(file.path(plots_dir, paste0(name_of_comparison, "_", db_name, "_", direction, "_network.png")), 
                 emap, width = 12, height = 10)
        }, silent = TRUE)
      }
      
      # Create and save cnetplot
      try({
        genes_to_use <- deg_results_list[[direction]]$gene_symbol
        cnet <- cnetplot(dir_results, showCategory = min(10, nrow(dir_results@result)), 
                         foldChange = deg_results_list[[direction]]$log2FoldChange,
                         categorySize = "pvalue")
        ggsave(file.path(plots_dir, paste0(name_of_comparison, "_", db_name, "_", direction, "_genenet.png")), 
               cnet, width = 12, height = 10)
      }, silent = TRUE)
    }
  }
  
  print(paste("Visualizations saved to", plots_dir))
}

# Print completion message
print(paste("Pathway enrichment analysis completed. Results saved to", out_path))

