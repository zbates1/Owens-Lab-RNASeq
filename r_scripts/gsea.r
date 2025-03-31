# ===== Setting up environment =====
# Clean environment
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)
set.seed(123456)

# Install missing packages
required_packages <- c("tidyverse", "RColorBrewer", "fgsea", "msigdbr","msigdbdf", "data.table")
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages) > 0) {
  if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_packages)
}

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(msigdbdf)
  library(data.table)
})

# Set paths
base_dir <- getwd()
in_path <- file.path(base_dir, "Datasets")
out_path <- file.path(base_dir, "GSEA/Results")
bg_path <- file.path(base_dir, "GSEA/Background_genes")

# Create directories
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
dir.create(bg_path, showWarnings = FALSE, recursive = TRUE)

# ===== Data Preparation =====
# Read DESeq2 results
deseq_file <- file.choose() # Select your DESeq2 results CSV
df <- read.csv(deseq_file, row.names = 1) %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(gene_symbol = coalesce(gene_symbol, rownames(.)))

# Create ranked list
rankings <- df %>%
  mutate(rank = sign(log2FoldChange) * (-log10(pvalue))) %>%
  dplyr::select(gene_symbol, rank) %>%
  deframe() %>%
  sort(decreasing = TRUE)

# Handle infinite values
finite_vals <- rankings[is.finite(rankings)]
rankings[rankings == Inf] <- max(finite_vals) * 10
rankings[rankings == -Inf] <- min(finite_vals) * 10

# ===== Gene Set Preparation =====
prepare_gmt <- function(gmt_file, genes_in_data) {
  pathways <- gmtPathways(gmt_file)
  filtered_pathways <- lapply(pathways, function(p) intersect(p, genes_in_data))
  Filter(function(p) length(p) >= 10 & length(p) <= 500, filtered_pathways)
}

# Get background genes (example with MSigDB Hallmarks)
# Get Hallmark gene sets for mouse
msig_h <- msigdbr(
  species = "Mus musculus",
  category = "H",              # Still works but will give warning
  collection = "H"             # Preferred new syntax
) %>%
  dplyr::select(gs_name, gene_symbol)

# Save and load GMT (for reproducibility)
gmt_file <- file.path(bg_path, "msigdb_Hallmarks.gmt")
if(!file.exists(gmt_file)) {
  msig_h %>% 
    group_by(gs_name) %>% 
    summarise(gene_symbol = paste(gene_symbol, collapse = "\t")) %>%
    mutate(combined = paste(gs_name, "\t", gene_symbol)) %>%
    pull(combined) %>%
    writeLines(gmt_file)
}

pathways <- prepare_gmt(gmt_file, names(rankings))

# ===== Run GSEA =====
gsea_res <- fgsea(
  pathways = pathways,
  stats = rankings,
  scoreType = "std",
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000,
  eps = 0
)

# ===== Results Processing =====
# Collapse redundant pathways
collapsed <- collapsePathways(gsea_res[order(pval)], pathways, rankings)
main_pathways <- gsea_res[pathway %in% collapsed$mainPathways]

# Format results
formatted_res <- main_pathways %>%
  arrange(-NES) %>%
  mutate(
    leadingEdge = map_chr(leadingEdge, ~paste(.x, collapse = ",")),
    pathway = gsub("HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway)
  )

# ===== Visualization =====
# 1. Enrichment plot for top pathway
top_pathway <- formatted_res$pathway[1]
plotEnrichment(pathways[[top_pathway]], rankings) +
  labs(title = top_pathway, subtitle = paste("NES =", round(formatted_res$NES[1], 2))) +
  theme_classic()

ggsave(file.path(out_path, "top_pathway_enrichment.png"), width = 8, height = 6)

# 2. GSEA Table
top_pathways <- formatted_res$pathway[1:6]
plotGseaTable(pathways[top_pathways], rankings, gsea_res, 
              colwidths = c(5, 3, 0.8, 1.2, 1.2),
              gseaParam = 0.5)

# 3. NES Barplot
nes_plot <- formatted_res %>%
  head(20) %>%
  mutate(pathway = fct_reorder(pathway, NES)) %>%
  ggplot(aes(NES, pathway, fill = NES > 0)) +
  geom_col() +
  scale_fill_manual(values = c("blue", "red"), guide = "none") +
  labs(x = "Normalized Enrichment Score (NES)", y = "") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_path, "NES_barplot.png"), nes_plot, width = 10, height = 8)

# ===== Save Results =====
fwrite(formatted_res, file.path(out_path, "GSEA_results.csv"))
saveRDS(gsea_res, file.path(out_path, "GSEA_results.rds"))

# Session info
cat("\n=== Analysis Completed ===\n")
sessionInfo()