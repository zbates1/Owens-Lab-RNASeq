# Install and load required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(biomaRt)

convert_ensembl_to_symbol <- function(df, ensembl_column) {
  # Try connecting to Ensembl BioMart
  mart <- tryCatch({
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")
  }, error = function(e) {
    message("Failed to connect to Ensembl BioMart. Retrying with a different mirror...")
    useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://uswest.ensembl.org")
  })
  
  # Get gene symbols for the provided Ensembl IDs
  results <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = df[[ensembl_column]],
    mart = mart
  )
  
  # Add the new columns to the original dataframe
  df$gene_symbol <- results$external_gene_name[match(df[[ensembl_column]], results$ensembl_gene_id)]
  df$gene_name <- results$description[match(df[[ensembl_column]], results$ensembl_gene_id)]
  
  return(df)
}


# Example usage:
# Assuming your dataframe is called 'res_df' and the column with Ensembl IDs is called 'gene'
# res_df <- convert_ensembl_to_symbol(res_df, "gene")