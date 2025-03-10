# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!require("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
if (!require("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")

library(biomaRt)
library(AnnotationDbi)
library(org.Mm.eg.db)

convert_ensembl_to_symbol <- function(df, ensembl_column) {
  # Define a list of potential Ensembl hosts to try
  hosts <- c(
    "https://www.ensembl.org",
    "https://useast.ensembl.org",
    "https://uswest.ensembl.org",
    "https://asia.ensembl.org",
    "https://ensembl.org"
  )
  
  # Try each host until one works
  for (host in hosts) {
    message(paste("Trying to connect to", host, "..."))
    mart <- tryCatch({
      useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    }, error = function(e) {
      message(paste("Failed to connect to", host))
      return(NULL)
    })
    
    if (!is.null(mart)) {
      message(paste("Successfully connected to", host))
      # Try to query the mart
      results <- tryCatch({
        getBM(
          attributes = c("ensembl_gene_id", "external_gene_name", "description"),
          filters = "ensembl_gene_id",
          values = df[[ensembl_column]],
          mart = mart
        )
      }, error = function(e) {
        message(paste("Failed to query", host, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(results) && nrow(results) > 0) {
        # Add the new columns to the original dataframe
        df$gene_symbol <- results$external_gene_name[match(df[[ensembl_column]], results$ensembl_gene_id)]
        df$gene_name <- results$description[match(df[[ensembl_column]], results$ensembl_gene_id)]
        return(df)
      }
    }
  }
  
  # If all BioMart attempts failed, try using org.Mm.eg.db as fallback
  message("All BioMart attempts failed. Trying local annotation database...")
  tryCatch({
    # Remove version numbers from Ensembl IDs if present (e.g., ENSMUSG00000000001.5 -> ENSMUSG00000000001)
    clean_ids <- sub("\\.[0-9]+$", "", df[[ensembl_column]])
    
    # Map Ensembl IDs to gene symbols
    gene_symbols <- mapIds(org.Mm.eg.db, 
                           keys = clean_ids,
                           column = "SYMBOL", 
                           keytype = "ENSEMBL",
                           multiVals = "first")
    
    # Map Ensembl IDs to gene names/descriptions
    gene_names <- mapIds(org.Mm.eg.db, 
                         keys = clean_ids,
                         column = "GENENAME", 
                         keytype = "ENSEMBL",
                         multiVals = "first")
    
    df$gene_symbol <- gene_symbols[match(clean_ids, names(gene_symbols))]
    df$gene_name <- gene_names[match(clean_ids, names(gene_names))]
    
    message("Successfully mapped IDs using local database.")
    return(df)
  }, error = function(e) {
    message(paste("Local database mapping failed:", e$message))
    message("Returning original dataframe without gene symbol mapping.")
    return(df)
  })
}