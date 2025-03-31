
        # Function to convert Ensembl IDs to gene symbols
        convert_ensembl_to_symbol <- function(data_frame, id_column) {
          # Check if biomaRt is installed
          if (!requireNamespace("biomaRt", quietly = TRUE)) {
            warning("biomaRt package is required but not installed.")
            return(data_frame)
          }
          
          # Connect to Ensembl
          ensembl <- tryCatch({
            biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
          }, error = function(e) {
            warning("Could not connect to Ensembl. Using offline method.")
            return(NULL)
          })
          
          if (!is.null(ensembl)) {
            # Get Ensembl IDs from data frame
            ensembl_ids <- data_frame[[id_column]]
            
            # Get gene information
            gene_info <- biomaRt::getBM(
              attributes = c("ensembl_gene_id", "external_gene_name", "description"),
              filters = "ensembl_gene_id",
              values = ensembl_ids,
              mart = ensembl
            )
            
            # Add information to data frame
            data_frame$gene_symbol <- NA
            data_frame$gene_description <- NA
            
            # Match and update
            for (i in 1:nrow(gene_info)) {
              idx <- which(data_frame[[id_column]] == gene_info$ensembl_gene_id[i])
              if (length(idx) > 0) {
                data_frame$gene_symbol[idx] <- gene_info$external_gene_name[i]
                data_frame$gene_description[idx] <- gene_info$description[i]
              }
            }
          } else {
            # Fallback method
            data_frame$gene_symbol <- data_frame[[id_column]]
            data_frame$gene_description <- NA
          }
          
          return(data_frame)
        }
        