source("function_bulkseq.R")

run_rnaseq_analysis(
  count_file = "../bulk_RNAseq/Col1GFP_Sepsis/genes.readcount.cntlD14.txt",
  columns_to_delete = "C_3_280N",
  conditions = c("control", "control", "treatment", "treatment", "treatment")
)

