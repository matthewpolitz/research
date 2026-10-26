# Dump structure of the Ernst et al. marker workbook so annotation columns can be set in config.yml
source("code/utils.R"); cfg <- load_config(); library(readxl)
sheets <- excel_sheets(cfg$marker_xlsx); print(sheets)
for (s in sheets) {
  d <- read_excel(cfg$marker_xlsx, sheet = s, n_max = 8)
  cat("\n=== sheet:", s, "===\n"); print(colnames(d)); print(as.data.frame(d))
}
cat("\nSet annotation$sheet, $cluster_col, $gene_col (and optionally $rank_col) in config/config.yml\n")
