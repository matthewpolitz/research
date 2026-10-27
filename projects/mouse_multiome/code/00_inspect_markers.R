# Summarise the Ernst et al. marker workbook (per reference cluster)
source("code/utils.R"); cfg <- load_config(); a <- cfg$annotation
mk <- read_ernst_markers(file.path(cfg$project_root, cfg$marker_xlsx))
print(table(mk$ref)[unique(mk$ref)])
sel <- subset(mk, FDR < a$fdr_max & Top <= a$max_top)
cat("\nGenes passing FDR < ", a$fdr_max, " & Top <= ", a$max_top, ":\n", sep = ""); print(table(sel$ref)[unique(mk$ref)])
print(head(sel[order(sel$ref, sel$Top, sel$FDR), ], 15))
