# 05: collect QC and composition tables
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
source("code/utils.R"); cfg <- load_config(); outd <- res_dir(cfg, "05_summary")
obj <- readRDS(file.path(cfg$data_path, "mouse_multiome.annotated.rds"))
tab <- as.data.frame(table(sample = obj$sample, cell_type = obj$cell_type)) |> group_by(sample) |> mutate(frac = Freq / sum(Freq))
write.csv(tab, file.path(outd, "composition_by_sample.csv"), row.names = FALSE)
p <- ggplot(tab, aes(sample, frac, fill = cell_type)) + geom_col() + theme_bw()
save_plot(p, file.path(outd, "composition.png"), 8, 6)
