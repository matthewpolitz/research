# 04: cell type annotation using Ernst et al. (staged mouse spermatogenesis) cluster marker genes
# Approach: (a) AddModuleScore per reference cluster, (b) average module score per query cluster,
# (c) best-scoring reference label (with margin), (d) cross-check with canonical germ/somatic markers.
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(pheatmap) })
source("code/utils.R"); cfg <- load_config(); a <- cfg$annotation
outd <- res_dir(cfg, "04_annotation")
obj <- readRDS(file.path(cfg$data_path, "mouse_multiome.integrated.rds"))
DefaultAssay(obj) <- if (cfg$rna$normalisation == "SCT") "SCT" else "RNA"

# ---- markers ----
mk <- read_excel(cfg$marker_xlsx, sheet = a$sheet)
stopifnot(a$cluster_col %in% colnames(mk), a$gene_col %in% colnames(mk))
mk <- mk |> transmute(ref = as.character(.data[[a$cluster_col]]), gene = as.character(.data[[a$gene_col]]),
                      rank = if (nzchar(a$rank_col)) .data[[a$rank_col]] else NA)
if (nzchar(a$rank_col)) mk <- mk |> group_by(ref) |> arrange(desc(rank), .by_group = TRUE)  # edit direction if needed
mk <- mk |> group_by(ref) |> slice_head(n = a$top_n) |> ungroup()
genes_ok <- rownames(obj)
sets <- split(mk$gene, mk$ref)
sets <- lapply(sets, function(g) intersect(unique(g), genes_ok))
cat("Genes found per reference cluster:\n"); print(sapply(sets, length))
sets <- sets[sapply(sets, length) >= 5]

# ---- module scores ----
obj <- AddModuleScore(obj, features = sets, name = "ref_", ctrl = 50, seed = cfg$clustering$seed)
sc_cols <- paste0("ref_", seq_along(sets)); names(sc_cols) <- names(sets)
sc <- obj@meta.data[, sc_cols, drop = FALSE]; colnames(sc) <- paste0("score_", names(sets))
obj <- AddMetaData(obj, sc)

# per-cell best label (diagnostic) and per-cluster label (used)
obj$ref_label_cell <- names(sets)[max.col(as.matrix(sc), ties.method = "first")]
clus <- obj$seurat_clusters
avg <- aggregate(sc, list(cluster = clus), mean)
rownames(avg) <- avg$cluster; avg$cluster <- NULL; colnames(avg) <- names(sets)
best <- apply(avg, 1, function(x) { o <- order(x, decreasing = TRUE); c(names(x)[o[1]], x[o[1]] - x[o[2]]) })
lab <- data.frame(cluster = rownames(avg), ref_label = best[1, ], margin = as.numeric(best[2, ]))
lab$confident <- lab$margin > 0.05
lab$label <- ifelse(lab$confident, lab$ref_label, paste0(lab$ref_label, "?"))
write.csv(lab, file.path(outd, "cluster_label_proposal.csv"), row.names = FALSE)
write.csv(avg, file.path(outd, "cluster_x_reference_scores.csv"))
pheatmap(t(scale(avg)), filename = file.path(outd, "cluster_vs_reference_heatmap.pdf"), width = 9, height = 7)

obj$cell_type_auto <- lab$label[match(as.character(clus), lab$cluster)]

# ---- canonical marker sanity check (mouse testis) ----
canon <- list(
  Spermatogonia = c("Zbtb16", "Gfra1", "Id4", "Sall4"), Differentiating_SG = c("Kit", "Stra8", "Dmrt1"),
  Spermatocyte = c("Sycp3", "Sycp1", "Ccna1", "Prdm9", "Hormad1"), Round_spermatid = c("Acrv1", "Spaca1", "Tex101", "Catsper3"),
  Elongating_spermatid = c("Prm1", "Prm2", "Tnp1", "Tnp2", "Akap4"),
  Sertoli = c("Sox9", "Amh", "Wt1", "Cldn11"), Leydig = c("Cyp17a1", "Insl3", "Star", "Hsd3b1"),
  Peritubular = c("Acta2", "Myh11", "Dcn"), Macrophage = c("Adgre1", "Cd14", "C1qa"),
  Endothelial = c("Pecam1", "Cldn5", "Vwf"))
canon <- lapply(canon, intersect, genes_ok)
p <- DotPlot(obj, features = canon, group.by = "seurat_clusters") + RotatedAxis() + theme(axis.text.x = element_text(size = 7))
save_plot(p, file.path(outd, "canonical_markers_dotplot.png"), 16, 7)
save_plot(DimPlot(obj, reduction = "umap.wnn", group.by = "cell_type_auto", label = TRUE, repel = TRUE) , file.path(outd, "umap_wnn_auto_labels.png"), 10, 7)

# ---- optional manual overrides: edit docs/manual_annotation.csv (cluster,cell_type) after reviewing the above ----
man <- file.path(cfg$project_root, "docs", "manual_annotation.csv")
obj$cell_type <- obj$cell_type_auto
if (file.exists(man)) {
  m <- read.csv(man, stringsAsFactors = FALSE)
  obj$cell_type <- m$cell_type[match(as.character(clus), as.character(m$cluster))]
  log_msg("Applied manual annotation from ", man)
}
Idents(obj) <- "cell_type"
saveRDS(obj, file.path(cfg$data_path, "mouse_multiome.annotated.rds"))
write.csv(obj@meta.data[, c("sample", "seurat_clusters", "cell_type_auto", "cell_type")], file.path(outd, "cell_annotations.csv"))
session_log(outd)
