# 03: merge samples, re-normalise jointly, Harmony on RNA PCA and ATAC LSI, WNN, clustering
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(harmony); library(ggplot2); library(patchwork) })
source("code/utils.R"); cfg <- load_config()
outd <- res_dir(cfg, "03_integration"); set.seed(cfg$clustering$seed)
rd <- parse_dims(cfg$rna$dims); ad <- parse_dims(cfg$atac$lsi_dims)
gb <- cfg$integration$group_by

objs <- lapply(cfg$samples, function(s) readRDS(file.path(cfg$data_path, paste0(s, ".norm.rds"))))
for (i in seq_along(objs)) { DefaultAssay(objs[[i]]) <- "RNA"; objs[[i]][["SCT"]] <- NULL
  objs[[i]] <- DietSeurat(objs[[i]], assays = c("RNA", "peaks"), dimreducs = NULL) }
obj <- merge(objs[[1]], objs[-1], add.cell.ids = cfg$samples)
rm(objs); gc()
obj[[gb]] <- obj$sample

# ---- RNA ----
DefaultAssay(obj) <- "RNA"
if (cfg$rna$normalisation == "SCT") {
  obj <- SCTransform(obj, vars.to.regress = unlist(cfg$rna$regress_out), variable.features.n = cfg$rna$n_variable_features, verbose = FALSE)
} else {
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = cfg$rna$n_variable_features) |> ScaleData()
}
obj <- RunPCA(obj, npcs = cfg$rna$npcs, reduction.name = "pca", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = rd, reduction.name = "umap.rna.unint", reduction.key = "rnaUMAPu_", verbose = FALSE)
obj <- harmony::RunHarmony(obj, group.by.vars = gb, reduction.use = "pca", dims.use = rd,
                           theta = cfg$integration$rna_theta, reduction.save = "harmony.rna", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony.rna", dims = rd, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE)

# ---- ATAC ----
DefaultAssay(obj) <- "peaks"
obj <- FindTopFeatures(obj, min.cutoff = cfg$atac$min_cells_feature)
obj <- RunTFIDF(obj) |> RunSVD(n = cfg$atac$n_lsi_components, reduction.name = "lsi")
obj <- RunUMAP(obj, reduction = "lsi", dims = ad, reduction.name = "umap.atac.unint", reduction.key = "atacUMAPu_", verbose = FALSE)
obj <- harmony::RunHarmony(obj, group.by.vars = gb, reduction.use = "lsi", dims.use = ad, assay.use = "peaks",
                           theta = cfg$integration$atac_theta, reduction.save = "harmony.atac", project.dim = FALSE, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony.atac", dims = seq_along(ad), reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE)

# ---- WNN ----
obj <- FindMultiModalNeighbors(obj, reduction.list = list("harmony.rna", "harmony.atac"),
                               dims.list = list(rd, seq_along(ad)), k.nn = 20, verbose = FALSE)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", verbose = FALSE)
for (r in cfg$clustering$resolutions) {
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = r, algorithm = cfg$clustering$algorithm,
                      random.seed = cfg$clustering$seed, verbose = FALSE, cluster.name = paste0("wsnn_res.", r))
}
obj$seurat_clusters <- obj[[paste0("wsnn_res.", cfg$clustering$default_resolution)]][, 1]
Idents(obj) <- "seurat_clusters"
DefaultAssay(obj) <- "RNA"
if (cfg$rna$normalisation == "SCT") DefaultAssay(obj) <- "SCT"

# ---- diagnostics ----
pl <- list(
  DimPlot(obj, reduction = "umap.rna.unint", group.by = gb) + ggtitle("RNA before"),
  DimPlot(obj, reduction = "umap.rna", group.by = gb) + ggtitle("RNA Harmony"),
  DimPlot(obj, reduction = "umap.atac.unint", group.by = gb) + ggtitle("ATAC before"),
  DimPlot(obj, reduction = "umap.atac", group.by = gb) + ggtitle("ATAC Harmony"))
save_plot(wrap_plots(pl, ncol = 2), file.path(outd, "integration_before_after.png"), 12, 10)
p <- DimPlot(obj, reduction = "umap.wnn", group.by = gb) | DimPlot(obj, reduction = "umap.wnn", label = TRUE)
save_plot(p, file.path(outd, "umap_wnn.png"), 14, 6)
write.csv(table(obj$seurat_clusters, obj[[gb]][, 1]), file.path(outd, "cluster_by_sample.csv"))
saveRDS(obj, file.path(cfg$data_path, "mouse_multiome.integrated.rds"))
session_log(outd)
