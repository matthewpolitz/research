# 02: merged peak set + per-sample normalisation. Input: step 01 outputs (SoupX-corrected RNA, ATAC on cellranger peaks)
# 1) optional MACS2 peak calling per sample, merge to a common peak set, requantify
# 2) RNA: SCTransform or LogNormalize; ATAC: TF-IDF
suppressPackageStartupMessages({ library(rtracklayer); library(Seurat); library(Signac); library(GenomicRanges); library(dplyr) })
source("code/utils.R"); cfg <- load_config()
outd <- res_dir(cfg, "02_normalise")
objs <- setNames(lapply(cfg$samples, function(s) readRDS(file.path(cfg$project_root, glue::glue(cfg$qc_rds_template, sample = s)))), cfg$samples)
annotation <- Annotation(objs[[1]][["ATAC"]])   # GENCODE vM38, chromosome style already matched to fragments

# ---- common peak set ----
if (isTRUE(cfg$atac$call_peaks)) {
  macs2 <- if (nzchar(cfg$atac$macs2_path)) cfg$atac$macs2_path else Sys.which("macs2")
  if (!nzchar(macs2)) stop("macs2 not found; install via conda env or set atac.macs2_path, or set atac.call_peaks: false")
  peaks_list <- lapply(names(objs), function(s) {
    log_msg("MACS2 ", s)
    DefaultAssay(objs[[s]]) <- "ATAC"
    CallPeaks(objs[[s]], macs2.path = macs2, outdir = file.path(cfg$data_path, "macs2", s), cleanup = FALSE)
  })
  combined <- reduce(do.call(c, peaks_list))
  combined <- keepStandardChromosomes(combined, pruning.mode = "coarse")
  combined <- combined[width(combined) < cfg$atac$peak_max_width & width(combined) > cfg$atac$peak_min_width]
  # drop blacklist
  if (nzchar(cfg$genome$blacklist_bed)) combined <- subsetByOverlaps(combined, rtracklayer::import(cfg$genome$blacklist_bed), invert = TRUE)
} else {
  # use cellranger-called peaks (union over samples)
  # union of the cellranger-arc peaks already stored in each sample's ATAC assay
  combined <- reduce(do.call(c, lapply(objs, function(o) granges(o[["ATAC"]]))))
}
log_msg("Combined peak set: ", length(combined))
saveRDS(combined, file.path(cfg$data_path, "combined_peaks.rds"))
write.table(as.data.frame(combined)[, 1:3], file.path(outd, "combined_peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- requantify each sample on common peaks and normalise ----
for (s in names(objs)) {
  log_msg("Requantify + normalise ", s)
  o <- objs[[s]]
  frag <- Fragments(o[["ATAC"]])
  cnt <- FeatureMatrix(fragments = frag, features = combined, cells = colnames(o))
  o[["peaks"]] <- CreateChromatinAssay(cnt, fragments = frag, annotation = annotation, min.cells = 0)
  o[["ATAC"]] <- NULL; o[["RNA_raw"]] <- NULL   # keep uncorrected counts in the QC rds only

  # RNA
  DefaultAssay(o) <- "RNA"
  if (cfg$rna$normalisation == "SCT") {
    o <- SCTransform(o, vars.to.regress = unlist(cfg$rna$regress_out), variable.features.n = cfg$rna$n_variable_features, verbose = FALSE)
  } else {
    o <- NormalizeData(o) |> FindVariableFeatures(nfeatures = cfg$rna$n_variable_features) |> ScaleData(vars.to.regress = unlist(cfg$rna$regress_out))
  }
  o <- RunPCA(o, npcs = cfg$rna$npcs, verbose = FALSE)

  # ATAC
  DefaultAssay(o) <- "peaks"
  o <- FindTopFeatures(o, min.cutoff = cfg$atac$min_cells_feature)
  o <- RunTFIDF(o) |> RunSVD(n = cfg$atac$n_lsi_components)

  rd <- parse_dims(cfg$rna$dims); ad <- parse_dims(cfg$atac$lsi_dims)
  o <- RunUMAP(o, reduction = "pca", dims = rd, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE)
  o <- RunUMAP(o, reduction = "lsi", dims = ad, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE)
  p <- (DimPlot(o, reduction = "umap.rna") + ggtitle(paste(s, "RNA"))) | (DimPlot(o, reduction = "umap.atac") + ggtitle(paste(s, "ATAC")))
  save_plot(p, file.path(outd, paste0(s, "_umaps_unintegrated.png")), 12, 5)
  p <- DepthCor(o, reduction = "lsi"); save_plot(p, file.path(outd, paste0(s, "_lsi_depth_correlation.png")), 6, 4)
  saveRDS(o, file.path(cfg$data_path, paste0(s, ".norm.rds")))
}
session_log(outd)
