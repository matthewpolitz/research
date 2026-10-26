# Verify cellranger-arc outputs exist for each sample before anything else runs.
source("code/utils.R"); cfg <- load_config()
need <- c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5", "atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi",
          "per_barcode_metrics.csv", "atac_peaks.bed")
ok <- TRUE
for (s in cfg$samples) {
  d <- outs_dir(cfg, s); cat("\n==", s, ":", d, "\n")
  for (f in need) { e <- file.exists(file.path(d, f)); cat(sprintf("  %-40s %s\n", f, if (e) "OK" else "MISSING")); ok <- ok && e }
}
cat("\nAvailable files in first sample outs:\n"); print(list.files(outs_dir(cfg, cfg$samples[1])))
if (!ok) stop("Some inputs missing - fix outs_template / cellranger_dir in config.yml")
