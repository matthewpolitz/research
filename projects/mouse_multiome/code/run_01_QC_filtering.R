#!/usr/bin/env Rscript
# Render 01_QC_filtering.Rmd for every sample in the sample sheet (each in a
# fresh R session), then combine the per-sample QC tables.
#
# Usage, from the project root:
#   Rscript code/run_01_QC_filtering.R                 # all samples
#   Rscript code/run_01_QC_filtering.R Run1_Sample1    # selected samples

rmd          <- "code/01_QC_filtering.Rmd"
sample_sheet <- "config/samples_QC.csv"
summary_dir  <- "data/qc/summary"

if (!file.exists(rmd)) stop("Run this script from the project root")

# Outside RStudio, fall back to the pandoc bundled with RStudio (Desktop on macOS,
# Server on Linux)
if (!rmarkdown::pandoc_available()) {
  rstudio_pandoc <- Sys.glob(c(
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/*/pandoc",
    "/usr/lib/rstudio-server/bin/quarto/bin/tools/*/pandoc",
    "/usr/lib/rstudio-server/bin/quarto/bin/tools/pandoc",
    "/usr/lib/rstudio-server/bin/pandoc/pandoc"
  ))
  arch_match <- grepl(R.version$arch, rstudio_pandoc)
  rstudio_pandoc <- c(rstudio_pandoc[arch_match], rstudio_pandoc[!arch_match])
  if (length(rstudio_pandoc)) Sys.setenv(RSTUDIO_PANDOC = dirname(rstudio_pandoc[1]))
}

samples <- read.csv(sample_sheet)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  unknown <- setdiff(args, samples$sample)
  if (length(unknown)) stop("Not in ", sample_sheet, ": ", paste(unknown, collapse = ", "))
  samples <- samples[samples$sample %in% args, ]
}

# --- Render each sample ---------------------------------------------------
ok <- vapply(samples$sample, function(s) {
  out_dir <- file.path(getwd(), samples$out_dir[samples$sample == s])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== ", s, " ===")
  tryCatch({
    xfun::Rscript_call(rmarkdown::render, list(
      input       = rmd,
      params      = list(sample = s, sample_sheet = sample_sheet),
      output_file = paste0("01_QC_filtering_", s, ".html"),
      output_dir  = out_dir,
      quiet       = TRUE
    ))
    TRUE
  }, error = function(e) {
    message("FAILED: ", s, "\n", conditionMessage(e))
    FALSE
  })
}, logical(1))

# --- Combine per-sample tables --------------------------------------------
done <- samples[ok, ]
if (nrow(done)) {
  read_all <- function(suffix) {
    do.call(rbind, lapply(seq_len(nrow(done)), function(i) {
      read.csv(file.path(done$out_dir[i], paste0(done$sample[i], "_", suffix)),
               check.names = FALSE)
    }))
  }
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  cells_removed <- read_all("qc_cells_removed.csv")
  write.csv(cells_removed, file.path(summary_dir, "qc_cells_removed_all_samples.csv"), row.names = FALSE)
  write.csv(read_all("qc_thresholds.csv"),
            file.path(summary_dir, "qc_thresholds_all_samples.csv"), row.names = FALSE)

  # Cells remaining after each step, one column per sample
  remaining <- reshape(cells_removed[, c("Sample", "Step", "Remaining")],
                       idvar = "Step", timevar = "Sample", direction = "wide")
  names(remaining) <- sub("^Remaining\\.", "", names(remaining))
  write.csv(remaining, file.path(summary_dir, "qc_cells_remaining_by_sample.csv"), row.names = FALSE)
  print(remaining, row.names = FALSE)
}

message("\nSucceeded: ", paste(samples$sample[ok], collapse = ", "))
if (any(!ok)) {
  message("Failed: ", paste(samples$sample[!ok], collapse = ", "))
  quit(status = 1)
}
