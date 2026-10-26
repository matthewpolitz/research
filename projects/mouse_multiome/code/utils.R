# Shared helpers for the mouse_multiome pipeline
suppressPackageStartupMessages({ library(yaml); library(glue) })

load_config <- function(path = Sys.getenv("MM_CONFIG", "config/config.yml")) {
  cfg <- yaml::read_yaml(path)
  cfg$data_path    <- file.path(cfg$project_root, cfg$data_dir)
  cfg$results_path <- file.path(cfg$project_root, cfg$results_dir)
  dir.create(cfg$data_path, FALSE, TRUE); dir.create(cfg$results_path, FALSE, TRUE)
  cfg
}

outs_dir <- function(cfg, sample) {
  glue::glue(cfg$outs_template, cellranger_dir = cfg$cellranger_dir, sample = sample)
}

res_dir <- function(cfg, step) { d <- file.path(cfg$results_path, step); dir.create(d, FALSE, TRUE); d }

parse_dims <- function(x) eval(parse(text = x))

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# MAD-based bounds
mad_bounds <- function(x, nmads = 3, log = TRUE) {
  v <- if (log) log10(x + 1) else x
  m <- median(v); s <- mad(v)
  b <- c(m - nmads * s, m + nmads * s)
  if (log) b <- 10^b - 1
  b
}

save_plot <- function(p, file, w = 8, h = 6) {
  ggplot2::ggsave(file, p, width = w, height = h, dpi = 200)
  ggplot2::ggsave(sub("\\.png$", ".pdf", file), p, width = w, height = h)
}

session_log <- function(dir) {
  writeLines(capture.output(sessionInfo()), file.path(dir, "sessionInfo.txt"))
  writeLines(Sys.getenv("MM_COMMIT", "unknown"), file.path(dir, "deployed_commit.txt"))
}
