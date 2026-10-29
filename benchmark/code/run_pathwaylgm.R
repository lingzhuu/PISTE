# Target-vs-rest PathwayLGM fits for all pathways.

# ---- config ---------------------------------------------------------------------
work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")

# ---- load -----------------------------------------------------------------------
library(PathwayLGM)

x <- readRDS(file.path(data_dir, "raw_data.rds"))
residual <- readRDS(file.path(data_dir, "processed", "pearson_residuals.rds"))
out <- file.path(work_dir, "PathwayLGM")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- run ------------------------------------------------------------------------
for (p in names(x$pathways)) {
  Y <- t(residual[x$pathways[[p]], , drop = FALSE])
  region <- ifelse(x$oracle_target[, p], "target", "rest")
  fit <- observed_fit(Y, region, cfg = list(baseline_region = "rest"))
  screen_write(screen_summarize(fit), file.path(out, p))
}
