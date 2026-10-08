# Target-vs-rest PathwayLGM fits for all pathways.

# ---- config ---------------------------------------------------------------------
work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")
set.seed(42)  

# ---- load -----------------------------------------------------------------------
library(PathwayLGM)
library(fmesher)

x <- readRDS(file.path(data_dir, "raw_data.rds"))
residual <- readRDS(file.path(data_dir, "processed", "pearson_residuals.rds"))

mesh <- fm_mesh_2d(loc = x$xy, max.edge = .25, cutoff = .08, offset = c(.1, .2))
X <- matrix(1, nrow(x$xy), 1, dimnames = list(rownames(x$xy), "intercept"))
out <- file.path(work_dir, "PathwayLGM")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- run ------------------------------------------------------------------------
for (p in names(x$pathways)) {
  Y <- t(residual[x$pathways[[p]], , drop = FALSE])
  region <- ifelse(x$oracle_target[, p], "target", "rest")
  fit <- screen_joint_fit(Y, x$xy, region, mesh, .3, family = "gaussian",
                          X = X, cfg = list(threads = "2:1", pi = TRUE,
                                            baseline_region = "rest"))
  screen_write(screen_summarize(fit), file.path(out, p))
}
