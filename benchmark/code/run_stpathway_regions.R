# Run the STpathway predefined-area test on every pathway of one simulation
# dataset and save the combined score table for downstream use:
#   stpathway_region_scores.csv one row per pathway
#
# Per pathway, the SCT Pearson residuals of its genes are split into spots
# inside vs outside the oracle target area and compared with
# energy::dcov.test: genes are the observations and the spots of one area
# form one sample.

# ---- config ---------------------------------------------------------------------
library(energy)

work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")

# ---- load -----------------------------------------------------------------------
x <- readRDS(file.path(data_dir, "raw_data.rds"))
residual <- readRDS(file.path(data_dir, "processed", "pearson_residuals.rds"))

# ---- run ------------------------------------------------------------------------
rows <- list()
for (p in names(x$pathways)) {
  genes <- x$pathways[[p]]
  mask <- as.logical(x$oracle_target[, p])
  set.seed(42)
  res <- dcov.test(residual[genes, mask, drop = FALSE],
                   residual[genes, !mask, drop = FALSE], R = 199)
  rows[[p]] <- data.frame(ID = p, n_gene = length(genes),
    pval = res$p.value, dcov = unname(res$estimates[1]),
    dcor = unname(res$estimates[2]),
    target_spots = sum(mask), rest_spots = sum(!mask))
  cat(p, "p =", rows[[p]]$pval, "dcor =", rows[[p]]$dcor, "\n")
}

# ---- save -----------------------------------------------------------------------
write.csv(do.call(rbind, rows),
          file.path(data_dir, "stpathway_region_scores.csv"),
          row.names = FALSE)
