# fgsea with genes ranked by Wilcoxon AUC (presto::wilcoxauc), all simulated_v4 seeds.
# Data = log1p of the SCT corrected counts, each
# pathway's oracle target vs the rest. Ranking = AUC - 0.5, ties broken by 1e-6 * logFC.
# Saves the full presto table, the ranking and the fgsea row per pathway; finished pathways
# are skipped, so a rerun resumes.
# Usage: Rscript run_findmarkers_fgsea.R <simulated_v4 root> <out_dir>
suppressPackageStartupMessages({ library(Matrix); library(presto); library(fgsea); library(BiocParallel) })
args <- commandArgs(TRUE)
root <- args[1]; out_root <- args[2]

for (s in 2020:2026) {
  seed_dir <- file.path(root, paste0("seed", s)); out_dir <- file.path(out_root, paste0("seed", s))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  d <- readRDS(file.path(seed_dir, "dataset.rds"))
  cc <- readRDS(file.path(seed_dir, "SCTransform_v2", "corrected_counts.rds"))
  stopifnot(identical(colnames(cc), rownames(d$oracle_target)))
  x <- log1p(cc)
  status_file <- file.path(out_dir, "status.csv")
  for (p in colnames(d$oracle_target)) {
    pdir <- file.path(out_dir, p)
    if (file.exists(file.path(pdir, "fgsea.rds"))) next
    dir.create(pdir, showWarnings = FALSE)
    t0 <- Sys.time()
    grp <- ifelse(d$oracle_target[, p], "target", "rest")
    w <- wilcoxauc(x, grp)
    w <- w[w$group == "target", ]
    t1 <- Sys.time()
    stat <- w$auc - .5 + 1e-6 * w$logFC
    names(stat) <- w$feature
    stat <- sort(stat, decreasing = TRUE)
    set.seed(s * 100 + match(p, colnames(d$oracle_target)))
    fg <- fgsea(pathways = setNames(list(d$pathways[[p]]), p), stats = stat, eps = 0,
                minSize = 1, maxSize = Inf, BPPARAM = SerialParam())
    t2 <- Sys.time()
    saveRDS(w, file.path(pdir, "wilcoxauc.rds"))
    saveRDS(stat, file.path(pdir, "ranking.rds"))
    saveRDS(fg, file.path(pdir, "fgsea.rds"))
    st <- data.frame(seed = s, pathway = p, n_target = sum(d$oracle_target[, p]), n_genes_ranked = length(stat),
                     auc_seconds = as.numeric(t1 - t0, units = "secs"),
                     fgsea_seconds = as.numeric(t2 - t1, units = "secs"), finished = format(t2, "%Y-%m-%dT%H:%M:%S"))
    write.table(st, status_file, sep = ",", row.names = FALSE, append = file.exists(status_file),
                col.names = !file.exists(status_file))
  }
  cat("seed", s, "done\n")
}
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
