# DLPFC 12 sections: WM (target) versus all other labelled spots (rest) in the
# composition-adjusted mode of compare(), 30 gene sets of
# benchmark/analysis/DLPFC/data/DLPFC12_gene_sets.gmt (the sets of Figure X).
# Inputs of each section, in benchmark/analysis/DLPFC/data/:
#   * raw UMI counts of <section>/<section>_filtered_feature_bc_matrix.h5;
#     gene-set symbols matched through <section>/genes.csv (first Ensembl ID
#     of a repeated symbol)
#   * spots with an author layer label that RCTD fitted
#   * RCTD full-mode proportions of the 7 LIBD types
#     (RCTD_LIBD/<section>/<section>_spot_celltype_proportions.csv) and the
#     reference profiles RCTD used (reference_profiles.rds)
# Each section is fitted on its own; results in
# benchmark/analysis/DLPFC/adjusted/<section>/. A finished section is skipped.
# Run from this folder.
# Usage: Rscript run_wm_adjusted.R [SECTION ...]   (default: all 12)
suppressPackageStartupMessages({library(Matrix); library(PISTE)})

work_dir <- "../.."
data_root <- file.path(work_dir, "data")
out_root <- file.path(work_dir, "adjusted")
gmt_file <- file.path(data_root, "DLPFC12_gene_sets.gmt")
types <- c("Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC")
sections <- c("151507", "151508", "151509", "151510", "151669", "151670",
              "151671", "151672", "151673", "151674", "151675", "151676")
args <- commandArgs(TRUE)
if (length(args)) {stopifnot(all(args %in% sections)); sections <- args}

profiles <- readRDS(file.path(data_root, "reference_profiles.rds"))$profiles   # genes x cell types

read_counts <- function(file) {
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  counts <- sparseMatrix(i = h5[["matrix/indices"]][] + 1, p = h5[["matrix/indptr"]][],
                         x = as.numeric(h5[["matrix/data"]][]), dims = h5[["matrix/shape"]][])
  dimnames(counts) <- list(h5[["matrix/features/id"]][], h5[["matrix/barcodes"]][])
  counts
}

for (section in sections) {
  out <- file.path(out_root, section)
  if (file.exists(file.path(out, "overall_scores.csv"))) next
  sec_dir <- file.path(data_root, section)

  # Raw counts, spots x genes, with gene symbols as column names.
  counts <- read_counts(file.path(sec_dir, paste0(section, "_filtered_feature_bc_matrix.h5")))
  genes_csv <- read.csv(file.path(sec_dir, "genes.csv"), stringsAsFactors = FALSE)
  genes_csv <- genes_csv[!duplicated(genes_csv$symbol), ]
  counts <- t(counts[genes_csv$gene_id, ])
  colnames(counts) <- genes_csv$symbol

  # Labelled spots that RCTD fitted, and their cell-type proportions.
  rctd <- read.csv(file.path(data_root, "RCTD_LIBD", section,
                             paste0(section, "_spot_celltype_proportions.csv")),
                   stringsAsFactors = FALSE)
  spots <- rctd[!is.na(rctd$layer) & nzchar(rctd$layer) & rctd$fit_status == "fitted", ]
  proportions <- as.matrix(spots[, types])
  rownames(proportions) <- spots$spot
  region <- ifelse(spots$layer == "WM", "WM", "rest")

  # Reference profiles, cell types x genes, with gene symbols.
  ids <- intersect(genes_csv$gene_id, rownames(profiles))
  reference <- t(profiles[ids, types, drop = FALSE])
  colnames(reference) <- genes_csv$symbol[match(ids, genes_csv$gene_id)]

  set.seed(42)
  compare(counts[spots$spot, ], region, gmt_file, mode = "adjusted",
          proportions = proportions, reference = reference,
          baseline = "rest", out = out)
}
