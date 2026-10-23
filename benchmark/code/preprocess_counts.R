# Preprocessing shared by every analysis in the paper (DLPFC, semi-synthetic benchmark,
# mouse brain, breast cancer). Applied to the Space Ranger filtered (in-tissue) matrix
# before any normalization; see README, Methods > Preprocessing.
#
#   1. Spots: keep the spots named in `spots` (those with a region/layer annotation);
#      no further spot-level filtering (Totty et al., 2025; Maynard et al., 2021).
#   2. Genes: drop mitochondrial genes (names starting MT- or mt-), whose share of the
#      library differs systematically between regions and would shift the normalized
#      values of every other gene in the same direction; then drop genes detected in
#      fewer than 5 of the retained spots (SCTransform's own minimum, min_cells = 5).
#      Duplicate gene symbols are summed when the matrix is read from the 10x h5 file.
#
# preprocess_counts(counts, spots) returns list(counts = genes x spots sparse matrix,
# manifest = data.frame of what was removed). As a script:
#   Rscript preprocess_counts.R <filtered_feature_bc_matrix.h5> <spots.csv> <out_prefix>
# writes <out_prefix>_counts.mtx(.gz), <out_prefix>_genes.txt, <out_prefix>_spots.txt and
# <out_prefix>_manifest.csv; spots.csv needs a column `spot` (barcode).
suppressPackageStartupMessages(library(Matrix))

preprocess_counts <- function(counts, spots = colnames(counts), drop_mito = TRUE, min_spots = 5,
                              mito_pattern = "^(MT|mt)-") {
  spots <- intersect(spots, colnames(counts))
  counts <- counts[, spots, drop = FALSE]
  is_mito <- grepl(mito_pattern, rownames(counts))
  keep <- rep(TRUE, nrow(counts))
  if (drop_mito) keep <- keep & !is_mito
  sparse <- Matrix::rowSums(counts > 0) < min_spots
  keep <- keep & !sparse
  manifest <- data.frame(
    spots_kept = length(spots), genes_in = nrow(counts),
    mito_genes_removed = if (drop_mito) sum(is_mito) else 0L,
    min_spots = min_spots,
    sparse_genes_removed = sum(sparse & !(drop_mito & is_mito)),
    genes_kept = sum(keep),
    mito_fraction_of_counts = sum(counts[is_mito, ]) / sum(counts))
  list(counts = counts[keep, , drop = FALSE], manifest = manifest)
}

read_10x_h5 <- function(file) {
  h5 <- hdf5r::H5File$new(file, mode = "r"); on.exit(h5$close_all())
  m <- sparseMatrix(i = h5[["matrix/indices"]][] + 1, p = h5[["matrix/indptr"]][],
                    x = as.numeric(h5[["matrix/data"]][]), dims = h5[["matrix/shape"]][])
  dimnames(m) <- list(h5[["matrix/features/name"]][], h5[["matrix/barcodes"]][])
  # duplicate symbols (24 in the DLPFC data) are summed, as in the paper's SCT input preparation
  if (anyDuplicated(rownames(m))) {
    f <- factor(rownames(m), levels = unique(rownames(m)))
    m <- as(Matrix::sparseMatrix(i = as.integer(f)[m@i + 1L], p = m@p, x = m@x, dims = c(nlevels(f), ncol(m)),
                                 dimnames = list(levels(f), colnames(m))), "CsparseMatrix")
  }
  m
}

if (sys.nframe() == 0 && length(commandArgs(trailingOnly = TRUE)) == 3) {
  a <- commandArgs(trailingOnly = TRUE)
  res <- preprocess_counts(read_10x_h5(a[1]), read.csv(a[2])$spot)
  writeMM(res$counts, paste0(a[3], "_counts.mtx"))
  writeLines(rownames(res$counts), paste0(a[3], "_genes.txt"))
  writeLines(colnames(res$counts), paste0(a[3], "_spots.txt"))
  write.csv(res$manifest, paste0(a[3], "_manifest.csv"), row.names = FALSE)
  print(res$manifest)
}
