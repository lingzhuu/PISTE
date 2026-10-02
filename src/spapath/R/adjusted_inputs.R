# =============================================================================
# adjusted_inputs.R
#
# Inputs of the composition-adjusted mode of screen_joint_fit():
#
#   composition_offset()         log expression predicted from cell-type
#                                proportions and a reference profile
#   estimate_dispersion_scale()  transcriptome-wide relative NB dispersions
#
# Typical workflow
#   offset <- composition_offset(proportions, reference)
#   scales <- estimate_dispersion_scale(counts, library_size)
# =============================================================================


#' Log expression predicted from cell-type composition
#'
#' For spot i and gene g, the predicted expression is the proportion-weighted
#' average of the cell-type reference profiles,
#' m[i, g] = sum_t p[i, t] nu[t, g], and the offset is log m[i, g] with no
#' pseudocount. Proportions are rescaled to sum to one in every spot.
#'
#' @param proportions Non-negative spots x cell types matrix of cell-type
#'   proportions, e.g. normalised RCTD weights, with spot rownames.
#' @param reference Non-negative cell types x genes matrix of reference mean
#'   expression on a linear scale, e.g. counts per 10,000. Its rownames must
#'   match `colnames(proportions)`.
#' @return Spots x genes matrix of log m[i, g]. Genes with zero predicted
#'   expression in some spot get `-Inf` there and must be removed before
#'   fitting.
#' @export
composition_offset <- function(proportions, reference) {
  proportions <- as.matrix(proportions)
  reference <- as.matrix(reference)
  if (!is.numeric(proportions) || !is.numeric(reference) ||
      any(!is.finite(proportions)) || any(proportions < 0) ||
      any(!is.finite(reference)) || any(reference < 0)) {
    stop("proportions and reference must be finite non-negative numeric matrices.")
  }
  if (is.null(rownames(reference)) ||
      !setequal(rownames(reference), colnames(proportions))) {
    stop("rownames(reference) must be the cell types in colnames(proportions).")
  }

  proportions <- proportions[, rownames(reference), drop = FALSE]
  totals <- rowSums(proportions)
  if (any(totals <= 0)) stop("Every spot needs a positive total proportion.")

  log((proportions / totals) %*% reference)
}


#' Relative negative binomial dispersion of every gene
#'
#' Estimates gene-wise NB dispersions alpha_g transcriptome-wide with DESeq2
#' (intercept-only design, trend fit and MAP shrinkage), then returns the
#' relative dispersion scale S_g = exp(mean(log alpha)) / alpha_g, so that
#' the scales have geometric mean one. Pass `dispersion_scale[genes]` to
#' [screen_joint_fit()]; never re-estimate on a single gene set.
#'
#' Gene-wise estimation and MAP shrinkage run in blocks of `block_size`
#' genes to bound memory, while the dispersion trend and its prior variance
#' are fitted once on all eligible genes. This gives the same estimates as
#' an unblocked DESeq2 run.
#'
#' @param counts Raw counts, spots x genes (dense or sparse), with unique
#'   spot rownames and gene colnames. Use all genes, not only pathway genes.
#' @param library_size Length-N library sizes, e.g. total UMI per spot.
#' @param min_total,min_detected A gene is eligible when its total count is
#'   at least `min_total` and it is detected in at least `min_detected` spots.
#'   At least 30 genes must be eligible.
#' @param fit_type DESeq2 dispersion trend type.
#' @param block_size Number of genes per block.
#' @return Data frame with one row per gene: `gene`, `total`, `detected`,
#'   `eligible`, `dispersion` (alpha_g) and `dispersion_scale` (S_g; NA for
#'   ineligible genes).
#' @export
estimate_dispersion_scale <- function(counts, library_size, min_total = 10L,
                                      min_detected = 3L,
                                      fit_type = "parametric",
                                      block_size = 256L) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("estimate_dispersion_scale() requires the Bioconductor package DESeq2.")
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("counts needs unique spot rownames and gene colnames.")
  }
  if (length(library_size) != nrow(counts) || any(!is.finite(library_size)) ||
      any(library_size <= 0)) {
    stop("library_size must be positive and finite, one value per spot.")
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values)) || any(values < 0) || any(values != round(values)) ||
      any(values > .Machine$integer.max)) {
    stop("counts must hold non-negative integer counts.")
  }
  rm(values)

  total <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  eligible <- total >= min_total & detected >= min_detected
  if (sum(eligible) < 30L) {
    stop("Need at least 30 eligible genes; estimate on the whole transcriptome, not one gene set.")
  }

  # DESeq2 works on genes x samples.
  gene_counts <- t(as.matrix(counts[, eligible, drop = FALSE]))
  storage.mode(gene_counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(
    gene_counts, data.frame(row.names = rownames(counts)), design = ~1)
  DESeq2::sizeFactors(dds) <- library_size / exp(mean(log(library_size)))
  rm(gene_counts)
  gc(FALSE)

  dds <- dispersion_blocked_deseq2(dds, fit_type, block_size)
  alpha <- DESeq2::dispersions(dds)
  if (any(!is.finite(alpha) | alpha <= 0)) {
    stop("DESeq2 returned invalid dispersions; inspect the counts.")
  }

  result <- data.frame(gene = colnames(counts), total = total,
                       detected = detected, eligible = eligible,
                       dispersion = NA_real_, dispersion_scale = NA_real_,
                       row.names = NULL)
  result$dispersion[eligible] <- alpha
  result$dispersion_scale[eligible] <- exp(mean(log(alpha)) - log(alpha))
  result
}


#' DESeq2 dispersions with gene-wise and MAP steps in blocks of genes
#'
#' Step 1 fits gene-wise dispersions block by block, keeping each block's
#' fitted means on disk. Step 2 fits one global trend and prior variance on
#' all genes. Step 3 applies MAP shrinkage block by block with that global
#' trend and prior; a separate trend per block is never fitted.
#' @noRd
dispersion_blocked_deseq2 <- function(dds, fit_type, block_size) {
  blocks <- split(seq_len(nrow(dds)), ceiling(seq_len(nrow(dds)) / block_size))
  mu_dir <- tempfile("dispersion_mu_")
  dir.create(mu_dir)
  on.exit(unlink(mu_dir, recursive = TRUE), add = TRUE)
  mu_file <- function(b) file.path(mu_dir, sprintf("block_%04d.rds", b))

  # Step 1: gene-wise estimates.
  gene_estimates <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    block <- DESeq2::estimateDispersionsGeneEst(dds[blocks[[b]], ], quiet = TRUE)
    gene_estimates[[b]] <- S4Vectors::mcols(block)
    saveRDS(SummarizedExperiment::assay(block, "mu"), mu_file(b), compress = FALSE)
    rm(block)
    gc(FALSE)
  }
  S4Vectors::mcols(dds) <- do.call(rbind, gene_estimates)

  # Step 2: global trend and prior variance.
  dds <- DESeq2::estimateDispersionsFit(dds, fitType = fit_type, quiet = TRUE)
  prior_variance <- DESeq2::estimateDispersionsPriorVar(dds)
  trend <- DESeq2::dispersionFunction(dds)
  attr(trend, "dispPriorVar") <- prior_variance
  # The public dispersionFunction<- re-estimates varLogDispEsts on the object
  # it is given; set the slot directly to keep the global value in each block.
  methods::slot(dds, "dispersionFunction") <- trend

  # Step 3: MAP shrinkage.
  shrunken <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    block <- dds[blocks[[b]], ]
    SummarizedExperiment::assay(block, "mu") <- readRDS(mu_file(b))
    methods::slot(block, "dispersionFunction") <- trend
    block <- DESeq2::estimateDispersionsMAP(block, dispPriorVar = prior_variance,
                                            quiet = TRUE)
    shrunken[[b]] <- S4Vectors::mcols(block)
    rm(block)
    gc(FALSE)
  }
  S4Vectors::mcols(dds) <- do.call(rbind, shrunken)
  dds
}
