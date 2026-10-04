# =============================================================================
# run_spapath.R
#
# One-call workflow: run SpaPath on every gene set of a panel.
#
#   res <- run_spapath(counts, xy, region, gene_sets)            # observed
#   res <- run_spapath(counts, xy, region, gene_sets,
#                      mode = "adjusted",
#                      proportions = omega, reference = nu)       # adjusted
#
# Steps
#   align spots -> normalise with SCTransform (observed) or build the
#   dispersion scales and composition offset (adjusted) -> scale coordinates
#   and build the mesh -> fit, summarise (and write) each gene set ->
#   combined summary tables
#
# The step-by-step functions it calls (build_mesh(), estimate_dispersion_scale(),
# composition_offset(), screen_joint_fit(), screen_summarize(), screen_write())
# remain available for custom workflows.
# =============================================================================


#' Run SpaPath on a panel of gene sets
#'
#' Runs the whole workflow in one call: normalisation (observed mode) or the
#' composition inputs (adjusted mode), the mesh, and one fit per gene set.
#' Coordinates may be in any unit; they are rescaled to the unit square.
#'
#' @param counts Spots x genes matrix (dense or sparse) of raw counts, with
#'   spot rownames and gene colnames. In observed mode with
#'   `normalized = TRUE`, normalised expression instead.
#' @param xy Spots x 2 matrix of spot coordinates, in the spot order of
#'   `counts`.
#' @param region Region label of each spot, in the spot order of `counts`.
#' @param gene_sets Named list of gene vectors, or the path of a GMT file.
#' @param mode `"observed"` (Gaussian model of SCTransform-normalised
#'   expression) or `"adjusted"` (negative binomial model of raw counts
#'   adjusted for cell-type composition).
#' @param out Optional output folder. When given, each gene set is written to
#'   its own subfolder with [screen_write()], together with
#'   `overall_scores.csv`, `pairwise_scores.csv` and `gene_sets.csv`.
#' @param normalized Observed mode: `counts` is already normalised, so
#'   SCTransform is skipped.
#' @param proportions Adjusted mode: spots x cell types matrix of cell-type
#'   proportions, e.g. normalised RCTD weights, with spot rownames.
#' @param reference Adjusted mode: cell types x genes matrix of reference
#'   mean expression on a linear scale, e.g. from [reference_profile()].
#' @param baseline Baseline region of the model (default: the first label).
#' @param h Spatial length scale in unit-square units.
#' @param min_genes Gene sets with fewer usable genes are skipped.
#' @param mesh_args List of arguments passed to [build_mesh()].
#' @param cfg List of options passed to [screen_joint_fit()]. Unless set
#'   here, sets of 150 genes or more use one thread and `low_memory`.
#' @param verbose Print progress messages.
#' @return Invisibly, a list with
#'   * `overall`, `pairwise`: the π tables of all gene sets, with columns
#'     `gene_set` and `genes` added;
#'   * `gene_sets`: one row per gene set with the genes used and whether it
#'     was fitted, skipped or failed;
#'   * `results`: the [screen_summarize()] output of each fitted gene set.
#' @export
run_spapath <- function(counts, xy, region, gene_sets,
                        mode = c("observed", "adjusted"), out = NULL,
                        normalized = FALSE, proportions = NULL,
                        reference = NULL, baseline = NULL, h = 0.3,
                        min_genes = 5L, mesh_args = list(), cfg = list(),
                        verbose = TRUE) {
  mode <- match.arg(mode)
  log_step <- function(...) if (verbose) message(format(Sys.time(), "[%H:%M:%S] "), ...)

  # --- Inputs ----------------------------------------------------------------
  if (is.character(gene_sets)) gene_sets <- read_gmt(gene_sets)
  stopifnot(is.list(gene_sets), length(gene_sets) > 0, !is.null(names(gene_sets)))
  xy <- as.matrix(xy)
  region <- as.character(region)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts needs spot rownames and gene colnames.")
  }
  if (nrow(xy) != nrow(counts) || length(region) != nrow(counts)) {
    stop("xy and region must have one entry per row of counts.")
  }
  if (anyNA(region)) stop("region has missing labels; remove those spots first.")
  if (mode == "adjusted") {
    if (is.null(proportions) || is.null(reference)) {
      stop("The adjusted mode needs proportions and reference.")
    }
    if (normalized) stop("normalized applies to the observed mode only.")
  }

  # --- Model inputs ------------------------------------------------------------
  pathway_genes <- unique(unlist(gene_sets))
  if (mode == "observed") {
    expression <- if (normalized) {
      as.matrix(counts[, intersect(pathway_genes, colnames(counts)), drop = FALSE])
    } else {
      log_step("Normalising with SCTransform")
      sctransform_residuals(counts, pathway_genes)
    }
    available <- colnames(expression)
  } else {
    log_step("Estimating gene dispersions and the composition offset")
    adjusted <- adjusted_inputs(counts, proportions, reference)
    available <- adjusted$usable
  }
  log_step(length(intersect(pathway_genes, available)), " of ",
           length(pathway_genes), " pathway genes are usable")

  # --- Mesh --------------------------------------------------------------------
  xy <- unit_square(xy)
  mesh <- do.call(build_mesh, c(list(xy), mesh_args))
  log_step("Mesh with ", mesh$n, " nodes")

  # --- Gene sets -----------------------------------------------------------------
  status <- data.frame(gene_set = names(gene_sets),
                       genes_in_set = lengths(gene_sets),
                       genes_used = NA_integer_, status = NA_character_,
                       message = NA_character_, stringsAsFactors = FALSE)
  results <- list()

  for (s in seq_along(gene_sets)) {
    set_name <- names(gene_sets)[s]
    genes <- intersect(gene_sets[[s]], available)
    status$genes_used[s] <- length(genes)
    progress <- paste0("[", s, "/", length(gene_sets), "] ", set_name)
    if (length(genes) < max(2L, min_genes)) {
      status$status[s] <- "skipped"
      status$message[s] <- paste("fewer than", max(2L, min_genes), "usable genes")
      log_step(progress, ": skipped (", length(genes), " genes)")
      next
    }

    log_step(progress, ": ", length(genes), " genes")
    set_cfg <- modifyList(gene_set_cfg(length(genes), baseline), cfg)
    summary <- tryCatch({
      fit <- if (mode == "observed") {
        screen_joint_fit(expression[, genes, drop = FALSE], xy, region, mesh, h,
                         mode = "observed", cfg = set_cfg)
      } else {
        screen_joint_fit(as.matrix(counts[, genes, drop = FALSE]), xy, region,
                         mesh, h, mode = "adjusted",
                         library_size = adjusted$library_size,
                         dispersion_scale = adjusted$dispersion_scale[genes],
                         composition_offset = adjusted$log_m[, genes, drop = FALSE],
                         cfg = set_cfg)
      }
      screen_summarize(fit)
    }, error = function(e) e)

    if (inherits(summary, "error")) {
      status$status[s] <- "failed"
      status$message[s] <- conditionMessage(summary)
      log_step("  failed: ", conditionMessage(summary))
      next
    }
    status$status[s] <- "done"
    results[[set_name]] <- summary
    if (!is.null(out)) screen_write(summary, file.path(out, make.names(set_name)))
  }

  # --- Combined tables -------------------------------------------------------------
  combine <- function(table) {
    rows <- lapply(names(results), function(set_name) {
      x <- results[[set_name]]$pi[[table]]
      if (is.null(x)) return(NULL)
      cbind(gene_set = set_name,
            genes = status$genes_used[status$gene_set == set_name], x)
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) do.call(rbind, rows) else NULL
  }
  res <- list(overall = combine("overall"), pairwise = combine("pairwise"),
              gene_sets = status, results = results)

  if (!is.null(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(status, file.path(out, "gene_sets.csv"), row.names = FALSE)
    if (!is.null(res$overall)) {
      utils::write.csv(res$overall, file.path(out, "overall_scores.csv"), row.names = FALSE)
      utils::write.csv(res$pairwise, file.path(out, "pairwise_scores.csv"), row.names = FALSE)
    }
  }
  log_step("Done: ", sum(status$status == "done"), " fitted, ",
           sum(status$status == "skipped"), " skipped, ",
           sum(status$status == "failed"), " failed")
  invisible(res)
}


#' Reference expression of each cell type from single-cell counts
#'
#' Pools the raw counts of all cells of each cell type and scales them to
#' counts per 10,000, giving the `reference` of [run_spapath()] and
#' [composition_offset()].
#'
#' @param counts Cells x genes matrix (dense or sparse) of raw counts with
#'   gene colnames.
#' @param cell_type Cell type of each cell; cells with `NA` are dropped.
#' @return Cell types x genes matrix.
#' @export
reference_profile <- function(counts, cell_type) {
  stopifnot(length(cell_type) == nrow(counts), !is.null(colnames(counts)))
  keep <- !is.na(cell_type)
  counts <- counts[keep, , drop = FALSE]
  cell_type <- as.character(cell_type[keep])

  pooled <- sapply(split(seq_along(cell_type), cell_type),
                   function(rows) Matrix::colSums(counts[rows, , drop = FALSE]))
  t(sweep(pooled, 2, colSums(pooled), "/") * 1e4)
}


#' Default per-gene-set options of screen_joint_fit()
#'
#' Sets of 150 genes or more use one thread and the low-memory reordering.
#' @noRd
gene_set_cfg <- function(n_genes, baseline) {
  large <- n_genes >= 150
  cfg <- list(threads = if (large) "1:1" else "2:1", low_memory = large)
  if (!is.null(baseline)) cfg$baseline_region <- baseline
  cfg
}


#' Coordinates shifted to start at zero and divided by their largest value
#' @noRd
unit_square <- function(xy) {
  xy <- sweep(xy, 2, apply(xy, 2, min))
  xy / max(xy)
}


#' SCTransform Pearson residuals for selected genes
#'
#' Runs sctransform v2 on all genes (as Seurat's SCTransform does) and keeps
#' the residuals of `genes`, clipped at +-sqrt(N / 30) like Seurat.
#' @return Dense spots x genes matrix.
#' @noRd
sctransform_residuals <- function(counts, genes) {
  n_spots <- nrow(counts)
  umi <- Matrix::t(methods::as(counts, "CsparseMatrix"))
  vst_out <- sctransform::vst(umi, vst.flavor = "v2", verbosity = 0,
                              res_clip_range = c(-1, 1) * sqrt(n_spots / 30))
  residuals <- vst_out$y
  genes <- intersect(genes, rownames(residuals))
  t(residuals[genes, , drop = FALSE])
}


#' Inputs of the composition-adjusted mode
#'
#' @return List with `library_size`, `dispersion_scale` (named by gene),
#'   `log_m` (spots x shared genes) and the `usable` genes: shared by the
#'   counts and the reference, with a dispersion scale and positive predicted
#'   expression in every spot.
#' @noRd
adjusted_inputs <- function(counts, proportions, reference) {
  proportions <- as.matrix(proportions)
  if (is.null(rownames(proportions))) {
    stopifnot(nrow(proportions) == nrow(counts))
    rownames(proportions) <- rownames(counts)
  }
  missing <- setdiff(rownames(counts), rownames(proportions))
  if (length(missing)) stop(length(missing), " spots have no cell-type proportions.")
  proportions <- proportions[rownames(counts), , drop = FALSE]

  library_size <- Matrix::rowSums(counts)
  dispersion <- estimate_dispersion_scale(counts, library_size)
  scale <- setNames(dispersion$dispersion_scale, dispersion$gene)

  shared <- intersect(colnames(counts), colnames(reference))
  log_m <- composition_offset(proportions, reference[, shared, drop = FALSE])
  usable <- shared[!is.na(scale[shared]) & colSums(!is.finite(log_m)) == 0]

  list(library_size = library_size, dispersion_scale = scale,
       log_m = log_m, usable = usable)
}
