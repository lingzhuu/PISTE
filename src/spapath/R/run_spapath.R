# =============================================================================
# run_spapath.R
#
# One-call workflow: run SpaPath on every gene set of a panel.
#
#   res <- run_spapath(counts, region, gene_sets)                  # observed
#   res <- run_spapath(counts, region, gene_sets, section = sec,
#                      joint_sections = TRUE)                      # adjacent sections
#   res <- run_spapath(counts, region, gene_sets, mode = "adjusted",
#                      proportions = omega, reference = nu)        # adjusted
#
# Steps
#   normalise with SCTransform (observed) or build the dispersion scales and
#   composition offset (adjusted) -> fit, summarise (and write) each gene set
#   -> combined summary tables
#
# The step-by-step functions it calls (observed_fit(),
# estimate_dispersion_scale(), composition_offset(),
# screen_joint_fit(), screen_summarize(), screen_write()) remain available for
# custom workflows.
# =============================================================================


#' Run SpaPath on a panel of gene sets
#'
#' Runs the whole workflow in one call: normalisation (observed mode) or the
#' composition inputs (adjusted mode), and one fit per gene set. The observed
#' mode is fitted for each section on its own or, with
#' `joint_sections = TRUE`, jointly for the sections given by `section`,
#' which share the regional effects (see [observed_fit()]). Spot positions
#' are not used.
#'
#' @param counts Spots x genes matrix (dense or sparse) of raw counts, with
#'   spot rownames and gene colnames. In observed mode with
#'   `normalized = TRUE`, normalised expression instead.
#' @param region Region label of each spot, in the spot order of `counts`.
#' @param gene_sets Named list of gene vectors, or the path of a GMT file.
#' @param mode `"observed"` (Gaussian model of SCTransform-normalised
#'   expression) or `"adjusted"` (negative binomial model of raw counts
#'   adjusted for cell-type composition).
#' @param out Optional output folder. When given, each gene set is written to
#'   its own subfolder with [screen_write()] (one subfolder per section when
#'   the sections are fitted separately), together with
#'   `overall_scores.csv`, `pairwise_scores.csv` and `gene_sets.csv`.
#' @param normalized Observed mode: `counts` is already normalised, so
#'   SCTransform is skipped.
#' @param proportions Adjusted mode: spots x cell types matrix of cell-type
#'   proportions, e.g. normalised RCTD weights, with spot rownames.
#' @param reference Adjusted mode: cell types x genes matrix of reference
#'   mean expression on a linear scale, e.g. from [reference_profile()].
#' @param baseline Baseline region of the model (default: the first label).
#' @param section Optional section label of each spot, in the spot order of
#'   `counts`. Sections are normalised separately. By default every section
#'   is fitted on its own and the result tables get a `section` column.
#' @param joint_sections Observed mode: fit the sections given by `section`
#'   jointly, sharing the regional means and the gene-level regional effects
#'   and estimating a section-specific variance of the gene-level regional
#'   effects (`nu2_r` in the parameter table); π is defined on the shared
#'   component. The sections must come from the same tissue block and carry
#'   the same region labels.
#' @param min_genes Gene sets with fewer usable genes are skipped.
#' @param cfg Options: `pi_alpha` (0.05; a region is called higher or lower
#'   when the evidence in that direction is at least 1 - pi_alpha / 2),
#'   `pi_report_p_value` (FALSE; add the `p_value` column, see
#'   [pathway_pi_integrate()]), `pi` (FALSE skips the π scores) and
#'   `threads` (INLA threads, adjusted mode).
#' @param verbose Print progress messages.
#' @return Invisibly, a list with
#'   * `overall`, `pairwise`: the π tables of all gene sets, with columns
#'     `gene_set` and `genes` added (and `section` when sections are fitted
#'     separately);
#'   * `gene_sets`: one row per gene set (and section) with the genes used
#'     and whether it was fitted, skipped or failed;
#'   * `results`: the [screen_summarize()] output of each fitted gene set,
#'     named by gene set (and `<section>/<gene set>` when fitted separately).
#' @export
run_spapath <- function(counts, region, gene_sets,
                        mode = c("observed", "adjusted"), out = NULL,
                        normalized = FALSE, proportions = NULL,
                        reference = NULL, baseline = NULL, section = NULL,
                        joint_sections = FALSE, min_genes = 5, cfg = list(),
                        verbose = TRUE) {
  mode <- match.arg(mode)
  log_step <- function(...) if (verbose) message(format(Sys.time(), "[%H:%M:%S] "), ...)

  # --- Inputs ----------------------------------------------------------------
  if (is.character(gene_sets)) gene_sets <- read_gmt(gene_sets)
  stopifnot(is.list(gene_sets), length(gene_sets) > 0, !is.null(names(gene_sets)))
  region <- as.character(region)
  stopifnot(length(region) == nrow(counts))
  if (is.null(section)) {
    section <- rep("section", nrow(counts))
  } else {
    section <- as.character(section)
    stopifnot(length(section) == nrow(counts))
  }
  if (joint_sections && mode == "adjusted") stop("joint_sections applies to the observed mode.")

  # Spots without a region or section label are dropped.
  keep <- !is.na(region) & !is.na(section)
  if (!all(keep)) {
    log_step("Dropped ", sum(!keep), " spots without region or section label")
    counts <- counts[keep, , drop = FALSE]
    region <- region[keep]
    section <- section[keep]
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts needs spot rownames and gene colnames.")
  }
  if (mode == "adjusted") {
    if (is.null(proportions) || is.null(reference)) {
      stop("The adjusted mode needs proportions and reference.")
    }
    if (normalized) stop("normalized applies to the observed mode only.")
  }
  sections <- unique(section)
  per_section <- split(seq_along(section), factor(section, levels = sections))
  if (length(sections) > 1) {
    log_step(length(sections), " sections: ", paste(sections, collapse = ", "),
             if (joint_sections) " (fitted jointly)" else " (fitted separately)")
  }

  # --- Model inputs, per section ------------------------------------------------
  pathway_genes <- unique(unlist(gene_sets))
  inputs <- lapply(sections, function(s) {
    rows <- per_section[[s]]
    if (mode == "observed") {
      expression <- if (normalized) {
        as.matrix(counts[rows, intersect(pathway_genes, colnames(counts)), drop = FALSE])
      } else {
        log_step("Normalising ", if (length(sections) > 1) paste0("section ", s, " ") else "",
                 "with SCTransform")
        sctransform_residuals(counts[rows, , drop = FALSE], pathway_genes)
      }
      list(expression = expression, available = colnames(expression),
           region = region[rows])
    } else {
      log_step("Estimating gene dispersions and the composition offset",
               if (length(sections) > 1) paste0(" of section ", s) else "")
      section_proportions <- if (is.null(rownames(proportions))) {
        proportions[rows, , drop = FALSE]
      } else {
        proportions[rownames(counts)[rows], , drop = FALSE]
      }
      adjusted <- adjusted_inputs(counts[rows, , drop = FALSE], section_proportions,
                                  reference, region[rows])
      list(adjusted = adjusted, available = adjusted$usable, region = region[rows],
           counts = counts[rows, , drop = FALSE])
    }
  })
  names(inputs) <- sections
  available <- Reduce(intersect, lapply(inputs, `[[`, "available"))
  log_step(length(intersect(pathway_genes, available)), " of ",
           length(pathway_genes), " pathway genes are usable",
           if (length(sections) > 1) " in every section" else "")
  if (mode == "observed") {
    log_step("Observed mode: closed-form engine, no INLA")
  }

  # --- Gene sets -----------------------------------------------------------------
  fit_units <- if (joint_sections || length(sections) == 1) list(NULL) else as.list(sections)
  status <- do.call(rbind, lapply(fit_units, function(s) {
    data.frame(section = if (is.null(s)) NA_character_ else s,
               gene_set = names(gene_sets), genes_in_set = lengths(gene_sets),
               genes_used = NA, status = NA, message = NA, stringsAsFactors = FALSE)
  }))
  if (joint_sections || length(sections) == 1) status$section <- NULL
  results <- list()

  fit_one <- function(genes, unit, set_cfg) {
    if (joint_sections) {
      return(observed_fit(do.call(rbind, lapply(inputs, function(x) x$expression[, genes, drop = FALSE])),
                          unlist(lapply(inputs, `[[`, "region")),
                          rep(sections, vapply(inputs, function(x) length(x$region), 1L)),
                          cfg = set_cfg))
    }
    x <- inputs[[if (is.null(unit)) 1 else unit]]
    if (mode == "observed") {
      observed_fit(x$expression[, genes, drop = FALSE], x$region, cfg = set_cfg)
    } else {
      screen_joint_fit(as.matrix(x$counts[, genes, drop = FALSE]), x$region,
                       mode = "adjusted",
                       library_size = x$adjusted$library_size,
                       dispersion_scale = x$adjusted$dispersion_scale[genes],
                       composition_offset = x$adjusted$log_m[, genes, drop = FALSE],
                       cfg = set_cfg)
    }
  }

  for (i in seq_len(nrow(status))) {
    set_name <- status$gene_set[i]
    unit <- if (is.null(status$section)) NULL else status$section[i]
    genes <- intersect(gene_sets[[set_name]], available)
    status$genes_used[i] <- length(genes)
    label <- if (is.null(unit)) set_name else paste0(unit, "/", set_name)
    progress <- paste0("[", i, "/", nrow(status), "] ", label)
    if (length(genes) < max(2, min_genes)) {
      status$status[i] <- "skipped"
      status$message[i] <- paste("fewer than", max(2, min_genes), "usable genes")
      log_step(progress, ": skipped (", length(genes), " genes)")
      next
    }

    log_step(progress, ": ", length(genes), " genes")
    set_cfg <- modifyList(if (is.null(baseline)) list() else list(baseline_region = baseline), cfg)
    summary <- tryCatch(screen_summarize(fit_one(genes, unit, set_cfg)), error = function(e) e)

    if (inherits(summary, "error")) {
      status$status[i] <- "failed"
      status$message[i] <- conditionMessage(summary)
      log_step("  failed: ", conditionMessage(summary))
      next
    }
    status$status[i] <- "done"
    results[[label]] <- summary
    if (!is.null(out)) {
      screen_write(summary, if (is.null(unit)) file.path(out, make.names(set_name)) else
                     file.path(out, make.names(unit), make.names(set_name)))
    }
  }

  # --- Combined tables -------------------------------------------------------------
  combine <- function(table) {
    rows <- lapply(which(status$status == "done"), function(i) {
      unit <- if (is.null(status$section)) NULL else status$section[i]
      label <- if (is.null(unit)) status$gene_set[i] else paste0(unit, "/", status$gene_set[i])
      x <- results[[label]]$pi[[table]]
      if (is.null(x)) return(NULL)
      front <- data.frame(gene_set = status$gene_set[i], genes = status$genes_used[i],
                          stringsAsFactors = FALSE)
      if (!is.null(unit)) front <- cbind(section = unit, front, stringsAsFactors = FALSE)
      cbind(front, x)
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


#' SCTransform Pearson residuals for selected genes
#'
#' Runs sctransform v2 on all genes (as Seurat's SCTransform does) and keeps
#' the residuals of `genes`, clipped at +-sqrt(N / 30) like Seurat. Unlike
#' Seurat's SCTransform() (do.center = TRUE) the residuals are not centred
#' per gene; the model's shared gene intercept a_g absorbs each gene's
#' overall level, so centring is unnecessary.
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
adjusted_inputs <- function(counts, proportions, reference, region = NULL) {
  proportions <- as.matrix(proportions)
  if (is.null(rownames(proportions))) {
    stopifnot(nrow(proportions) == nrow(counts))
    rownames(proportions) <- rownames(counts)
  }
  missing <- setdiff(rownames(counts), rownames(proportions))
  if (length(missing)) stop(length(missing), " spots have no cell-type proportions.")
  proportions <- proportions[rownames(counts), , drop = FALSE]

  library_size <- Matrix::rowSums(counts)
  dispersion <- estimate_dispersion_scale(counts, library_size, region = region)
  scale <- setNames(dispersion$dispersion_scale, dispersion$gene)

  shared <- intersect(colnames(counts), colnames(reference))
  log_m <- composition_offset(proportions, reference[, shared, drop = FALSE])
  usable <- shared[!is.na(scale[shared]) & colSums(!is.finite(log_m)) == 0]

  list(library_size = library_size, dispersion_scale = scale,
       log_m = log_m, usable = usable)
}
