# =============================================================================
# map_activity.R
#
# One-call workflow: map the spatial activity of every gene set of a panel.
#
#   res <- map_activity(counts, coords, gene_sets)
#   res <- map_activity(counts, coords, gene_sets, gene_field = FALSE)   # fast
#
# Steps
#   field_fit() and field_regions() for each
#   gene set -> combined summary table
# =============================================================================


#' Map the spatial activity of a panel of gene sets
#'
#' Fits the spatial field of each gene set with [field_fit()] and calls its regions of higher and lower activity
#' with [field_regions()]. No region labels are needed.
#'
#' @param counts Spots x genes matrix (dense or sparse) of normalised
#'   expression, with spot rownames and gene colnames.
#' @param coords Spots x 2 matrix of spot coordinates, in the spot order of
#'   `counts`.
#' @param gene_sets Named list of gene vectors, or the path of a GMT file.
#' @param gene_field Add a spatial field of each gene (see [field_fit()]);
#'   `FALSE` is about 50 times faster but less accurate.
#' @param d,alpha Minimum effect and joint posterior level of the regions
#'   (see [field_regions()]).
#' @param out Optional output folder. When given, each gene set is written to
#'   its own subfolder (`spots.csv`, `hyperparameters.csv`), together with
#'   `gene_sets.csv`.
#' @param min_genes Gene sets with fewer usable genes are skipped.
#' @param threads INLA `num.threads`.
#' @param verbose Print progress messages.
#' @return Invisibly, a list with `gene_sets` (one row per gene set: genes
#'   used, status, sd of z, spots called higher and lower) and `results` (per
#'   fitted gene set, the [field_regions()] table and the hyperparameters).
#' @export
map_activity <- function(counts, coords, gene_sets, gene_field = TRUE, d = .1,
                         alpha = .05, out = NULL, min_genes = 5,
                         threads = "8:1", verbose = TRUE) {
  log_step <- function(...) if (verbose) message(format(Sys.time(), "[%H:%M:%S] "), ...)
  if (is.character(gene_sets)) gene_sets <- read_gmt(gene_sets)
  stopifnot(nrow(coords) == nrow(counts), !is.null(rownames(counts)), !is.null(colnames(counts)))

  pathway_genes <- unique(unlist(gene_sets))
  expression <- as.matrix(counts[, intersect(pathway_genes, colnames(counts)), drop = FALSE])
  log_step(ncol(expression), " of ", length(pathway_genes), " pathway genes are usable",
           if (gene_field) "; gene fields on" else "")

  status <- data.frame(gene_set = names(gene_sets), genes_in_set = lengths(gene_sets),
                       genes_used = NA, status = NA, message = NA, sd_z = NA,
                       spots_higher = NA, spots_lower = NA, fit_seconds = NA,
                       stringsAsFactors = FALSE)
  results <- list()
  for (i in seq_len(nrow(status))) {
    set_name <- status$gene_set[i]
    genes <- intersect(gene_sets[[set_name]], colnames(expression))
    status$genes_used[i] <- length(genes)
    progress <- paste0("[", i, "/", nrow(status), "] ", set_name)
    if (length(genes) < max(2, min_genes)) {
      status$status[i] <- "skipped"
      status$message[i] <- paste("fewer than", max(2, min_genes), "usable genes")
      log_step(progress, ": skipped (", length(genes), " genes)")
      next
    }
    log_step(progress, ": ", length(genes), " genes")
    fit <- tryCatch(field_fit(expression[, genes, drop = FALSE], coords, gene_field, threads),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      status$status[i] <- "failed"
      status$message[i] <- conditionMessage(fit)
      log_step("  failed: ", conditionMessage(fit))
      next
    }
    spots <- field_regions(fit, d, alpha)
    status[i, c("status", "sd_z", "spots_higher", "spots_lower", "fit_seconds")] <-
      list("done", fit$hyperparameters["sd_z", "mean"], sum(spots$call == "higher"),
           sum(spots$call == "lower"), fit$metadata$fit_seconds)
    results[[set_name]] <- list(spots = spots, hyperparameters = fit$hyperparameters)
    if (!is.null(out)) {
      folder <- file.path(out, make.names(set_name))
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(spots, file.path(folder, "spots.csv"), row.names = FALSE)
      utils::write.csv(fit$hyperparameters, file.path(folder, "hyperparameters.csv"))
    }
  }

  if (!is.null(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(status, file.path(out, "gene_sets.csv"), row.names = FALSE)
  }
  log_step("Done: ", sum(status$status == "done"), " fitted, ",
           sum(status$status == "skipped"), " skipped, ",
           sum(status$status == "failed"), " failed")
  invisible(list(gene_sets = status, results = results))
}
