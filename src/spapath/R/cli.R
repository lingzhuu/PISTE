# =============================================================================
# cli.R
#
# Command-line interface: run SpaPath on every gene set of a GMT file from a
# shell, without writing R code.
#
#   spapath run --input sample.h5ad --region-key layer --gmt sets.gmt \
#               --out results/
#
# The `spapath` executable (inst/bin/spapath) calls spapath_cli(); the same
# function can be called from R with a character vector of arguments.
#
# Pipeline
#   read inputs -> align spots -> normalise (observed) or build the
#   composition offset (adjusted) -> scale coordinates and build the mesh ->
#   fit, summarise and write each gene set -> combined summary tables
# =============================================================================


# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

#' Command-line options of `spapath run`
#'
#' One row per option: name, default (NA = none), whether it is a flag and
#' its help text.
#' @noRd
cli_option_table <- function() {
  rows <- list(
    c("input", NA, FALSE, "Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized."),
    c("layer", NA, FALSE, ".h5ad only: read layers/<LAYER> instead of X."),
    c("genes-in-rows", "false", TRUE, ".csv/.tsv only: rows are genes and columns are spots."),
    c("coords", NA, FALSE, "Spot coordinates: table with spot ID and x, y, or Space Ranger tissue_positions.csv. Not needed for .h5ad with obsm['spatial']."),
    c("spatial-key", "spatial", FALSE, ".h5ad only: obsm key of the coordinates."),
    c("regions", NA, FALSE, "Table with spot ID and region label."),
    c("region-key", NA, FALSE, ".h5ad only: obs column with the region labels."),
    c("gmt", NA, FALSE, "Gene sets in GMT format."),
    c("out", NA, FALSE, "Output directory."),
    c("mode", "observed", FALSE, "observed (Gaussian model of normalised expression) or adjusted (NB model of raw counts adjusted for cell-type composition)."),
    c("normalized", "false", TRUE, "Observed mode: the input is already normalised; skip SCTransform."),
    c("proportions", NA, FALSE, "Adjusted mode: table of cell-type proportions, spots x cell types, e.g. RCTD weights."),
    c("reference", NA, FALSE, "Adjusted mode: reference expression, either a table of cell types x genes on a linear scale (e.g. counts per 10,000) or a single-cell .h5ad/.h5/matrix directory of raw counts with --celltype-key."),
    c("celltype-key", NA, FALSE, "Adjusted mode: obs column (in a .h5ad reference) or table with cell ID and cell type."),
    c("target", NA, FALSE, "Compare this region with all other spots, which are pooled into 'rest'."),
    c("baseline", NA, FALSE, "Baseline region of the model (default: the first label found)."),
    c("min-genes", "5", FALSE, "Skip gene sets with fewer usable genes than this (at least 2)."),
    c("h", "0.3", FALSE, "Spatial length scale, in units of the coordinates scaled to the unit square."),
    c("max-edge", "0.25", FALSE, "Largest mesh triangle edge."),
    c("cutoff", "0.08", FALSE, "Minimum distance between mesh nodes."),
    c("mesh-offset", "0.1,0.2", FALSE, "Inner and outer mesh extension around the spots."),
    c("threads", NA, FALSE, "INLA threads, e.g. 4:1 (default: 2:1, or 1:1 with low memory for sets of 150 genes or more)."),
    c("no-pi", "false", TRUE, "Skip pi scores and evidence; report regional means only."),
    c("help", "false", TRUE, "Show this help.")
  )
  table <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(table) <- c("name", "default", "flag", "help")
  table$flag <- table$flag == "TRUE"
  table
}


#' Usage text of `spapath run`
#' @noRd
cli_usage <- function() {
  options <- cli_option_table()
  lines <- vapply(seq_len(nrow(options)), function(i) {
    o <- options[i, ]
    label <- if (o$flag) paste0("--", o$name) else paste0("--", o$name, " ", toupper(gsub("-", "_", o$name)))
    default <- if (!o$flag && !is.na(o$default)) paste0(" [default: ", o$default, "]") else ""
    paste0("  ", formatC(label, width = -30), " ", o$help, default)
  }, character(1))
  c("Usage: spapath run [options]",
    "",
    "Run SpaPath on every gene set of a GMT file.",
    "",
    "Options:",
    lines,
    "",
    "Other commands:",
    "  spapath --version             Print the package version.")
}


#' Parse `--name value`, `--name=value` and `--flag` arguments
#'
#' @return Named list with every option of cli_option_table(); unset options
#'   without a default are NULL, flags are logical.
#' @noRd
cli_parse_args <- function(args) {
  options <- cli_option_table()
  values <- setNames(vector("list", nrow(options)), options$name)
  for (i in seq_len(nrow(options))) {
    default <- options$default[i]
    values[i] <- list(if (options$flag[i]) FALSE else if (is.na(default)) NULL else default)
  }

  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (!startsWith(arg, "--")) stop("Unexpected argument: ", arg)
    name <- sub("^--", "", arg)
    value <- NULL
    if (grepl("=", name, fixed = TRUE)) {
      value <- sub("^[^=]*=", "", name)
      name <- sub("=.*$", "", name)
    }
    row <- match(name, options$name)
    if (is.na(row)) stop("Unknown option --", name, ". See spapath run --help.")

    if (options$flag[row]) {
      values[[name]] <- TRUE
    } else {
      if (is.null(value)) {
        i <- i + 1L
        if (i > length(args)) stop("Option --", name, " needs a value.")
        value <- args[i]
      }
      values[[name]] <- value
    }
    i <- i + 1L
  }
  values
}


#' Check option combinations and convert numeric options
#' @noRd
cli_check_options <- function(opts) {
  required <- c("input", "gmt", "out")
  missing <- required[vapply(opts[required], is.null, logical(1))]
  if (length(missing)) stop("Missing required option(s): ", paste0("--", missing, collapse = ", "))

  opts$mode <- match.arg(opts$mode, c("observed", "adjusted"))
  if (opts$mode == "adjusted") {
    if (is.null(opts$proportions) || is.null(opts$reference)) {
      stop("--mode adjusted needs --proportions and --reference.")
    }
    if (isTRUE(opts$normalized)) stop("--normalized applies to the observed mode only; the adjusted mode needs raw counts.")
  }
  if (is.null(opts$regions) && is.null(opts$`region-key`)) {
    stop("Give the region labels with --regions or --region-key.")
  }

  for (name in c("h", "max-edge", "cutoff", "min-genes")) {
    opts[[name]] <- as.numeric(opts[[name]])
    if (!is.finite(opts[[name]]) || opts[[name]] <= 0) stop("--", name, " must be a positive number.")
  }
  opts$`mesh-offset` <- as.numeric(strsplit(opts$`mesh-offset`, ",")[[1]])
  if (length(opts$`mesh-offset`) != 2L || any(!is.finite(opts$`mesh-offset`))) {
    stop("--mesh-offset must be two numbers separated by a comma, e.g. 0.1,0.2.")
  }
  opts$`min-genes` <- max(2, opts$`min-genes`)
  opts
}


#' Timestamped progress message on stderr
#' @noRd
cli_log <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), ...)
}


# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

#' Read expression, coordinates and regions and align them by spot
#'
#' @return List with `counts` (spots x genes), `xy` (spots x 2) and
#'   `region` (character), all in the same spot order.
#' @noRd
cli_load_inputs <- function(opts) {
  cli_log("Reading expression from ", opts$input)
  data <- read_expression(opts$input, layer = opts$layer,
                          genes_in_rows = isTRUE(opts$`genes-in-rows`))
  counts <- data$counts
  cli_log(nrow(counts), " spots x ", ncol(counts), " genes")

  xy <- if (!is.null(opts$coords)) {
    read_coordinates(opts$coords)
  } else if (!is.null(data$obsm[[opts$`spatial-key`]])) {
    coords <- data$obsm[[opts$`spatial-key`]][, 1:2, drop = FALSE]
    dimnames(coords) <- list(rownames(counts), c("x", "y"))
    coords
  } else {
    stop("No coordinates: give --coords, or an .h5ad with obsm['", opts$`spatial-key`, "'].")
  }

  region <- if (!is.null(opts$regions)) {
    read_spot_table(opts$regions)
  } else {
    key <- opts$`region-key`
    if (is.null(data$obs) || !key %in% names(data$obs)) {
      stop("--region-key ", key, " is not an obs column of ", opts$input)
    }
    setNames(as.character(data$obs[[key]]), rownames(data$obs))
  }
  region <- region[!is.na(region) & nzchar(region) & region != "NA"]

  spots <- Reduce(intersect, list(rownames(counts), rownames(xy), names(region)))
  if (length(spots) < 10L) {
    stop("Only ", length(spots), " spots have expression, coordinates and a region label; check that the spot IDs match.")
  }
  dropped <- nrow(counts) - length(spots)
  if (dropped > 0) cli_log("Dropped ", dropped, " spots without coordinates or region label")

  region <- region[spots]
  if (!is.null(opts$target)) {
    if (!opts$target %in% region) stop("--target ", opts$target, " is not a region label.")
    region <- ifelse(region == opts$target, opts$target, "rest")
  }
  if (length(unique(region)) < 2L) stop("At least two regions are needed.")
  cli_log("Regions: ", paste(names(table(region)), table(region), sep = "=", collapse = ", "))

  list(counts = counts[spots, , drop = FALSE],
       xy = xy[spots, , drop = FALSE],
       region = unname(region))
}


#' Coordinates shifted to start at zero and scaled to the unit square
#' @noRd
cli_scale_coordinates <- function(xy) {
  xy <- sweep(xy, 2, apply(xy, 2, min))
  xy / max(xy)
}


#' Mesh for the spatial fields
#' @noRd
cli_build_mesh <- function(xy, opts) {
  mesh <- fmesher::fm_mesh_2d(loc = xy, max.edge = opts$`max-edge`,
                              cutoff = opts$cutoff, offset = opts$`mesh-offset`)
  cli_log("Mesh with ", mesh$n, " nodes")
  mesh
}


#' SCTransform Pearson residuals for selected genes
#'
#' Runs sctransform v2 on all genes (as Seurat's SCTransform does) and keeps
#' the residuals of `genes`, clipped at +-sqrt(N / 30) like Seurat.
#' @return Dense spots x genes matrix.
#' @noRd
cli_sctransform <- function(counts, genes) {
  cli_log("Normalising with SCTransform (sctransform v2)")
  n_spots <- nrow(counts)
  umi <- Matrix::t(counts)
  vst_out <- sctransform::vst(umi, vst.flavor = "v2", verbosity = 0,
                              res_clip_range = c(-1, 1) * sqrt(n_spots / 30))
  residuals <- vst_out$y
  genes <- intersect(genes, rownames(residuals))
  t(residuals[genes, , drop = FALSE])
}


#' Reference expression of each cell type (cell types x genes)
#'
#' A table is used as is. Single-cell counts are pooled per cell type and
#' scaled to counts per 10,000.
#' @noRd
cli_reference_profile <- function(path, celltype_key) {
  name <- tolower(basename(path))
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", name)) {
    return(as.matrix(read_delimited(path, row_names = TRUE)))
  }

  if (is.null(celltype_key)) stop("A single-cell --reference needs --celltype-key.")
  reference <- read_expression(path)
  cell_type <- if (!is.null(reference$obs) && celltype_key %in% names(reference$obs)) {
    setNames(as.character(reference$obs[[celltype_key]]), rownames(reference$obs))
  } else if (file.exists(celltype_key)) {
    read_spot_table(celltype_key)
  } else {
    stop("--celltype-key must be an obs column of the reference .h5ad or a table of cell ID and cell type.")
  }
  cells <- intersect(rownames(reference$counts), names(cell_type)[!is.na(cell_type)])
  counts <- reference$counts[cells, , drop = FALSE]
  cell_type <- cell_type[cells]

  pooled <- sapply(split(seq_along(cell_type), cell_type),
                   function(rows) Matrix::colSums(counts[rows, , drop = FALSE]))
  t(sweep(pooled, 2, colSums(pooled), "/") * 1e4)
}


#' Inputs of the composition-adjusted mode
#'
#' @return List with `library_size`, `dispersion_scale` (named by gene),
#'   `log_m` (spots x shared genes) and the `usable` genes.
#' @noRd
cli_adjusted_inputs <- function(counts, opts) {
  proportions <- as.matrix(read_delimited(opts$proportions, row_names = TRUE))
  missing <- setdiff(rownames(counts), rownames(proportions))
  if (length(missing)) stop(length(missing), " spots have no cell-type proportions.")
  proportions <- proportions[rownames(counts), , drop = FALSE]

  nu <- cli_reference_profile(opts$reference, opts$`celltype-key`)
  if (!setequal(rownames(nu), colnames(proportions))) {
    stop("Cell types differ between --proportions (", paste(colnames(proportions), collapse = ", "),
         ") and --reference (", paste(rownames(nu), collapse = ", "), ").")
  }

  library_size <- Matrix::rowSums(counts)
  cli_log("Estimating gene dispersions with DESeq2 on ", ncol(counts), " genes")
  dispersion <- estimate_dispersion_scale(counts, library_size)
  scale <- setNames(dispersion$dispersion_scale, dispersion$gene)

  shared <- intersect(colnames(counts), colnames(nu))
  log_m <- composition_offset(proportions, nu[, shared, drop = FALSE])
  usable <- shared[!is.na(scale[shared]) & colSums(!is.finite(log_m)) == 0]
  cli_log(length(usable), " genes are usable in the adjusted mode")

  list(library_size = library_size, dispersion_scale = scale,
       log_m = log_m, usable = usable)
}


# -----------------------------------------------------------------------------
# Running gene sets
# -----------------------------------------------------------------------------

#' Fit, summarise and write one gene set
#' @return The summary from screen_summarize().
#' @noRd
cli_fit_gene_set <- function(genes, inputs, mesh, opts, out) {
  cfg <- list(pi = !isTRUE(opts$`no-pi`))
  if (!is.null(opts$baseline)) cfg$baseline_region <- opts$baseline
  large <- length(genes) >= 150
  cfg$threads <- if (!is.null(opts$threads)) opts$threads else if (large) "1:1" else "2:1"
  if (large) cfg$low_memory <- TRUE

  fit <- if (opts$mode == "observed") {
    screen_joint_fit(inputs$expression[, genes, drop = FALSE], inputs$xy,
                     inputs$region, mesh, opts$h, mode = "observed", cfg = cfg)
  } else {
    adjusted <- inputs$adjusted
    screen_joint_fit(as.matrix(inputs$counts[, genes, drop = FALSE]), inputs$xy,
                     inputs$region, mesh, opts$h, mode = "adjusted",
                     library_size = adjusted$library_size,
                     dispersion_scale = adjusted$dispersion_scale[genes],
                     composition_offset = adjusted$log_m[, genes, drop = FALSE],
                     cfg = cfg)
  }
  summary <- screen_summarize(fit)
  screen_write(summary, out)
  summary
}


#' Run every gene set and write the combined tables
#' @noRd
cli_run <- function(opts) {
  opts <- cli_check_options(opts)
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

  inputs <- cli_load_inputs(opts)
  gene_sets <- read_gmt(opts$gmt)
  cli_log(length(gene_sets), " gene sets in ", opts$gmt)
  pathway_genes <- unique(unlist(gene_sets))

  if (opts$mode == "observed") {
    inputs$expression <- if (isTRUE(opts$normalized)) {
      as.matrix(inputs$counts[, intersect(pathway_genes, colnames(inputs$counts)), drop = FALSE])
    } else {
      cli_sctransform(inputs$counts, pathway_genes)
    }
    available <- colnames(inputs$expression)
  } else {
    inputs$adjusted <- cli_adjusted_inputs(inputs$counts, opts)
    available <- inputs$adjusted$usable
  }

  inputs$xy <- cli_scale_coordinates(inputs$xy)
  mesh <- cli_build_mesh(inputs$xy, opts)

  status <- data.frame(gene_set = names(gene_sets),
                       genes_in_set = lengths(gene_sets),
                       genes_used = NA_integer_, status = NA_character_,
                       message = NA_character_, stringsAsFactors = FALSE)
  overall <- list()
  pairwise <- list()

  for (s in seq_along(gene_sets)) {
    set_name <- names(gene_sets)[s]
    genes <- intersect(gene_sets[[s]], available)
    status$genes_used[s] <- length(genes)
    if (length(genes) < opts$`min-genes`) {
      status$status[s] <- "skipped"
      status$message[s] <- paste("fewer than", opts$`min-genes`, "usable genes")
      cli_log("[", s, "/", length(gene_sets), "] ", set_name, ": skipped (", length(genes), " genes)")
      next
    }

    cli_log("[", s, "/", length(gene_sets), "] ", set_name, ": ", length(genes), " genes")
    folder <- file.path(opts$out, make.names(set_name))
    result <- tryCatch(cli_fit_gene_set(genes, inputs, mesh, opts, folder),
                       error = function(e) e)
    if (inherits(result, "error")) {
      status$status[s] <- "failed"
      status$message[s] <- conditionMessage(result)
      cli_log("  failed: ", conditionMessage(result))
      next
    }
    status$status[s] <- "done"
    if (!is.null(result$pi)) {
      overall[[set_name]] <- cbind(gene_set = set_name, genes = length(genes), result$pi$overall)
      pairwise[[set_name]] <- cbind(gene_set = set_name, genes = length(genes), result$pi$pairwise)
    }
  }

  out_file <- function(name) file.path(opts$out, name)
  utils::write.csv(status, out_file("gene_sets.csv"), row.names = FALSE)
  if (length(overall)) {
    utils::write.csv(do.call(rbind, overall), out_file("overall_scores.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, pairwise), out_file("pairwise_scores.csv"), row.names = FALSE)
  }
  cli_log("Done: ", sum(status$status == "done"), " fitted, ",
          sum(status$status == "skipped"), " skipped, ",
          sum(status$status == "failed"), " failed. Results in ", opts$out)
  invisible(status)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

#' Command-line entry point of SpaPath
#'
#' Runs SpaPath on every gene set of a GMT file. Called by the `spapath`
#' executable installed with the package; see `spapath run --help`.
#'
#' Outputs in `--out`: one folder per gene set with the files of
#' [screen_write()], `overall_scores.csv` and `pairwise_scores.csv` combining
#' all gene sets, and `gene_sets.csv` with the status of each set.
#'
#' @param args Character vector of command-line arguments, e.g.
#'   `c("run", "--input", "sample.h5ad", ...)`.
#' @return Invisibly, the exit status (0 on success).
#' @export
spapath_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[1] %in% c("-h", "--help", "help")) {
    writeLines(cli_usage())
    return(invisible(0L))
  }
  if (args[1] %in% c("--version", "version")) {
    writeLines(paste("SpaPath (PathwayLGM)", as.character(utils::packageVersion("PathwayLGM"))))
    return(invisible(0L))
  }
  if (args[1] != "run") stop("Unknown command '", args[1], "'. Use: spapath run [options]")

  opts <- cli_parse_args(args[-1])
  if (isTRUE(opts$help)) {
    writeLines(cli_usage())
    return(invisible(0L))
  }
  status <- cli_run(opts)
  invisible(if (any(status$status == "done")) 0L else 1L)
}
