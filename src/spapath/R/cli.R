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
#   read inputs -> align spots -> run_spapath() (see run_spapath.R)
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
    c("max-edge", NA, FALSE, "Largest mesh triangle edge, one value or inner,outer (default: h/3,h)."),
    c("cutoff", NA, FALSE, "Minimum distance between mesh nodes (default: h/9)."),
    c("mesh-offset", NA, FALSE, "Inner and outer mesh extension around the spots (default: h/2,1.5h)."),
    c("threads", NA, FALSE, "INLA threads, e.g. 4:1 (default: 1:1)."),
    c("delta", "0.05", FALSE, "Minimal effect: a region is called higher when P(pi > 0.5 + DELTA) > LEVEL and lower when P(pi < 0.5 - DELTA) > LEVEL."),
    c("level", "0.9", FALSE, "Posterior probability required for a call (see --delta)."),
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

  i <- 1
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
        i <- i + 1
        if (i > length(args)) stop("Option --", name, " needs a value.")
        value <- args[i]
      }
      values[[name]] <- value
    }
    i <- i + 1
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
  }
  if (is.null(opts$regions) && is.null(opts$`region-key`)) {
    stop("Give the region labels with --regions or --region-key.")
  }

  for (name in c("h", "min-genes")) {
    opts[[name]] <- as.numeric(opts[[name]])
  }
  if (!is.finite(opts$h) || opts$h <= 0) stop("--h must be a positive number.")
  # Mesh parameters: NULL means "derive from h" in build_mesh().
  parse_numbers <- function(value, name, n_allowed) {
    if (is.null(value)) return(NULL)
    x <- as.numeric(strsplit(as.character(value), ",")[[1]])
    if (any(!is.finite(x)) || !length(x) %in% n_allowed) {
      stop("--", name, " must be ", paste(n_allowed, collapse = " or "),
           " number(s) separated by a comma.")
    }
    x
  }
  opts$`max-edge` <- parse_numbers(opts$`max-edge`, "max-edge", 1:2)
  opts$cutoff <- parse_numbers(opts$cutoff, "cutoff", 1)
  if (!is.finite(opts$`min-genes`) || opts$`min-genes` <= 0) {
    stop("--min-genes must be a positive number.")
  }
  opts$`mesh-offset` <- parse_numbers(opts$`mesh-offset`, "mesh-offset", 2)
  opts$`min-genes` <- max(2, opts$`min-genes`)
  opts$delta <- as.numeric(opts$delta)
  opts$level <- as.numeric(opts$level)
  if (!is.finite(opts$delta) || opts$delta < 0 || opts$delta >= 0.5) {
    stop("--delta must satisfy 0 <= delta < 0.5.")
  }
  if (!is.finite(opts$level) || opts$level <= 0.5 || opts$level >= 1) {
    stop("--level must satisfy 0.5 < level < 1.")
  }
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
  if (length(spots) < 10) {
    stop("Only ", length(spots), " spots have expression, coordinates and a region label; check that the spot IDs match.")
  }
  dropped <- nrow(counts) - length(spots)
  if (dropped > 0) cli_log("Dropped ", dropped, " spots without coordinates or region label")

  region <- region[spots]
  if (!is.null(opts$target)) {
    if (!opts$target %in% region) stop("--target ", opts$target, " is not a region label.")
    region <- ifelse(region == opts$target, opts$target, "rest")
  }
  cli_log("Regions: ", paste(names(table(region)), table(region), sep = "=", collapse = ", "))

  list(counts = counts[spots, , drop = FALSE],
       xy = xy[spots, , drop = FALSE],
       region = unname(region))
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
  reference_profile(reference$counts[cells, , drop = FALSE], cell_type[cells])
}


#' Read the cell-type proportions and the reference of the adjusted mode
#' @noRd
cli_adjusted_tables <- function(opts) {
  proportions <- as.matrix(read_delimited(opts$proportions, row_names = TRUE))
  reference <- cli_reference_profile(opts$reference, opts$`celltype-key`)
  list(proportions = proportions, reference = reference)
}


# -----------------------------------------------------------------------------
# Running gene sets
# -----------------------------------------------------------------------------

#' Run every gene set with run_spapath() and write the results
#' @return The `gene_sets` status table of run_spapath().
#' @noRd
cli_run <- function(opts) {
  opts <- cli_check_options(opts)
  inputs <- cli_load_inputs(opts)
  gene_sets <- read_gmt(opts$gmt)
  cli_log(length(gene_sets), " gene sets in ", opts$gmt)

  tables <- if (opts$mode == "adjusted") cli_adjusted_tables(opts) else list()
  cfg <- list()
  if (isTRUE(opts$`no-pi`)) cfg$pi <- FALSE
  cfg$pi_delta <- opts$delta
  cfg$pi_level <- opts$level
  if (!is.null(opts$threads)) cfg$threads <- opts$threads

  res <- run_spapath(inputs$counts, inputs$xy, inputs$region, gene_sets,
                     mode = opts$mode, out = opts$out,
                     normalized = isTRUE(opts$normalized),
                     proportions = tables$proportions,
                     reference = tables$reference,
                     baseline = opts$baseline, h = opts$h,
                     min_genes = opts$`min-genes`,
                     mesh_args = Filter(Negate(is.null),
                                        list(max_edge = opts$`max-edge`,
                                             cutoff = opts$cutoff,
                                             offset = opts$`mesh-offset`)),
                     cfg = cfg)
  cli_log("Results in ", opts$out)
  invisible(res$gene_sets)
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
    return(invisible(0))
  }
  if (args[1] %in% c("--version", "version")) {
    writeLines(paste("SpaPath (PathwayLGM)", as.character(utils::packageVersion("PathwayLGM"))))
    return(invisible(0))
  }
  if (args[1] != "run") stop("Unknown command '", args[1], "'. Use: spapath run [options]")

  opts <- cli_parse_args(args[-1])
  if (isTRUE(opts$help)) {
    writeLines(cli_usage())
    return(invisible(0))
  }
  status <- cli_run(opts)
  invisible(if (any(status$status == "done")) 0 else 1)
}
