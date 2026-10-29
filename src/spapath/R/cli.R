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
#
# Several sections can be given in one run with --sections or --section-key;
# they are normalised separately and fitted on their own, or jointly with
# --joint-sections when they are adjacent sections of the same tissue.
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
    c("regions", NA, FALSE, "Table with spot ID and region label."),
    c("region-key", NA, FALSE, ".h5ad only: obs column with the region labels."),
    c("sections", NA, FALSE, "Optional table with spot ID and section label; sections are normalised separately and fitted on their own."),
    c("section-key", NA, FALSE, ".h5ad only: obs column with the section labels."),
    c("joint-sections", "false", TRUE, "Observed mode: fit the sections jointly, sharing the regional effects (adjacent sections of the same tissue)."),
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
    c("threads", NA, FALSE, "Adjusted mode: INLA threads, e.g. 4:1 (default: 1:1)."),
    c("alpha", "0.05", FALSE, "A region is called higher or lower when the posterior probability of that direction is at least 1 - ALPHA/2."),
    c("p-value", "false", TRUE, "Also report p_value = 2 min(P(pi > 0.5), P(pi < 0.5)), the doubled posterior probability of the less likely direction (not a frequentist p-value)."),
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
  if (isTRUE(opts$`joint-sections`)) {
    if (is.null(opts$sections) && is.null(opts$`section-key`)) {
      stop("--joint-sections needs the section labels (--sections or --section-key).")
    }
    if (opts$mode == "adjusted") stop("--joint-sections applies to the observed mode.")
  }

  opts$`min-genes` <- as.numeric(opts$`min-genes`)
  if (!is.finite(opts$`min-genes`) || opts$`min-genes` <= 0) {
    stop("--min-genes must be a positive number.")
  }
  opts$`min-genes` <- max(2, opts$`min-genes`)
  opts$alpha <- as.numeric(opts$alpha)
  if (!is.finite(opts$alpha) || opts$alpha <= 0 || opts$alpha >= 1) {
    stop("--alpha must satisfy 0 < alpha < 1.")
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

#' Read expression, regions and optional sections and align them by spot
#'
#' @return List with `counts` (spots x genes), `region` (character) and
#'   `section` (character or NULL), all in the same spot order.
#' @noRd
cli_load_inputs <- function(opts) {
  cli_log("Reading expression from ", opts$input)
  data <- read_expression(opts$input, layer = opts$layer,
                          genes_in_rows = isTRUE(opts$`genes-in-rows`))
  counts <- data$counts
  cli_log(nrow(counts), " spots x ", ncol(counts), " genes")

  spot_labels <- function(path, key, what) {
    labels <- if (!is.null(path)) {
      read_spot_table(path)
    } else {
      if (is.null(data$obs) || !key %in% names(data$obs)) {
        stop("--", what, "-key ", key, " is not an obs column of ", opts$input)
      }
      setNames(as.character(data$obs[[key]]), rownames(data$obs))
    }
    labels[!is.na(labels) & nzchar(labels) & labels != "NA"]
  }
  region <- spot_labels(opts$regions, opts$`region-key`, "region")
  section <- if (!is.null(opts$sections) || !is.null(opts$`section-key`)) {
    spot_labels(opts$sections, opts$`section-key`, "section")
  }

  spots <- intersect(rownames(counts), names(region))
  if (!is.null(section)) spots <- intersect(spots, names(section))
  if (length(spots) < 10) {
    stop("Only ", length(spots), " spots have expression and a region", if (!is.null(section)) " and section",
         " label; check that the spot IDs match.")
  }
  dropped <- nrow(counts) - length(spots)
  if (dropped > 0) cli_log("Dropped ", dropped, " spots without region or section label")

  region <- region[spots]
  if (!is.null(opts$target)) {
    if (!opts$target %in% region) stop("--target ", opts$target, " is not a region label.")
    region <- ifelse(region == opts$target, opts$target, "rest")
  }
  cli_log("Regions: ", paste(names(table(region)), table(region), sep = "=", collapse = ", "))
  if (!is.null(section)) {
    section <- section[spots]
    cli_log("Sections: ", paste(names(table(section)), table(section), sep = "=", collapse = ", "))
  }

  list(counts = counts[spots, , drop = FALSE], region = unname(region),
       section = if (is.null(section)) NULL else unname(section))
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
  cfg$pi_alpha <- opts$alpha
  if (isTRUE(opts$`p-value`)) cfg$pi_report_p_value <- TRUE
  if (!is.null(opts$threads)) cfg$threads <- opts$threads

  res <- run_spapath(inputs$counts, inputs$region, gene_sets,
                     mode = opts$mode, out = opts$out,
                     normalized = isTRUE(opts$normalized),
                     proportions = tables$proportions,
                     reference = tables$reference,
                     baseline = opts$baseline, section = inputs$section,
                     joint_sections = isTRUE(opts$`joint-sections`),
                     min_genes = opts$`min-genes`, cfg = cfg)
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
