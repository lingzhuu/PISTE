# =============================================================================
# cli.R
#
# Command-line interface: run PISTE on every gene set of a GMT file from a
# shell, without writing R code.
#
#   piste compare --input sample.h5ad --region layer --gmt sets.gmt \
#                   --out results/
#   piste compare -joint_sections --input sections.h5ad --region layer --section sample \
#                         --gmt sets.gmt --out results/
#   piste map --input sample.h5ad --coords positions.csv --gmt sets.gmt \
#               --out results/
#
# The `piste` executable (inst/bin/piste) calls piste_cli(); the same
# function can be called from R with a character vector of arguments.
#
# Pipeline
#   compare: read inputs -> align spots -> compare() (compare.R)
#   compare -joint_sections: the same, with the sections fitted jointly (observed mode)
#   map:     read inputs -> align spots -> map_activity() (map_activity.R)
#
# Several sections can be given to compare with --section; they are
# normalised separately and fitted on their own. compare -joint_sections fits adjacent sections
# of the same tissue jointly.
# =============================================================================


# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

#' Command-line options of `piste compare`, `piste compare -joint_sections` and `piste map`
#'
#' One row per option: name, default (NA = none), whether it is a flag, the
#' commands it belongs to and its help text.
#' @noRd
cli_option_table <- function(command) {
  rows <- list(
    c("input", NA, FALSE, "compare,joint,map", "Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized."),
    c("layer", NA, FALSE, "compare,joint,map", ".h5ad only: read layers/<LAYER> instead of X."),
    c("genes-in-rows", "false", TRUE, "compare,joint,map", ".csv/.tsv only: rows are genes and columns are spots."),
    c("region", NA, FALSE, "compare,joint", "Region labels: an obs column of the .h5ad, or a table with spot ID and region label."),
    c("coords", NA, FALSE, "map", "Spot coordinates: table with spot ID, x and y, or a Space Ranger tissue_positions file."),
    c("section", NA, FALSE, "compare", "Optional section labels: an obs column of the .h5ad, or a table with spot ID and section label; sections are normalised separately and fitted on their own."),
    c("section", NA, FALSE, "joint", "Section labels: an obs column of the .h5ad, or a table with spot ID and section label. The sections are fitted jointly, sharing the regional effects."),
    c("gmt", NA, FALSE, "compare,joint,map", "Gene sets in GMT format."),
    c("out", NA, FALSE, "compare,joint,map", "Output directory."),
    c("mode", "observed", FALSE, "compare", "observed (Gaussian model of normalised expression) or adjusted (NB model of raw counts adjusted for cell-type composition)."),
    c("normalized", "false", TRUE, "compare", "Observed mode: the input is already normalised; skip SCTransform."),
    c("normalized", "false", TRUE, "joint", "The input is already normalised; skip SCTransform."),
    c("normalized", "false", TRUE, "map", "The input is already normalised; skip SCTransform."),
    c("proportions", NA, FALSE, "compare", "Adjusted mode: table of cell-type proportions, spots x cell types, e.g. RCTD weights."),
    c("reference", NA, FALSE, "compare", "Adjusted mode: reference expression, either a table of cell types x genes on a linear scale (e.g. counts per 10,000) or a single-cell .h5ad/.h5/matrix directory of raw counts with --celltype-key."),
    c("celltype-key", NA, FALSE, "compare", "Adjusted mode: obs column (in a .h5ad reference) or table with cell ID and cell type."),
    c("target", NA, FALSE, "compare,joint", "Compare this region with all other spots, which are pooled into 'rest'."),
    c("baseline", NA, FALSE, "compare,joint", "Baseline region of the model (default: the first label found)."),
    c("no-gene-field", "false", TRUE, "map", "Leave out the spatial field of each gene: about 50 times faster but less accurate, as patterns of single genes can appear as pathway regions."),
    c("d", "0.1", FALSE, "map", "Minimum effect of a region, in units of the noise sd."),
    c("min-genes", "5", FALSE, "compare,joint,map", "Skip gene sets with fewer usable genes than this (at least 2)."),
    c("threads", NA, FALSE, "compare", "Adjusted mode: INLA threads, e.g. 4:1 (default: 1:1)."),
    c("threads", "8:1", FALSE, "map", "INLA threads."),
    c("alpha", "0.05", FALSE, "compare,joint", "A region is called higher or lower when the posterior probability of that direction is at least 1 - ALPHA/2."),
    c("alpha", "0.05", FALSE, "map", "A region of higher or lower activity has joint posterior probability at least 1 - ALPHA."),
    c("p-value", "false", TRUE, "compare,joint", "Also report p_value = 2 min(P(pi > 0.5), P(pi < 0.5)), the doubled posterior probability of the less likely direction (not a frequentist p-value)."),
    c("no-pi", "false", TRUE, "compare,joint", "Skip pi scores and evidence; report regional means only."),
    c("help", "false", TRUE, "compare,joint,map", "Show this help.")
  )
  table <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(table) <- c("name", "default", "flag", "command", "help")
  table$flag <- table$flag == "TRUE"
  table <- table[vapply(strsplit(table$command, ","), function(x) command %in% x, logical(1)), ]
  rownames(table) <- NULL
  table
}


#' Usage text of one command
#' @noRd
cli_usage <- function(command) {
  options <- cli_option_table(command)
  lines <- vapply(seq_len(nrow(options)), function(i) {
    o <- options[i, ]
    label <- if (o$flag) paste0("--", o$name) else paste0("--", o$name, " ", toupper(gsub("-", "_", o$name)))
    default <- if (!o$flag && !is.na(o$default)) paste0(" [default: ", o$default, "]") else ""
    paste0("  ", formatC(label, width = -30), " ", o$help, default)
  }, character(1))
  name <- if (command == "joint") "compare -joint_sections" else command
  c(paste0("Usage: piste ", name, " [options]"),
    "",
    switch(command,
           compare = "Compare regions on every gene set of a GMT file.",
           joint = "Compare regions on every gene set of a GMT file, fitting adjacent sections jointly.",
           map = "Map the spatial activity of every gene set of a GMT file."),
    "",
    "Options:",
    lines,
    "",
    "Other commands:",
    c(compare = "  piste compare --help                   Compare labelled regions.",
      joint = "  piste compare -joint_sections --help   Compare labelled regions, fitting adjacent sections jointly.",
      map = "  piste map --help                       Map spatial activity without region labels.")[
        setdiff(c("compare", "joint", "map"), command)],
    "  piste --version                        Print the package version.")
}


#' Parse `--name value`, `--name=value` and `--flag` arguments
#'
#' @return Named list with every option of cli_option_table(command); unset
#'   options without a default are NULL, flags are logical.
#' @noRd
cli_parse_args <- function(args, command) {
  options <- cli_option_table(command)
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
    if (is.na(row)) stop("Unknown option --", name, ". See piste ", command, " --help.")

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
cli_check_options <- function(opts, command) {
  required <- c("input", "gmt", "out", switch(command, map = "coords", joint = c("region", "section")))
  missing <- required[vapply(opts[required], is.null, logical(1))]
  if (length(missing)) stop("Missing required option(s): ", paste0("--", missing, collapse = ", "))

  opts$`min-genes` <- as.numeric(opts$`min-genes`)
  if (!is.finite(opts$`min-genes`) || opts$`min-genes` <= 0) {
    stop("--min-genes must be a positive number.")
  }
  opts$`min-genes` <- max(2, opts$`min-genes`)
  opts$alpha <- as.numeric(opts$alpha)
  if (!is.finite(opts$alpha) || opts$alpha <= 0 || opts$alpha >= 1) {
    stop("--alpha must satisfy 0 < alpha < 1.")
  }
  if (command == "map") {
    opts$d <- as.numeric(opts$d)
    if (!is.finite(opts$d) || opts$d < 0) stop("--d must be a non-negative number.")
    return(opts)
  }

  if (command == "joint") return(opts)
  opts$mode <- match.arg(opts$mode, c("observed", "adjusted"))
  if (opts$mode == "adjusted") {
    if (is.null(opts$proportions) || is.null(opts$reference)) {
      stop("--mode adjusted needs --proportions and --reference.")
    }
  }
  if (is.null(opts$region)) stop("Give the region labels with --region.")
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

  # A label option names an obs column of the .h5ad input or a table file.
  spot_labels <- function(value, what) {
    labels <- if (!is.null(data$obs) && value %in% names(data$obs)) {
      setNames(as.character(data$obs[[value]]), rownames(data$obs))
    } else if (file.exists(value)) {
      read_spot_table(value)
    } else {
      stop("--", what, " ", value, " is neither an obs column of ", opts$input, " nor a file")
    }
    labels[!is.na(labels) & nzchar(labels) & labels != "NA"]
  }
  region <- spot_labels(opts$region, "region")
  section <- if (!is.null(opts$section)) spot_labels(opts$section, "section")

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

#' Read expression and coordinates and align them by spot
#'
#' @return List with `counts` (spots x genes) and `coords` (spots x 2), in the
#'   same spot order.
#' @noRd
cli_load_map_inputs <- function(opts) {
  cli_log("Reading expression from ", opts$input)
  data <- read_expression(opts$input, layer = opts$layer,
                          genes_in_rows = isTRUE(opts$`genes-in-rows`))
  coords <- read_coordinates(opts$coords)
  spots <- intersect(rownames(data$counts), rownames(coords))
  if (length(spots) < 10) {
    stop("Only ", length(spots), " spots have expression and coordinates; check that the spot IDs match.")
  }
  cli_log(length(spots), " spots x ", ncol(data$counts), " genes")
  list(counts = data$counts[spots, , drop = FALSE], coords = coords[spots, , drop = FALSE])
}


#' Compare regions on every gene set with compare() and write the results
#'
#' `joint = TRUE` (`piste compare -joint_sections`) fits the sections jointly in observed mode.
#' @return The `gene_sets` status table of compare().
#' @noRd
cli_compare <- function(opts, joint = FALSE) {
  opts <- cli_check_options(opts, if (joint) "joint" else "compare")
  if (joint) opts$mode <- "observed"
  inputs <- cli_load_inputs(opts)
  gene_sets <- read_gmt(opts$gmt)
  cli_log(length(gene_sets), " gene sets in ", opts$gmt)

  tables <- if (opts$mode == "adjusted") cli_adjusted_tables(opts) else list()
  cfg <- list()
  if (isTRUE(opts$`no-pi`)) cfg$pi <- FALSE
  cfg$pi_alpha <- opts$alpha
  if (isTRUE(opts$`p-value`)) cfg$pi_report_p_value <- TRUE
  if (!is.null(opts$threads)) cfg$threads <- opts$threads

  res <- compare(inputs$counts, inputs$region, gene_sets,
                 mode = opts$mode, out = opts$out,
                 normalized = isTRUE(opts$normalized),
                 proportions = tables$proportions,
                 reference = tables$reference,
                 baseline = opts$baseline, section = inputs$section,
                 joint_sections = joint,
                 min_genes = opts$`min-genes`, cfg = cfg)
  cli_log("Results in ", opts$out)
  invisible(res$gene_sets)
}


#' Map every gene set with map_activity() and write the results
#' @return The `gene_sets` status table of map_activity().
#' @noRd
cli_map <- function(opts) {
  opts <- cli_check_options(opts, "map")
  inputs <- cli_load_map_inputs(opts)
  gene_sets <- read_gmt(opts$gmt)
  cli_log(length(gene_sets), " gene sets in ", opts$gmt)

  res <- map_activity(inputs$counts, inputs$coords, gene_sets,
                      gene_field = !isTRUE(opts$`no-gene-field`), d = opts$d,
                      alpha = opts$alpha, out = opts$out,
                      normalized = isTRUE(opts$normalized),
                      min_genes = opts$`min-genes`, threads = opts$threads)
  cli_log("Results in ", opts$out)
  invisible(res$gene_sets)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

#' Command-line entry point of PISTE
#'
#' Runs PISTE on every gene set of a GMT file. Called by the `piste`
#' executable installed with the package; see `piste compare --help`,
#' `piste compare -joint_sections --help` and `piste map --help`.
#'
#' Outputs in `--out`: one folder per gene set and `gene_sets.csv` with the
#' status of each set. `compare` and `compare -joint_sections` write the files of [screen_write()] and
#' `overall_scores.csv` and `pairwise_scores.csv` combining all gene sets;
#' `map` writes `spots.csv` and `hyperparameters.csv` (see [map_activity()]).
#'
#' @param args Character vector of command-line arguments, e.g.
#'   `c("compare", "--input", "sample.h5ad", ...)`.
#' @return Invisibly, the exit status (0 on success).
#' @export
piste_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[1] %in% c("-h", "--help", "help")) {
    writeLines(c(cli_usage("compare"), "", cli_usage("joint"), "", cli_usage("map")))
    return(invisible(0))
  }
  if (args[1] %in% c("--version", "version")) {
    writeLines(paste("PISTE", as.character(utils::packageVersion("PISTE"))))
    return(invisible(0))
  }
  command <- args[1]
  if (!command %in% c("compare", "map")) {
    stop("Unknown command '", command, "'. Use: piste compare [-joint_sections] or piste map [options]")
  }
  args <- args[-1]
  # `piste compare -joint_sections` fits adjacent sections jointly.
  if (command == "compare" && length(args) && args[1] == "-joint_sections") {
    command <- "joint"
    args <- args[-1]
  }

  opts <- cli_parse_args(args, command)
  if (isTRUE(opts$help)) {
    writeLines(cli_usage(command))
    return(invisible(0))
  }
  status <- switch(command, compare = cli_compare(opts), joint = cli_compare(opts, joint = TRUE),
                   map = cli_map(opts))
  invisible(if (any(status$status == "done")) 0 else 1)
}
