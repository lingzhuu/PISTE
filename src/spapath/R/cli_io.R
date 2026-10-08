# =============================================================================
# cli_io.R
#
# Readers used by the command-line interface (see cli.R).
#
#   read_expression()   counts or normalised expression from .h5ad, 10x .h5,
#                       a 10x matrix directory, .csv/.tsv or .rds
#   read_coordinates()  spot coordinates from a table or Space Ranger
#                       tissue_positions file
#   read_spot_table()   two-column spot -> value table (e.g. region labels)
#   read_gmt()          gene sets
#
# Every expression reader returns spots x genes as a sparse dgCMatrix with
# spot rownames and gene colnames.
# =============================================================================


# -----------------------------------------------------------------------------
# Expression
# -----------------------------------------------------------------------------

#' Read expression from any supported format
#'
#' @param path File or directory.
#' @param layer For .h5ad: read `layers/<layer>` instead of `X`.
#' @param genes_in_rows For .csv/.tsv: TRUE when rows are genes and columns
#'   are spots (default: rows are spots).
#' @return List with `counts` (spots x genes dgCMatrix), `obs` (data frame of
#'   spot annotations or NULL) and `obsm` (named list of spot matrices or
#'   NULL).
#' @noRd
read_expression <- function(path, layer = NULL, genes_in_rows = FALSE) {
  if (dir.exists(path)) return(read_10x_mtx(path))
  if (!file.exists(path)) stop("Input not found: ", path)

  name <- tolower(basename(path))
  if (grepl("\\.h5ad$", name)) return(read_h5ad(path, layer))
  if (grepl("\\.h5$", name)) return(read_10x_h5(path))
  if (grepl("\\.rds$", name)) {
    counts <- readRDS(path)
    return(list(counts = as_spot_matrix(counts), obs = NULL, obsm = NULL))
  }
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", name)) {
    table <- read_delimited(path, row_names = TRUE)
    counts <- as.matrix(table)
    if (genes_in_rows) counts <- t(counts)
    return(list(counts = as_spot_matrix(counts), obs = NULL, obsm = NULL))
  }
  stop("Unsupported expression format: ", path,
       " (use .h5ad, 10x .h5, a 10x matrix directory, .csv/.tsv or .rds).")
}


#' Coerce a spots x genes matrix to a named dgCMatrix
#' @noRd
as_spot_matrix <- function(counts) {
  counts <- methods::as(methods::as(methods::as(counts, "dMatrix"),
                                    "generalMatrix"), "CsparseMatrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("The expression matrix needs spot rownames and gene colnames.")
  }
  if (anyDuplicated(colnames(counts))) {
    colnames(counts) <- make.unique(colnames(counts))
  }
  if (anyDuplicated(rownames(counts))) stop("Spot names must be unique.")
  counts
}


#' Delimited table with the first column as row names
#' @noRd
read_delimited <- function(path, row_names = FALSE, header = TRUE) {
  sep <- if (grepl("\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$", tolower(path))) "\t" else ","
  table <- utils::read.table(path, sep = sep, header = header,
                             check.names = FALSE, quote = "\"",
                             comment.char = "", stringsAsFactors = FALSE)
  if (row_names) {
    rownames(table) <- as.character(table[[1]])
    table <- table[-1]
  }
  table
}


#' 10x Genomics matrix directory (matrix.mtx, features/genes.tsv, barcodes.tsv)
#' @noRd
read_10x_mtx <- function(dir) {
  find <- function(stems) {
    candidates <- file.path(dir, c(stems, paste0(stems, ".gz")))
    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) stop("Missing ", stems[1], " in ", dir)
    hit[1]
  }
  open_file <- function(path) if (grepl("\\.gz$", path)) gzfile(path) else path

  matrix_file <- find("matrix.mtx")
  counts <- Matrix::readMM(open_file(matrix_file))       # genes x spots
  features <- utils::read.delim(open_file(find(c("features.tsv", "genes.tsv"))),
                                header = FALSE, stringsAsFactors = FALSE)
  barcodes <- readLines(open_file(find("barcodes.tsv")))

  gene_names <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  dimnames(counts) <- list(gene_names, barcodes)
  list(counts = as_spot_matrix(Matrix::t(counts)), obs = NULL, obsm = NULL)
}


#' Require the hdf5r package for HDF5-based formats
#' @noRd
require_hdf5r <- function() {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Reading .h5ad and 10x .h5 files requires the hdf5r package: ",
         "install.packages(\"hdf5r\").")
  }
}


#' Sparse matrix from CSR/CSC arrays stored in HDF5 (0-based indices)
#'
#' @param compressed "csr" (rows indexed by indptr) or "csc".
#' @noRd
hdf5_sparse <- function(data, indices, indptr, shape, compressed) {
  pointer_dim <- if (compressed == "csr") 1L else 2L
  index <- rep(seq_len(shape[pointer_dim]), diff(indptr))
  if (compressed == "csr") {
    Matrix::sparseMatrix(i = index, j = indices + 1L, x = as.numeric(data),
                         dims = shape)
  } else {
    Matrix::sparseMatrix(i = indices + 1L, j = index, x = as.numeric(data),
                         dims = shape)
  }
}


#' 10x Genomics feature-barcode matrix in HDF5 (Cell Ranger / Space Ranger)
#' @noRd
read_10x_h5 <- function(path) {
  require_hdf5r()
  file <- hdf5r::H5File$new(path, mode = "r")
  on.exit(file$close_all())
  if (!file$exists("matrix")) stop("Not a 10x feature-barcode h5 file: ", path)
  group <- file[["matrix"]]

  shape <- as.integer(group[["shape"]]$read())                # genes, spots
  counts <- hdf5_sparse(group[["data"]]$read(), group[["indices"]]$read(),
                        group[["indptr"]]$read(), shape, "csc")
  genes <- if (group$exists("features")) group[["features/name"]]$read() else
    group[["gene_names"]]$read()
  dimnames(counts) <- list(genes, group[["barcodes"]]$read())
  list(counts = as_spot_matrix(Matrix::t(counts)), obs = NULL, obsm = NULL)
}


#' AnnData .h5ad file: X or a layer, obs columns and obsm matrices
#' @noRd
read_h5ad <- function(path, layer = NULL) {
  require_hdf5r()
  file <- hdf5r::H5File$new(path, mode = "r")
  on.exit(file$close_all())

  obs <- h5ad_dataframe(file[["obs"]])
  var <- h5ad_dataframe(file[["var"]])
  key <- if (is.null(layer)) "X" else file.path("layers", layer)
  if (!file$exists(key)) stop("No ", key, " in ", path)
  counts <- h5ad_matrix(file[[key]], c(nrow(obs), nrow(var)))
  dimnames(counts) <- list(rownames(obs), rownames(var))

  obsm <- list()
  if (file$exists("obsm")) {
    for (name in names(file[["obsm"]])) {
      item <- file[["obsm"]][[name]]
      if (inherits(item, "H5D")) obsm[[name]] <- h5ad_dense(item)
    }
  }
  list(counts = as_spot_matrix(counts), obs = obs, obsm = obsm)
}


#' Dense HDF5 dataset in the row-major orientation written by h5py
#'
#' hdf5r returns dimensions in reverse order, so a 2D dataset is transposed.
#' @noRd
h5ad_dense <- function(dataset) {
  x <- dataset$read()
  if (is.matrix(x)) t(x) else x
}


#' AnnData matrix (dense or CSR/CSC sparse)
#' @noRd
h5ad_matrix <- function(item, shape) {
  if (inherits(item, "H5D")) return(Matrix::Matrix(h5ad_dense(item), sparse = TRUE))

  encoding <- if (item$attr_exists("encoding-type")) {
    hdf5r::h5attr(item, "encoding-type")
  } else if (item$attr_exists("h5sparse_format")) {
    paste0(hdf5r::h5attr(item, "h5sparse_format"), "_matrix")
  } else {
    stop("Unrecognised AnnData matrix encoding.")
  }
  compressed <- sub("_matrix$", "", encoding)
  if (!compressed %in% c("csr", "csc")) stop("Unsupported matrix encoding: ", encoding)
  hdf5_sparse(item[["data"]]$read(), item[["indices"]]$read(),
              item[["indptr"]]$read(), shape, compressed)
}


#' AnnData obs/var data frame, including categorical columns
#' @noRd
h5ad_dataframe <- function(group) {
  index_name <- if (group$attr_exists("_index")) hdf5r::h5attr(group, "_index") else "_index"
  index <- as.character(group[[index_name]]$read())
  columns <- setdiff(names(group), c(index_name, "__categories"))

  # anndata < 0.8 stores categories in a separate group.
  old_categories <- if (group$exists("__categories")) group[["__categories"]] else NULL

  frame <- data.frame(row.names = index)
  for (name in columns) {
    item <- group[[name]]
    if (inherits(item, "H5Group")) {
      if (!item$exists("codes")) next
      codes <- item[["codes"]]$read()
      categories <- as.character(item[["categories"]]$read())
      values <- categories[ifelse(codes < 0L, NA_integer_, codes + 1L)]   # -1 = missing
    } else {
      values <- item$read()
      if (!is.null(old_categories) && old_categories$exists(name)) {
        values <- as.character(old_categories[[name]]$read())[
          ifelse(values < 0L, NA_integer_, values + 1L)]
      }
    }
    if (length(values) == length(index)) frame[[name]] <- values
  }
  frame
}


# -----------------------------------------------------------------------------
# Spot tables and gene sets
# -----------------------------------------------------------------------------

#' Spot coordinates from a table or a Space Ranger tissue_positions file
#'
#' Accepted layouts: columns named `x` and `y`; Space Ranger columns
#' `pxl_col_in_fullres` and `pxl_row_in_fullres` (x = column, y = row); the
#' headerless six-column `tissue_positions_list.csv`; otherwise the first two
#' numeric columns after the spot ID.
#' @return N x 2 matrix with spot rownames and columns x, y.
#' @noRd
read_coordinates <- function(path) {
  table <- read_delimited(path)
  if (ncol(table) == 6L && !"pxl_row_in_fullres" %in% names(table) &&
      is.numeric(table[[5]]) && is.numeric(table[[6]]) &&
      !any(c("x", "y") %in% names(table))) {
    # Headerless tissue_positions_list.csv: the first row was read as header.
    table <- read_delimited(path, header = FALSE)
    names(table) <- c("barcode", "in_tissue", "array_row", "array_col",
                      "pxl_row_in_fullres", "pxl_col_in_fullres")
  }

  spots <- as.character(table[[1]])
  if (all(c("x", "y") %in% names(table))) {
    xy <- cbind(table$x, table$y)
  } else if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% names(table))) {
    xy <- cbind(table$pxl_col_in_fullres, table$pxl_row_in_fullres)
  } else {
    numeric_columns <- which(vapply(table[-1], is.numeric, logical(1))) + 1L
    if (length(numeric_columns) < 2L) stop("No coordinate columns found in ", path)
    xy <- as.matrix(table[numeric_columns[1:2]])
  }
  xy <- matrix(as.numeric(xy), ncol = 2, dimnames = list(spots, c("x", "y")))
  xy
}


#' Two-column table (spot ID, value) as a named character vector
#' @noRd
read_spot_table <- function(path) {
  table <- read_delimited(path)
  if (ncol(table) < 2L) stop(path, " needs a spot ID column and a value column.")
  setNames(as.character(table[[2]]), as.character(table[[1]]))
}


#' Gene sets from a GMT file (name, description, genes...)
#' @noRd
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t")
  fields <- fields[lengths(fields) >= 3L]
  sets <- lapply(fields, function(x) unique(x[-(1:2)][nzchar(x[-(1:2)])]))
  names(sets) <- vapply(fields, `[`, "", 1)
  sets
}
