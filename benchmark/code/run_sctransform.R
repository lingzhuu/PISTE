# Run Seurat SCTransform (vst v2) on the dataset and save the results for
# downstream use:
#   pearson_residuals.rds / _f32.bin   dense residual matrix (+ float32 binary)
#   corrected_counts.rds               sparse UMI-corrected counts
#   genes.txt / spots.txt              row/column names of the matrices
#   preprocess_manifest.rds            provenance, settings and summary stats

# ---- config ---------------------------------------------------------------------
library(Matrix)
library(Seurat)

work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")

# ---- run ------------------------------------------------------------------------
out <- file.path(data_dir, "processed")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(file.path(data_dir, "raw_data.rds"))

# Retain low-count spots that may reflect sparse tissue regions and 
# use do.correct.umi to correct for sequencing depth differences.
seu <- CreateSeuratObject(counts = Matrix(x$counts, sparse = TRUE),
                          min.cells = 0, min.features = 0)
seu <- SCTransform(seu, assay = "RNA", new.assay.name = "SCT",
                   variable.features.n = 10000, return.only.var.genes = FALSE, 
                   do.correct.umi = TRUE, seed.use = 42, verbose = FALSE)

# ---- extract ------------------------------------------------------------------------
residual <- GetAssayData(seu, assay = "SCT", layer = "scale.data")
corrected <- GetAssayData(seu, assay = "SCT", layer = "counts")

# Align rows to the input gene order; residuals become a dense matrix.
residual <- as.matrix(residual[rownames(x$counts), , drop = FALSE])
corrected <- corrected[rownames(x$counts), , drop = FALSE]

# ---- save ----------------------------------------------------------------------------
saveRDS(residual, file.path(out, "pearson_residuals.rds"), compress = "xz")
saveRDS(corrected, file.path(out, "corrected_counts.rds"), compress = "xz")
# Column-major float32 dump for non-R readers; genes.txt/spots.txt give the order.
writeBin(as.numeric(residual), file.path(out, "pearson_residuals_f32.bin"),
         size = 4, endian = "little")
writeLines(rownames(residual), file.path(out, "genes.txt"))
writeLines(colnames(residual), file.path(out, "spots.txt"))
