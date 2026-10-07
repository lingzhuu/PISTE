# Run Seurat SCTransform (vst v2) on the dataset and save the results for
# downstream use:
#   pearson_residuals.rds / _f32.bin   dense residual matrix (+ float32 binary)
#   corrected_counts.rds               sparse UMI-corrected counts
#   sct_models.rds                     fitted SCT models
#   genes.txt / spots.txt              row/column names of the matrices
#   preprocess_manifest.rds            provenance, settings and summary stats

# ---- config ---------------------------------------------------------------------
library(Matrix)
library(Seurat)

work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")
sct_seed <- 42L

# ---- run ------------------------------------------------------------------------
out <- file.path(data_dir, "SCTransform_v2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(file.path(data_dir, "processed_data.rds"))

set.seed(sct_seed)
start <- proc.time()[["elapsed"]]
seu <- CreateSeuratObject(counts = Matrix(x$counts, sparse = TRUE),
                          min.cells = 0, min.features = 0)
seu <- SCTransform(seu, assay = "RNA", new.assay.name = "SCT",
                   vst.flavor = "v2", variable.features.n = 10000L,
                   return.only.var.genes = FALSE, do.correct.umi = TRUE,
                   seed.use = sct_seed, verbose = FALSE)

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
         size = 4L, endian = "little")
writeLines(rownames(residual), file.path(out, "genes.txt"))
writeLines(colnames(residual), file.path(out, "spots.txt"))
saveRDS(seu[["SCT"]]@SCTModel.list, file.path(out, "sct_models.rds"),
        compress = "xz")

meta <- list(input = file.path(data_dir, "processed_data.rds"),
             dataset_seed = x$seed, n_genes = nrow(residual),
             n_spots = ncol(residual),
             Seurat = as.character(packageVersion("Seurat")),
             sctransform = as.character(packageVersion("sctransform")),
             method = "Seurat::SCTransform", vst_flavor = "v2",
             variable_features_n = 10000L, return_only_var_genes = FALSE,
             do_correct_umi = TRUE, seed = sct_seed,
             residual_min = min(residual), residual_max = max(residual),
             residual_mean = mean(residual),
             corrected_nnz = length(corrected@x),
             seconds = proc.time()[["elapsed"]] - start)
saveRDS(meta, file.path(out, "preprocess_manifest.rds"))
