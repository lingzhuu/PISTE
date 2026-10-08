# Run gsdensity on SCTransform output, as in the authors' 10x Visium example,
# and save:
#   gsdensity_deviation.csv             compute.kld result for all pathways
#   gsdensity_PAL_binary_labels.csv     spots x deviated pathways, binarized
#
# compute.mca() applies NormalizeData() to the counts layer of the SCT assay
# (the corrected counts) and then runs MCA; Pearson residuals are not used.

# ---- config ---------------------------------------------------------------------
library(Matrix)
library(Seurat)
library(gsdensity)
set.seed(42)

work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")

# ---- load -----------------------------------------------------------------------
x <- readRDS(file.path(data_dir, "raw_data.rds"))

# Same SCTransform settings as run_sctransform.R.
seu <- CreateSeuratObject(counts = Matrix(x$counts, sparse = TRUE),
                          min.cells = 0, min.features = 0)
seu <- SCTransform(seu, assay = "RNA", new.assay.name = "SCT",
                   variable.features.n = 10000, return.only.var.genes = FALSE,
                   do.correct.umi = TRUE, seed.use = 42, verbose = FALSE)

# ---- run ------------------------------------------------------------------------
ce <- compute.mca(object = seu, assay = "SCT")
genes_for_grid <- intersect(rownames(seu), rownames(ce))
sets <- lapply(x$pathways, intersect, y = genes_for_grid)
res <- compute.kld(coembed = ce, genes.use = genes_for_grid,
                   n.grids = 100, gene.set.list = sets,
                   gene.set.cutoff = 20, n.times = 100)
write.csv(res, file.path(data_dir, "gsdensity_deviation.csv"), row.names = FALSE)
gene.set.deviated <- res[!is.na(res$p.adj) & res$p.adj < 0.05, ]$gene.set

el <- compute.nn.edges(coembed = ce, nn.use = 300)
cv.df <- as.matrix(run.rwr.list(el = el, gene_set_list = sets[gene.set.deviated],
                                cells = colnames(seu)))
cv.df <- cv.df[colnames(seu), gene.set.deviated, drop = FALSE]
# Dividing each pathway's PAL by its mean keeps ranks and modes unchanged and
# avoids the multimode default tolerance dropping labels on the raw scale.
cl.df <- matrix(compute.cell.label.df(sweep(cv.df, 2, colMeans(cv.df), "/")),
                nrow = nrow(cv.df), ncol = ncol(cv.df), dimnames = dimnames(cv.df))

# ---- save -----------------------------------------------------------------------
write.csv(data.frame(spot = rownames(cl.df), cl.df, check.names = FALSE),
          file.path(data_dir, "gsdensity_PAL_binary_labels.csv"), row.names = FALSE)
cat("Selected", length(gene.set.deviated), "of", length(sets), "pathways\n")
