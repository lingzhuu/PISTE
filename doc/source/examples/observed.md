# Observed mode on DLPFC

This example tests pathways in white matter (WM) against the cortical layers
of the human dorsolateral prefrontal cortex (DLPFC) Visium section 151673
from spatialLIBD (Maynard et al., 2021). The observed mode fits a Gaussian
model to normalized expression and reports the total regional difference in
pathway activity, including differences that come from cell composition.

Packages used in addition to SpaPath:

```r
install.packages(c("Seurat", "msigdbr"))
BiocManager::install("spatialLIBD")
```

## 1. Load the data

```r
library(spatialLIBD)
library(SpatialExperiment)
library(Seurat)
library(PathwayLGM)

spe <- fetch_data(type = "spe")
spe <- spe[, spe$sample_id == "151673"]

# Keep spots with a layer annotation (Layer1 to Layer6, WM).
spe <- spe[, !is.na(spe$layer_guess_reordered)]

counts <- counts(spe)
rownames(counts) <- make.unique(rowData(spe)$gene_name)
```

## 2. Normalize with SCTransform

SpaPath models normalized expression as Gaussian. We recommend SCTransform and
use its Pearson residuals. Set `return.only.var.genes = FALSE` so that
residuals are kept for every gene, not only the variable genes.

```r
obj <- CreateSeuratObject(counts = counts)
obj <- SCTransform(obj, return.only.var.genes = FALSE)

expr <- t(as.matrix(GetAssayData(obj, assay = "SCT", layer = "scale.data")))
```

`expr` is a spots x genes matrix of Pearson residuals.

## 3. Coordinates and regions

```r
spots  <- match(rownames(expr), colnames(spe))
xy     <- spatialCoords(spe)[spots, ]
region <- as.character(spe$layer_guess_reordered[spots])
table(region)
```

All inputs must follow the same spot order. To make WM the baseline region,
pass `cfg = list(baseline_region = "WM")` when fitting.

## 4. Gene sets

```r
library(msigdbr)

go_cc <- msigdbr(species = "Homo sapiens", collection = "C5",
                 subcollection = "GO:CC")
gene_sets <- split(go_cc$gene_symbol, go_cc$gs_name)[
  c("GOCC_MYELIN_SHEATH", "GOCC_SYNAPTIC_VESICLE")
]
```

Older versions of msigdbr use `category` and `subcategory` instead of
`collection` and `subcollection`.

## 5. Build the mesh

Build the mesh once for the section and reuse it for every gene set.

```r
h <- 1000   # placeholder: spatial length scale in pixels
mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(300, 1000),
                                 cutoff = 100)   # placeholders
```

`spatialCoords()` returns full-resolution pixel coordinates, so `h` and the
mesh sizes are in pixels. With the default settings `h` is the practical range
of the spatial fields, the distance at which their correlation falls to about
0.1.

## 6. Run SpaPath

```r
for (set_name in names(gene_sets)) {
  genes <- intersect(gene_sets[[set_name]], colnames(expr))
  fit <- screen_joint_fit(expr[, genes], xy, region, mesh, h,
                          cfg = list(baseline_region = "WM"))
  res <- screen_summarize(fit)
  screen_write(res, file.path("results_observed", set_name))
}
```

## 7. Read the results

```r
res$pi$overall
```

The `overall` table compares each region with the spot-weighted average of all
other regions. For the WM row, `score` is the posterior mean of π, how often a
pathway gene's expression state in WM exceeds its state in the cortical
layers, and `evidence` is the posterior probability that π > 0.5. Evidence
near 1 indicates higher activity in WM and evidence near 0 higher activity in
the cortical layers. `res$pi$pairwise` gives the same quantities for every
pair of regions. See [API](../api.md#screen_summarize) for all outputs.
