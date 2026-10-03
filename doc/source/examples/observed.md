# Observed mode on DLPFC

This example tests pathways in white matter (WM) against the cortical layers
of the human dorsolateral prefrontal cortex (DLPFC) Visium section 151673
from spatialLIBD (Maynard et al., 2021). The observed mode fits a Gaussian
model to normalized expression and reports the total regional difference in
pathway activity, including differences that come from cell composition.

## Requirements

Besides SpaPath ([Installation](../installation.md)), this example uses:

| Package | Source | Used for |
|---------|--------|----------|
| spatialLIBD | Bioconductor | Provides the DLPFC Visium data and its layer annotation, and plots spots on the H&E image |
| Seurat | CRAN | SCTransform normalization |

```r
install.packages(c("BiocManager", "Seurat"))
BiocManager::install("spatialLIBD")
```

The gene sets are read from a file in the SpaPath repository, so no gene-set
package is needed.

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

The section and its manual layer annotation, Layer1 to Layer6 and white
matter, can be shown next to the H&E image:

```r
library(ggplot2)

he <- imgRaster(spe, sample_id = "151673", image_id = "lowres")
png("151673_HE.png", width = 600, height = 600)
plot(he)
dev.off()

layers <- vis_clus(spe, sampleid = "151673",
                   clustervar = "layer_guess_reordered",
                   colors = libd_layer_colors)
ggsave("151673_layers.png", layers, width = 6, height = 6)
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

Coordinates are shifted to start at zero and divided by their largest value,
so the section fits in the unit square. The mesh and `h` are in these units.

```r
spots  <- match(rownames(expr), colnames(spe))
xy     <- spatialCoords(spe)[spots, ]
xy     <- sweep(xy, 2, apply(xy, 2, min))
xy     <- xy / max(xy)
region <- as.character(spe$layer_guess_reordered[spots])
table(region)
```

All inputs must follow the same spot order. To make WM the baseline region,
pass `cfg = list(baseline_region = "WM")` when fitting.

## 4. Gene sets

The DLPFC panel of 20 MSigDB gene sets is stored as a GMT file in the SpaPath
repository:

```r
gmt <- paste0("https://raw.githubusercontent.com/lingzhuu/pathwayLGM/main/",
              "benchmark/analysis/DLPFC/data/DLPFC_brain_and_reference.gmt")
fields    <- strsplit(readLines(gmt), "\t")
gene_sets <- setNames(lapply(fields, function(x) x[-(1:2)]),
                      vapply(fields, `[`, "", 1))
names(gene_sets)
```

## 5. Build the mesh

Build the mesh once for the section and reuse it for every gene set.

```r
h    <- 0.3
mesh <- fmesher::fm_mesh_2d(loc = xy, max.edge = 0.25, cutoff = 0.08,
                            offset = c(0.1, 0.2))
```

With the default settings `h` is the practical range of the spatial fields,
the distance at which their correlation falls to about 0.1.

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
