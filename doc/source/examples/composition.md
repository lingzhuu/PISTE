# Composition-adjusted mode on DLPFC

The composition-adjusted mode reports the regional difference in pathway
activity that cannot be explained by differences in cell-type composition. It
models raw counts and uses, for each spot, the expression predicted from its
deconvolved cell-type proportions and a single-nucleus RNA-seq reference.

This example uses the same Visium section 151673 as the
[observed mode](observed.md) and the LIBD DLPFC snRNA-seq reference
(Huuki-Myers et al., 2024). Cell-type proportions are estimated with RCTD,
chosen here because it is fast. White matter (WM) is tested against all
other annotated spots.

## Requirements

Besides SpaPath ([Installation](../installation.md)), which already includes
DESeq2 for the dispersion step, this example uses:

| Package | Source | Used for |
|---------|--------|----------|
| spatialLIBD | Bioconductor | Provides the DLPFC Visium data and the LIBD DLPFC snRNA-seq reference |
| HDF5Array | Bioconductor | Loads the snRNA-seq reference, which is stored on disk in HDF5 format |
| spacexr | GitHub | RCTD deconvolution |

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("spatialLIBD", "HDF5Array"))
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

## 1. Load the Visium data and the snRNA-seq reference

```r
library(spatialLIBD)
library(SpatialExperiment)

spe <- fetch_data(type = "spe")
spe <- spe[, spe$sample_id == "151673"]
spe <- spe[, !is.na(spe$layer_guess_reordered)]
rownames(spe) <- make.unique(rowData(spe)$gene_name)

sce_zip <- fetch_data(type = "spatialDLPFC_snRNAseq")
unzip(sce_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)
sce <- sce[, !is.na(sce$cellType_broad_hc)]
table(sce$cellType_broad_hc)
```

The reference has seven broad cell types: Astro, EndoMural, Micro, Oligo, OPC,
Excit and Inhib. Check `colnames(colData(sce))` if the cell-type column has a
different name in your version.

## 2. Deconvolve with RCTD

RCTD is provided by the spacexr package. Build the reference and the spatial object from raw counts, then run RCTD in
full mode, which allows any number of cell types per spot:

```r
library(spacexr)
library(Matrix)

ref_counts <- as(counts(sce), "dgCMatrix")
cell_types <- factor(sce$cellType_broad_hc)
names(cell_types) <- colnames(sce)
reference <- Reference(ref_counts, cell_types, nUMI = colSums(ref_counts))

vis_counts <- as(counts(spe), "dgCMatrix")
coords <- as.data.frame(spatialCoords(spe))
colnames(coords) <- c("x", "y")
rownames(coords) <- colnames(spe)
puck <- SpatialRNA(coords, vis_counts, nUMI = colSums(vis_counts))

rctd <- create.RCTD(puck, reference, max_cores = 8)
rctd <- run.RCTD(rctd, doublet_mode = "full")

omega <- as.matrix(normalize_weights(rctd@results$weights))   # spots x cell types
```

Each row of `omega` holds the cell-type proportions of one spot and sums to one.

Alternatively, cell2location can be used in place of RCTD; it is more accurate but considerably slower, and a GPU is recommended (REF).

## 3. Reference mean expression

For each cell type, pool the raw counts of all its nuclei and scale them to
counts per 10,000:

```r
library(Matrix)

ref_counts <- as(counts(sce), "dgCMatrix")
rownames(ref_counts) <- make.unique(rowData(sce)$gene_name)
cell_type  <- sce$cellType_broad_hc

pooled <- sapply(split(seq_along(cell_type), cell_type),
                 function(cells) rowSums(ref_counts[, cells, drop = FALSE]))
nu <- t(sweep(pooled, 2, colSums(pooled), "/") * 1e4)   # cell types x genes
```

## 4. Spot inputs

```r
library(PathwayLGM)

C      <- t(as.matrix(counts(spe)))            # spots x genes, raw counts
region <- ifelse(spe$layer_guess_reordered == "WM", "WM", "rest")
depth  <- rowSums(C)                           # total UMI over all genes

omega <- omega[rownames(C), rownames(nu)]
```

`omega` holds the cell-type proportions of each spot, `nu` the reference mean
expression of each cell type, `C` the raw counts and `depth` the library
size.

Estimate the dispersion of every gene once, on the whole transcriptome:

```r
dispersion <- estimate_dispersion_scale(C, depth)
S <- setNames(dispersion$dispersion_scale, dispersion$gene)
```

Then compute the composition offset log m<sub>ig</sub> for the genes shared
by the section and the reference:

```r
shared <- intersect(colnames(C), colnames(nu))
log_m  <- composition_offset(omega, nu[, shared])
```

## 5. Coordinates, mesh and gene sets

Coordinates are shifted to start at zero and divided by their largest value,
so the section fits in the unit square. The mesh and `h` are in these units.

```r
xy <- spatialCoords(spe)
xy <- sweep(xy, 2, apply(xy, 2, min))
xy <- xy / max(xy)

h    <- 0.3
mesh <- fmesher::fm_mesh_2d(loc = xy, max.edge = 0.25, cutoff = 0.08,
                            offset = c(0.1, 0.2))
```

Use the gene sets of the [observed mode](observed.md#4-gene-sets).
The layer annotation and H&E image of the section are shown in the
[observed mode](observed.md#1-load-the-data).

## 6. Run SpaPath

Keep the pathway genes that pass the dispersion filter and have positive
predicted expression in every spot, then fit with `mode = "adjusted"`:

```r
for (set_name in names(gene_sets)) {
  genes <- intersect(gene_sets[[set_name]], shared)
  genes <- genes[!is.na(S[genes]) & colSums(!is.finite(log_m[, genes])) == 0]

  fit <- screen_joint_fit(
    C[, genes], xy, region, mesh, h,
    mode               = "adjusted",
    library_size       = depth,
    dispersion_scale   = S[genes],
    composition_offset = log_m[, genes],
    cfg = list(baseline_region = "rest",
               threads    = if (length(genes) >= 150) "1:1" else "2:1",
               low_memory = length(genes) >= 150)
  )
  res <- screen_summarize(fit)
  screen_write(res, file.path("results_adjusted", set_name))
}
```

Large gene sets use one thread and `low_memory`, which reorders the same
model to reduce memory.

## 7. Read the results

```r
res$pi$pairwise[res$pi$pairwise$region == "WM", ]
res$gene_coefficients
```

The WM row is read as in the [observed mode](observed.md#7-read-the-results),
but `score` and `evidence` now describe pathway activity after adjusting for
cell-type composition. `gene_coefficients` lists each gene's intercept `a_g`
and composition slope `lambda_g`; `lambda_g` near 1 means the gene follows
the reference prediction. See [API](../api.md#screen_summarize) for all
outputs.
