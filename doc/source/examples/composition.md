# Composition-adjusted mode on DLPFC

The composition-adjusted mode reports the regional difference in pathway
activity that cannot be explained by differences in cell-type composition. It
models raw counts and uses, for each spot, the expression predicted from its
deconvolved cell-type proportions and a single-nucleus RNA-seq reference.

This example uses the same Visium section 151673 as the
[observed mode](observed.md) and the LIBD DLPFC snRNA-seq reference
(Huuki-Myers et al., 2024). Cell-type proportions are estimated with RCTD,
chosen here because it is fast.

!!! note
    The composition-adjusted model is not yet included in the released
    package. Steps 1 to 4 prepare its inputs; step 5 will be completed when
    the mode is released.

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

RCTD is provided by the spacexr package:

```r
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

Build the reference and the spatial object from raw counts, then run RCTD in
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

cell2location can be used instead of RCTD. It was more accurate in earlier
benchmarks (REF), but it is much slower and a GPU is recommended.

## 3. Reference mean expression

For each cell type, compute the depth-normalized mean expression of each gene
in the reference, as counts per 10,000:

```r
library(Matrix)

ref_counts <- as(counts(sce), "dgCMatrix")
ref_norm   <- t(t(ref_counts) / colSums(ref_counts)) * 1e4
cell_type  <- sce$cellType_broad_hc

nu <- sapply(split(seq_along(cell_type), cell_type),
             function(cells) rowMeans(ref_norm[, cells, drop = FALSE]))
nu <- t(nu)   # cell types x genes
```

## 4. Spot inputs

```r
C      <- t(as.matrix(counts(spe)))            # spots x genes, raw counts
xy     <- spatialCoords(spe)
region <- as.character(spe$layer_guess_reordered)

omega <- omega[rownames(C), rownames(nu)]
```

`omega` holds the cell-type proportions of each spot, `nu` the reference mean
expression of each cell type, and `C` the raw counts. Use the same gene sets
and mesh as in the [observed mode](observed.md#4-gene-sets).

## 5. Run SpaPath

To be added when the composition-adjusted mode is released.
