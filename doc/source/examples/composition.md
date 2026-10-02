# Composition-adjusted mode on DLPFC

The composition-adjusted mode reports the regional difference in pathway
activity that cannot be explained by differences in cell-type composition. It
models raw counts and uses, for each spot, the expression predicted from its
deconvolved cell-type proportions and a single-nucleus RNA-seq reference.

This example uses the same Visium section 151673 as the
[observed mode](observed.md), deconvolved with cell2location using the LIBD
DLPFC snRNA-seq reference (Huuki-Myers et al., 2024).

!!! note
    The composition-adjusted model is not yet included in the released
    package. Steps 1 to 4 prepare its inputs; step 5 will be completed when
    the mode is released.

## 1. Load the Visium data and the snRNA-seq reference

```r
library(spatialLIBD)
library(SpatialExperiment)
library(zellkonverter)

spe <- fetch_data(type = "spe")
spe <- spe[, spe$sample_id == "151673"]
spe <- spe[, !is.na(spe$layer_guess_reordered)]
rownames(spe) <- make.unique(rowData(spe)$gene_name)

sce_zip <- fetch_data(type = "spatialDLPFC_snRNAseq")
unzip(sce_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)
table(sce$cellType_broad_hc)
```

The reference has seven broad cell types: Astro, EndoMural, Micro, Oligo, OPC,
Excit and Inhib. Check `colnames(colData(sce))` if the cell-type column has a
different name in your version.

Export both objects with raw counts for cell2location:

```r
spe_out <- SingleCellExperiment(assays = list(counts = counts(spe)),
                                colData = colData(spe))
writeH5AD(spe_out, "dlpfc_151673.h5ad", X_name = "counts")

sce <- sce[, !is.na(sce$cellType_broad_hc)]
writeH5AD(sce, "dlpfc_snrna.h5ad", X_name = "counts")
```

## 2. Deconvolve with cell2location

cell2location is a Python package (`pip install cell2location`). First
estimate the expression signature of each cell type from the reference, then
the cell-type abundance in each spot.

```python
import numpy as np
import scanpy as sc
import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

# Reference signatures
adata_ref = sc.read_h5ad("dlpfc_snrna.h5ad")
selected = filter_genes(adata_ref, cell_count_cutoff=5,
                        cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

RegressionModel.setup_anndata(adata_ref, batch_key="Sample",
                              labels_key="cellType_broad_hc")
ref_model = RegressionModel(adata_ref)
ref_model.train(max_epochs=250)
adata_ref = ref_model.export_posterior(
    adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
)

factors = adata_ref.uns["mod"]["factor_names"]
inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{f}" for f in factors]
].copy()
inf_aver.columns = factors

# Spatial mapping
adata_vis = sc.read_h5ad("dlpfc_151673.h5ad")
genes = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, genes].copy()
inf_aver = inf_aver.loc[genes, :]

cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
model = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=5,   # expected cells per spot; set for your tissue
    detection_alpha=20,
)
model.train(max_epochs=30000, batch_size=None, train_size=1)
adata_vis = model.export_posterior(
    adata_vis, sample_kwargs={"num_samples": 1000,
                              "batch_size": model.adata.n_obs}
)

# Cell-type proportions per spot
abundance = adata_vis.obsm["q05_cell_abundance_w_sf"]
abundance.columns = factors
proportions = abundance.div(abundance.sum(axis=1), axis=0)
proportions.to_csv("proportions_151673.csv")
```

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
omega  <- as.matrix(read.csv("proportions_151673.csv", row.names = 1))
C      <- t(as.matrix(counts(spe)))            # spots x genes, raw counts
xy     <- spatialCoords(spe)
region <- as.character(spe$layer_guess_reordered)

stopifnot(identical(rownames(omega), rownames(C)))
omega <- omega[, rownames(nu)]
```

`omega` holds the cell-type proportions of each spot, `nu` the reference mean
expression of each cell type, and `C` the raw counts. Use the same gene sets
and mesh as in the [observed mode](observed.md#4-gene-sets).

## 5. Run SpaPath

To be added when the composition-adjusted mode is released.
