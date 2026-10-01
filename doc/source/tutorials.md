# Tutorials

## Pathway activity across tissue regions

This tutorial tests one pathway across the annotated regions of a spatial
transcriptomics section, starting from a Seurat object `obj` that holds the
raw counts, the spot coordinates and a region label for each spot in
`obj$region`.

### Normalize the expression

SpaPath models expression as Gaussian, so the counts must be normalized first.
We recommend SCTransform and use its Pearson residuals. By default SCTransform
keeps residuals only for variable genes; set `return.only.var.genes = FALSE`
so that every pathway gene is available.

```r
library(Seurat)
library(PathwayLGM)

obj  <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE)
expr <- t(as.matrix(GetAssayData(obj, assay = "SCT", layer = "scale.data")))
```

`expr` is a spots x genes matrix.

### Select the pathway genes

```r
pathway_genes <- c("MBP", "PLP1", "MOBP", "MAG", "MOG")   # your gene set
Y <- expr[, intersect(pathway_genes, colnames(expr))]
```

`Y` must have at least two genes and no missing values. Its column names are
used as gene labels in the output.

### Coordinates and regions

```r
coords <- GetTissueCoordinates(obj)
xy     <- as.matrix(coords[rownames(Y), c("x", "y")])
region <- as.character(obj$region[rownames(Y)])
```

All three inputs must follow the same spot order. Column names of the
coordinate table differ between Seurat versions; any two numeric coordinate
columns work.

The first region label found in `region` is the baseline region of the model.
To choose another one, pass `cfg = list(baseline_region = "WM")` to
`screen_joint_fit()`.

### Build the mesh

Each gene has a spatial field represented on a triangulated mesh. Build the
mesh once per section and reuse it for every pathway.

```r
h <- 300   # spatial length scale, in the units of xy
mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(100, 300), cutoff = 20)
```

The numbers above are placeholders; choose them for your coordinates and spot
spacing. With the default settings `h` is the practical range of the spatial
fields, the distance at which spatial correlation falls to about 0.1. A finer
mesh gives a more detailed field but a slower fit, since the model has
`mesh$n * K` spatial variables for K genes.

### Fit the model

```r
fit <- screen_joint_fit(Y, xy, region, mesh, h)
```

The fit integrates over the model hyperparameters with INLA's CCD design.
Increase INLA threads with `cfg = list(threads = "8:1")`.

### Summarise and save

```r
res <- screen_summarize(fit)
screen_write(res, "results/myelin")
```

### Read the results

```r
res$pi$overall
res$pi$pairwise
```

`overall` compares each region with the spot-weighted average of all other
regions, and `pairwise` compares every ordered pair of regions. Both tables
have the columns

| Column | Meaning |
|--------|---------|
| `score` | Posterior mean of π, the probability that a pathway gene's state in the region exceeds its state in the reference. 0.5 means no difference. |
| `evidence` | Posterior probability that π > 0.5. Near 1, the pathway is more active in the region; near 0, less active. |
| `boundary_mass` | Posterior mass at exactly π = 0.5, normally 0. |

`res$parameters` lists the regional means (`absolute_mu_r1`, ...), the
between-gene variance of each region (`tau2_r1`, ...), the spatial variance
`sigma_s2` and the noise variance `sigma_e2`, with posterior means and 95 %
intervals. Region indices follow `res$metadata$region_labels`, baseline first.

`res$mean_contrasts` and `res$domains` compare the regional means only;
`domains` calls a region `"higher"` or `"lower"` than the rest when the
posterior probability of a positive difference is at least 0.95 or at most
0.05.

`screen_write()` saves these tables as `overall_scores.csv`,
`pairwise_scores.csv`, `parameters.csv`, `mean_contrasts.csv` and
`domains.csv`, together with the full result in `screen_summary.rds`,
integration diagnostics and run times.
