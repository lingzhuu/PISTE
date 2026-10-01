# Preparing input

`screen_joint_fit()` takes one pathway at a time. All spot-level inputs must
use the same spot order.

## Expression `Y`

- Numeric matrix, **spots x genes** (N rows, K columns), with K ≥ 2 genes.
- Only the genes of the pathway being tested.
- All values finite: no `NA`, `NaN` or `Inf`.
- Column names are used as gene labels in the output, so set them to gene
  symbols.

The expression must be **normalized beforehand**, because the model treats it
as Gaussian. We recommend SCTransform. With Seurat, the SCTransform Pearson
residuals are stored in the `scale.data` layer of the `SCT` assay. By default
SCTransform keeps residuals only for variable genes, so make sure residuals are
available for all pathway genes, for example by running it with
`return.only.var.genes = FALSE`:

```r
library(Seurat)
obj  <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE)
expr <- t(as.matrix(GetAssayData(obj, assay = "SCT", layer = "scale.data")))
Y    <- expr[, intersect(pathway_genes, colnames(expr))]
```

## Coordinates `xy`

- Numeric matrix with N rows and 2 columns, all values finite.
- The units of `xy` are the units of `h` and of the mesh. Pixel coordinates
  and micrometres both work as long as they are used consistently.

## Region labels `region`

- Vector of length N with no `NA`, at least two distinct labels.
- Labels are converted to character. Factor levels are ignored; regions are
  ordered by their first appearance in `region`.
- The **baseline region** is the first label found, unless you set
  `cfg$baseline_region`. Its mean is the model intercept, and the other
  regions are parameterised as offsets from it.

## Mesh and length scale

The model gives every gene its own Matérn spatial field, represented with the
SPDE approach on a mesh shared by all genes.

- `mesh` is an INLA / fmesher 2D mesh covering the spots, for example from
  `fmesher::fm_mesh_2d_inla()`. The number of latent spatial variables is
  `mesh$n * K`, so a finer mesh costs more time and memory.
- `h` is a positive spatial length scale. The SPDE scale is fixed at
  `kappa = cfg$kappa_h / h`, and the practical range is
  `rho = sqrt(8) / kappa`. With the default `kappa_h = sqrt(8)`, the practical
  range equals `h`. The range is fixed rather than estimated, which keeps the
  number of hyperparameters at R + 2 regardless of the number of genes.

## Covariates `X` (optional)

- Numeric matrix with N rows, unique non-empty column names and finite values.
- An `intercept` column of ones is added automatically if absent; it is the
  mean of the baseline region. No other column may be constant one, so drop
  the `(Intercept)` column when you build `X` with `model.matrix()`:

```r
X <- model.matrix(~ batch, data = spot_meta)[, -1, drop = FALSE]
```

- Covariate effects are shared by all genes. Columns named in
  `cfg$gene_specific_covariates` instead get an independent slope for each gene
  (see [Options](options.md)).

## Offset (optional)

`offset` is a known additive offset of length N. The same value is added to
every gene of a spot. Leave it `NULL` for SCTransform-normalized input.
