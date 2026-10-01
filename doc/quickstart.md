# Quick start

This page runs SpaPath on one pathway. It assumes you already have normalized
expression, spot coordinates and region labels; see
[Preparing input](input.md) for how to build them.

## 1. Load the package and inputs

```r
library(PathwayLGM)

# Y:      spots x genes matrix of normalized expression for ONE pathway,
#         column names = gene symbols
# xy:     spots x 2 matrix of spot coordinates, same spot order as Y
# region: vector of region labels, same spot order as Y
stopifnot(nrow(Y) == nrow(xy), nrow(Y) == length(region))
```

## 2. Build the mesh and choose the length scale

The spatial fields are represented on a triangulated mesh. Build it once per
tissue section and reuse it for every pathway.

```r
# Placeholders: choose these for your coordinates and spot spacing.
h              <- 300   # spatial length scale, in the units of xy
max_edge_inner <- 100   # largest triangle edge inside the tissue
max_edge_outer <- 300   # largest triangle edge in the outer extension
cutoff         <- 20    # minimum distance between mesh nodes

mesh <- fmesher::fm_mesh_2d_inla(
  loc = xy,
  max.edge = c(max_edge_inner, max_edge_outer),
  cutoff = cutoff
)
mesh$n                # number of mesh nodes
```

These numbers are not recommendations; they only show the shape of the call.
With the
default options, `h` is the practical range of the spatial fields, the distance
at which spatial correlation drops to about 0.1. See
[Preparing input](input.md#mesh-and-length-scale).

## 3. Fit, summarise and write the results

```r
fit     <- screen_joint_fit(Y, xy, region, mesh, h)
summary <- screen_summarize(fit)
screen_write(summary, "results/my_pathway")
```

## 4. Read the directional results

```r
summary$pi$overall    # each region against the rest
summary$pi$pairwise   # every ordered pair of regions
```

Example of the `overall` table layout:

| region | spots | score | evidence | boundary_mass |
|--------|-------|-------|----------|---------------|
| WM     | ...   | ...   | ...      | ...           |
| L1     | ...   | ...   | ...      | ...           |

Evidence near 1 means the pathway is more active in that region than in the
rest; evidence near 0 means it is less active. See [Output](output.md) for all
tables and files.

## Running many pathways

Each call to `screen_joint_fit()` fits one pathway. To screen a gene-set
collection, loop over the sets, subset the expression matrix to the genes of
each set, and write each result to its own folder:

```r
for (set_name in names(gene_sets)) {
  genes <- intersect(gene_sets[[set_name]], colnames(expr))
  fit <- screen_joint_fit(expr[, genes], xy, region, mesh, h)
  screen_write(screen_summarize(fit), file.path("results", set_name))
}
```

Here `expr` is the spots x genes matrix of all genes and `gene_sets` a named
list of gene vectors. Each set needs at least two genes.
