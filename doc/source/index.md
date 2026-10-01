# SpaPath

SpaPath tests whether a gene pathway is more active in one tissue region than in
another in spatial transcriptomics data. It models all genes of the pathway
jointly, gives every gene a smooth spatial random field to account for spatial
dependence between neighbouring spots, and fits the resulting latent Gaussian
model with INLA. The posterior is turned into a directional score and a
posterior evidence for every pair of regions and for every region against the
rest of the tissue.

SpaPath is the release name of the R package **PathwayLGM**. The package is
loaded with `library(PathwayLGM)`.

## What you need

| Input | Description |
|-------|-------------|
| Expression | Normalized expression of the pathway's genes, spots x genes. SCTransform is recommended. |
| Coordinates | Two-column matrix of spot coordinates. |
| Region labels | One label per spot, at least two regions. |
| Mesh and length scale | A 2D mesh over the tissue and the spatial length scale `h` used by the SPDE fields. |
| Covariates (optional) | Spot-level covariates such as batch. |

## What you get

For each ordered pair of regions *r* and *s*, and for each region against the
rest:

- **score**: the posterior mean of π<sub>rs</sub>, the probability that a pathway
  gene's state in region *r* exceeds its state in region *s*;
- **evidence**: the posterior probability that π<sub>rs</sub> > 0.5, i.e. that
  the pathway is more active in *r* than in *s*. Evidence near 1 means higher
  activity in *r*, evidence near 0 higher activity in *s*.

The package also reports regional means, their contrasts, variance components
and per-region calls. See [Output](output.md).

## Workflow in three lines

```r
fit     <- screen_joint_fit(Y, xy, region, mesh, h)
summary <- screen_summarize(fit)
screen_write(summary, "results/")
```

Start with [Installation](installation.md) and then the
[Quick start](quickstart.md).

!!! note "Planned features"
    The composition-adjusted mode described in the manuscript, which uses
    deconvolved cell-type proportions and a single-cell reference, is not yet
    implemented in the package. The current release fits the observed-expression
    model only.
