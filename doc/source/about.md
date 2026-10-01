# About

SpaPath is a Bayesian method for pathway analysis of spatial transcriptomics
data. The genes of a pathway can respond with different magnitudes and in
opposite directions, so SpaPath models all genes of a pathway jointly and
describes each region by the distribution of its gene-level expression states.

## Model

For spot *i* in region *r(i)* and gene *g* of the pathway, the normalized
expression is modelled as

```
y[i, g] = covariate effects
        + mu_r(i)          mean pathway level of the region
        + b_r(i)[g]        gene effect, independent across genes, variance tau2_r
        + w_g(s_i)         gene-specific Matérn spatial field
        + e[i, g]          Gaussian noise
```

The between-gene variance τ²<sub>r</sub> measures how consistently the pathway's
genes behave in region *r*. The spatial fields are represented with the SPDE
approach on a shared mesh, with the range fixed by the user, so the model has
R + 2 hyperparameters for R regions regardless of the number of genes. The
model is a latent Gaussian model fitted with INLA, integrating over the
hyperparameters with a CCD design.

A gene's state in region *r* is Z<sub>r</sub> ~ N(μ<sub>r</sub>, τ²<sub>r</sub>).
The directional score is the posterior mean of
π<sub>rs</sub> = P(Z<sub>r</sub> > Z<sub>s</sub>), and the evidence is the
posterior probability that π<sub>rs</sub> > 0.5.

## Planned features

A composition-adjusted mode, which accounts for differences in cell-type
composition between regions using deconvolved cell-type proportions and a
single-cell reference, is described in the manuscript but not yet available in
the package.

## Source code

[github.com/lingzhuu/pathwayLGM](https://github.com/lingzhuu/pathwayLGM),
package source in `src/spapath/`.
