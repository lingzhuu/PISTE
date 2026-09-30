# PathwayLGM (Gaussian, 0.5.0)

R package that fits a Gaussian spatial pathway model with INLA and scores
directional pathway enrichment between tissue regions.

## Model

For N spots, K pathway genes and R regions (region 1 is the baseline):

```
y[i, g] = x_i' beta          shared covariates incl. baseline intercept
        + mu_r(i)            relative regional offset (mu_1 = 0)
        + gene_r(i)[g]       region-specific iid gene effect, variance tau2_r
        + x_i' slope_g       optional gene-specific covariate slopes
        + spatial_g(s_i)     per-gene SPDE Matern field (range fixed by h)
        + e[i, g]            Gaussian noise
```

Hyperparameters are integrated with INLA's CCD design. The directional score
for regions r and s is `pi_rs = P(Z_r > Z_s)` with `Z_r ~ N(mu_r, tau2_r)`;
the package reports its posterior mean (score) and `P(pi_rs > 0.5)`
(evidence), both pairwise and for each region against the rest.

## Usage

```r
library(PathwayLGM)

fit     <- screen_joint_fit(Y, xy, region, mesh, h)   # Y: spots x genes
summary <- screen_summarize(fit)
screen_write(summary, "results/")
```

## Source layout

| File                   | Contents                                                     |
|------------------------|--------------------------------------------------------------|
| `R/inla_screen.R`      | model building, INLA fit, summaries, file output             |
| `R/pi_engine.R`        | directional probabilities, CCD mixture extraction, integration |
| `R/positive_summary.R` | checked log-scale summaries of positive marginals            |

## Dependencies

INLA (from https://inla.r-inla-download.org/R/stable), Matrix, qrng;
fmesher is suggested for building meshes.
