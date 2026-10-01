# Options

Options are passed to `screen_joint_fit()` as a list in `cfg`. Only the entries
you set are changed; all others keep their defaults.

```r
fit <- screen_joint_fit(Y, xy, region, mesh, h,
                        cfg = list(threads = "8:1", baseline_region = "WM"))
```

## Directional probability

| Option | Default | Description |
|--------|---------|-------------|
| `pi` | `TRUE` | Compute directional scores and evidence. Set `FALSE` to fit the model and report only regional means, contrasts and variance components. |
| `pi_nodes` | `4096` | Initial number of quadrature nodes for the evidence with three or more regions (power of two, ≥ 16). |
| `pi_max_nodes` | `65536` | Maximum number of quadrature nodes (power of two, ≥ 2 × `pi_nodes`). |
| `pi_tolerance` | `0.002` | Convergence tolerance between successive quadrature refinements. |

With two regions the evidence is computed analytically and the quadrature
options are not used.

## INLA

| Option | Default | Description |
|--------|---------|-------------|
| `threads` | `"2:1"` | Passed to INLA as `num.threads`. |
| `strategy` | `"gaussian"` | INLA approximation strategy (`control.inla$strategy`). |
| `kappa_h` | `sqrt(8)` | Product of the SPDE scale and `h`; must lie in [2.3, 3.2]. |

Hyperparameters are always integrated with INLA's CCD design.

## Prior on the regional gene-effect variances

Each region *r* has a between-gene variance τ²<sub>r</sub>. The prior is set on
the log precision.

| Option | Default | Description |
|--------|---------|-------------|
| `tau_prior` | `"normal"` | One of `"normal"`, `"inla_default"`, `"pc"`, `"flat"`. |
| `tau_logprecision_mean` | `0` | Mean of the normal prior on the log precision. |
| `tau_logprecision_sd` | `3` | Standard deviation of the normal prior on the log precision. |
| `tau_pc_u`, `tau_pc_alpha` | `2`, `0.05` | Parameters of the PC prior, used when `tau_prior = "pc"`. |

`"inla_default"` uses INLA's log-gamma(1, 5e-5) prior on the precision and
`"flat"` a flat prior on the log precision.

## Regions and covariates

| Option | Default | Description |
|--------|---------|-------------|
| `baseline_region` | first label | Label of the baseline region. |
| `other_weights` | spot counts | Weights of the other regions in the mean contrast of the baseline region against the rest (`dmu_1_vs_rest`). Length R − 1, non-negative. |
| `gene_specific_covariates` | none | Names of columns of `X` that get an independent slope for each gene. |
| `gene_covariate_prior_sd` | `10` | Standard deviation of the fixed normal prior on the gene-specific slopes. |

`other_weights` only affects the mean contrast. The directional "overall"
results always weight the other regions by their spot counts.

## Fixed priors

The following priors are fixed in the code and are not options: a PC prior
with parameters (2, 0.05) on the Gaussian observation precision, N(0, 100)
priors on the fixed effects, and a N(0, 1) prior on the SPDE variance
parameter.
