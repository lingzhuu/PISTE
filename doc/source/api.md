# API

```r
library(PathwayLGM)
```

| Function | Description |
|----------|-------------|
| [`screen_joint_fit()`](#screen_joint_fit) | Fit the spatial pathway model for one gene set. |
| [`screen_summarize()`](#screen_summarize) | Compute π scores, evidence and posterior summaries from a fit. |
| [`screen_write()`](#screen_write) | Write a summary to CSV and RDS files. |
| [`pathway_pi_integrate()`](#pathway_pi_integrate) | Recompute π scores and evidence from a fit. |
| [`pathway_pi()`](#pathway_pi) | π for given regional means and variances. |

---

## screen_joint_fit

```r
screen_joint_fit(Y, xy, region, mesh, h, family = "gaussian",
                 X = NULL, offset = NULL, cfg = list())
```

Fits the SpaPath model for one gene set with INLA, integrating over the
hyperparameters with a CCD design.

### Arguments

| Argument | Description |
|----------|-------------|
| `Y` | Numeric matrix of normalized expression, spots x genes, with at least two genes and no missing or infinite values. Column names are used as gene labels. |
| `xy` | Numeric matrix of spot coordinates with two columns, in the same spot order as `Y`. |
| `region` | Region label of each spot, at least two distinct labels, no `NA`. Labels are converted to character and ordered by first appearance. |
| `mesh` | 2D mesh for the spatial fields, e.g. from `fmesher::fm_mesh_2d_inla()`. |
| `h` | Positive spatial length scale in the units of `xy`. The SPDE scale is fixed at `cfg$kappa_h / h`; with the default `kappa_h`, `h` is the practical range of the spatial fields. |
| `family` | Likelihood. Only `"gaussian"` is supported. |
| `X` | Optional numeric covariate matrix with N rows and unique column names. An `intercept` column is added if absent, and no other column may be constant one. |
| `offset` | Optional known additive offset of length N, added to every gene of a spot. |
| `cfg` | List of options; see below. Entries not given keep their defaults. |

### Options in `cfg`

| Option | Default | Description |
|--------|---------|-------------|
| `baseline_region` | first label | Region whose mean is the model intercept. |
| `pi` | `TRUE` | Compute π scores and evidence. |
| `pi_nodes` | `4096` | Initial number of quadrature nodes for the evidence with three or more regions (power of two, ≥ 16). |
| `pi_max_nodes` | `65536` | Maximum number of quadrature nodes (power of two, ≥ 2 × `pi_nodes`). |
| `pi_tolerance` | `0.002` | Convergence tolerance between quadrature refinements. |
| `threads` | `"2:1"` | INLA `num.threads`. |
| `strategy` | `"gaussian"` | INLA approximation strategy. |
| `kappa_h` | `sqrt(8)` | SPDE scale times `h`; must lie in [2.3, 3.2]. |
| `tau_prior` | `"normal"` | Prior on the log precision of the regional gene effects: `"normal"`, `"inla_default"` (log-gamma(1, 5e-5) on the precision), `"pc"` or `"flat"`. |
| `tau_logprecision_mean` | `0` | Mean of the normal prior. |
| `tau_logprecision_sd` | `3` | Standard deviation of the normal prior. |
| `tau_pc_u`, `tau_pc_alpha` | `2`, `0.05` | Parameters of the PC prior. |
| `other_weights` | spot counts | Weights of the other regions in the mean contrast of the baseline region against the rest. Length R − 1, non-negative. |
| `gene_specific_covariates` | none | Columns of `X` that get an independent slope for each gene. |
| `gene_covariate_prior_sd` | `10` | Standard deviation of the fixed normal prior on the gene-specific slopes. |

### Value

A list with

- `fit`: the INLA fit;
- `joint`: the posterior of the regional means and variances used for π, or
  `NULL` when `cfg$pi = FALSE`;
- `metadata`: dimensions, region labels (baseline first) and spot counts, `h`,
  the options used, the INLA version, the model formula and the fit time.

---

## screen_summarize

```r
screen_summarize(result)
```

| Argument | Description |
|----------|-------------|
| `result` | Output of `screen_joint_fit()`. |

### Value

A list with

- `pi$overall`: one row per region, compared with the spot-weighted average of
  all other regions. Columns `region`, `spots`, `score`, `evidence`,
  `boundary_mass`.
- `pi$pairwise`: one row per ordered pair of regions. Columns `region`,
  `reference_region`, `score`, `evidence`, `boundary_mass`.
- `pi$diagnostics`: quadrature convergence, number of nodes and refinement
  history.
- `parameters`: posterior mean, 2.5 %, 50 % and 97.5 % quantiles of
  `baseline_intercept`, `relative_mu_r2`, …, covariate effects `cov_<name>`,
  regional means `absolute_mu_r1`, …, between-gene variances `tau2_r1`, …,
  spatial variance `sigma_s2` and noise variance `sigma_e2`.
- `mean_contrasts`: differences of regional means, each region against the
  rest (`dmu_<r>_vs_rest`) and every pair (`dmu_<r>_vs_<s>`), with
  `p_gt_zero` and `lfsr`.
- `domains`: each region called `"higher"`, `"lower"` or `"undecided"`
  against the rest from `p_gt_zero` ≥ 0.95 or ≤ 0.05 of its mean contrast.
- `gene_covariate_coefficients`: per-gene slopes, or `NULL`.
- `fixed_marginals`, `hyper_marginals`, `contrast_marginals`, `metadata`,
  `warnings`.

`score` is the posterior mean of π, the probability that a gene's expression
state in the region exceeds that in the reference; 0.5 means no difference.
`evidence` is the posterior probability that π > 0.5. Region indices in
parameter and contrast names follow `metadata$region_labels`.

---

## screen_write

```r
screen_write(summary, out)
```

| Argument | Description |
|----------|-------------|
| `summary` | Output of `screen_summarize()`. |
| `out` | Output directory, created if needed. |

Writes `screen_summary.rds`, `overall_scores.csv`, `pairwise_scores.csv`,
`parameters.csv`, `mean_contrasts.csv`, `domains.csv`, `runtime.csv`,
`pi_diagnostics.rds`, `pi_refinement.csv`, `joint_marginal_diagnostic.csv`
and `complete.txt`. The π files are written only when π was computed.

---

## pathway_pi_integrate

```r
pathway_pi_integrate(joint, nodes = 4096L, max_nodes = 65536L,
                     tolerance = 0.002)
```

| Argument | Description |
|----------|-------------|
| `joint` | The `joint` element of `screen_joint_fit()`. |
| `nodes` | Initial number of quadrature nodes (power of two, ≥ 16). |
| `max_nodes` | Maximum number of nodes (power of two, ≥ 2 × `nodes`). |
| `tolerance` | Convergence tolerance between refinements. |

Returns a list with `pairwise`, `overall` and `diagnostics`, as in
`screen_summarize()`. With two regions the evidence is computed analytically;
with more regions it uses deterministic quasi-Monte Carlo quadrature and warns
if `max_nodes` is reached without convergence.

---

## pathway_pi

```r
pathway_pi(mu_r, mu_s, tau2_r, tau2_s)
```

| Argument | Description |
|----------|-------------|
| `mu_r`, `mu_s` | Regional means. |
| `tau2_r`, `tau2_s` | Regional between-gene variances. |

Returns π = P(Z<sub>r</sub> > Z<sub>s</sub>) for independent
Z<sub>r</sub> ~ N(`mu_r`, `tau2_r`) and Z<sub>s</sub> ~ N(`mu_s`, `tau2_s`).
Vectorised over its arguments.
