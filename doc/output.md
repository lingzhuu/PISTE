# Output

`screen_summarize()` returns a list, and `screen_write()` saves its main parts
as CSV and RDS files.

## Directional results

`summary$pi` holds the directional results (it is `NULL` when
`cfg$pi = FALSE`).

For a pathway gene in region *r*, its state is modelled as
Z<sub>r</sub> ~ N(μ<sub>r</sub>, τ²<sub>r</sub>), where μ<sub>r</sub> is the
regional mean and τ²<sub>r</sub> the between-gene variance in that region. For
two regions, π<sub>rs</sub> = P(Z<sub>r</sub> > Z<sub>s</sub>).

| Column | Meaning |
|--------|---------|
| `score` | Posterior mean of π<sub>rs</sub>. 0.5 means no difference. |
| `evidence` | Posterior probability that π<sub>rs</sub> > 0.5, i.e. that the pathway is more active in *r* than in *s*. |
| `boundary_mass` | Posterior mass at exactly π<sub>rs</sub> = 0.5 (ties). Normally 0. |

`summary$pi$pairwise` has one row per ordered pair (`region`,
`reference_region`). `summary$pi$overall` has one row per region and compares it
with the spot-count-weighted average of all other regions; it also lists the
number of `spots`.

`summary$pi$diagnostics` records how the evidence was integrated: whether the
quadrature converged, the final number of nodes, the refinement history and
the number of CCD configurations. A warning is raised when the quadrature
reaches `pi_max_nodes` without converging.

## Parameters

`summary$parameters` has columns `parameter`, `mean`, `lower`, `median`,
`upper` (2.5 % and 97.5 % posterior quantiles). Rows:

| Parameter | Meaning |
|-----------|---------|
| `baseline_intercept` | Mean pathway level of the baseline region. |
| `relative_mu_r2`, ... | Offset of each other region from the baseline. |
| `cov_<name>` | Shared covariate effects. |
| `absolute_mu_r1`, ... | Mean pathway level of each region. |
| `tau2_r1`, ... | Between-gene variance in each region: how consistently the pathway's genes behave there. |
| `sigma_s2` | Variance of the spatial fields. |
| `sigma_e2` | Residual noise variance. |

Region indices follow `summary$metadata$region_labels`, with the baseline
region first.

## Mean contrasts and region calls

`summary$mean_contrasts` lists contrasts of the regional means μ<sub>r</sub>:
`dmu_<r>_vs_rest` (region *r* minus the weighted average of the others) and
`dmu_<r>_vs_<s>` for every pair. Columns: `quantity`, `mean`, `lower`,
`median`, `upper`, `p_gt_zero` (posterior probability that the contrast is
positive) and `lfsr` (local false sign rate, the smaller of `p_gt_zero` and
1 − `p_gt_zero`).

`summary$domains` takes the `dmu_<r>_vs_rest` rows, adds the region label and a
`direction`: `"higher"` when `p_gt_zero` ≥ 0.95, `"lower"` when it is ≤ 0.05,
otherwise `"undecided"`.

These tables concern the regional means only. The directional score and
evidence also account for the between-gene variances.

## Other elements

- `gene_covariate_coefficients`: per-gene slopes of the gene-specific
  covariates, or `NULL`.
- `fixed_marginals`, `hyper_marginals`, `contrast_marginals`: INLA posterior
  marginals.
- `metadata`: dimensions, region labels and spot counts, `h`, the options
  used, the INLA version, the model formula and fit time.
- `warnings`: warnings reported by INLA.

## Files written by `screen_write()`

| File | Content |
|------|---------|
| `screen_summary.rds` | The full summary list. Read with `readRDS()`. |
| `overall_scores.csv` | `summary$pi$overall` |
| `pairwise_scores.csv` | `summary$pi$pairwise` |
| `pi_diagnostics.rds`, `pi_refinement.csv` | Integration diagnostics. |
| `joint_marginal_diagnostic.csv` | Agreement between the CCD mixture and INLA's own marginals of the regional means. |
| `parameters.csv` | `summary$parameters` |
| `mean_contrasts.csv` | `summary$mean_contrasts` |
| `domains.csv` | `summary$domains` |
| `runtime.csv` | Fit and post-processing time in seconds. |
| `complete.txt` | Completion marker. |

The `pi` files are written only when directional inference is enabled.
