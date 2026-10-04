# API

```r
library(PathwayLGM)
```

SpaPath can be run in two ways. Quick mode runs the whole workflow for a
panel of gene sets in one call to `run_spapath()`. Step-by-step mode calls
the functions it is built from, for custom workflows: building the mesh once
with other settings, preparing the adjusted-mode inputs yourself, or fitting
a single gene set with other model options.

### Quick mode

| Function | Description |
|----------|-------------|
| [`run_spapath()`](#run_spapath) | Run SpaPath on every gene set of a panel, from raw counts to result tables. |
| [`reference_profile()`](#reference_profile) | Reference expression of each cell type from single-cell counts, for the adjusted mode. |

### Step-by-step mode

| Function | Description |
|----------|-------------|
| [`build_mesh()`](#build_mesh) | Build the mesh for the spatial fields of a section. |
| [`estimate_dispersion_scale()`](#estimate_dispersion_scale) | Relative negative binomial dispersion of every gene, for the adjusted mode. |
| [`composition_offset()`](#composition_offset) | Log expression predicted from cell-type proportions, for the adjusted mode. |
| [`screen_joint_fit()`](#screen_joint_fit) | Fit the spatial pathway model for one gene set. |
| [`screen_summarize()`](#screen_summarize) | Compute π scores, evidence and posterior summaries from a fit. |
| [`screen_write()`](#screen_write) | Write a summary to CSV and RDS files. |
| [`pathway_pi_integrate()`](#pathway_pi_integrate) | Recompute π scores and evidence from a fit. |
| [`pathway_pi()`](#pathway_pi) | π for given regional means and variances. |

---

## run_spapath

```r
run_spapath(counts, xy, region, gene_sets,
            mode = c("observed", "adjusted"), out = NULL,
            normalized = FALSE, proportions = NULL, reference = NULL,
            baseline = NULL, h = 0.3, min_genes = 5,
            mesh_args = list(), cfg = list(), verbose = TRUE)
```

Runs SpaPath on every gene set of a panel in one call:

1. In observed mode, normalizes the counts with SCTransform (sctransform v2
   Pearson residuals, clipped at ±√(N/30) as in Seurat). In adjusted mode,
   estimates the dispersion scales with `estimate_dispersion_scale()` and the
   composition offset with `composition_offset()`, and keeps the genes with a
   dispersion scale and positive predicted expression in every spot.
2. Rescales the coordinates to the unit square (shifted to start at zero and
   divided by their largest value) and builds the mesh with `build_mesh()`.
3. Fits each gene set with `screen_joint_fit()` and summarizes it with
   `screen_summarize()`. Gene sets with fewer than `min_genes` usable genes
   are skipped; a gene set whose fit fails is reported and the run goes on.

### Arguments

| Argument | Description |
|----------|-------------|
| `counts` | Spots x genes matrix (dense or sparse) of raw counts with spot rownames and gene colnames. With `normalized = TRUE` in observed mode, normalized expression instead. |
| `xy` | Spot coordinates with two columns, in the spot order of `counts`, in any unit. |
| `region` | Region label of each spot, in the spot order of `counts`, no `NA`. |
| `gene_sets` | Named list of gene vectors, or the path of a GMT file. |
| `mode` | `"observed"` (default) or `"adjusted"`. |
| `out` | Optional output folder; see Value. |
| `normalized` | Observed mode. `TRUE` when `counts` is already normalized, to skip SCTransform. |
| `proportions` | Adjusted mode. Spots x cell types matrix of cell-type proportions, e.g. normalized RCTD weights, with spot rownames. |
| `reference` | Adjusted mode. Cell types x genes matrix of reference mean expression on a linear scale, e.g. from `reference_profile()`. |
| `baseline` | Baseline region of the model. Default: the first label in `region`. |
| `h` | Spatial length scale in unit-square units; the practical range of the spatial fields with the default `kappa_h`. |
| `min_genes` | Minimum number of usable genes for a gene set to be fitted (at least 2). |
| `mesh_args` | List of arguments passed to `build_mesh()`, e.g. `list(cutoff = 0.05)`. |
| `cfg` | List of options passed to `screen_joint_fit()`; see [Options in `cfg`](#options-in-cfg). |
| `verbose` | Print progress messages. |

### Value

Invisibly, a list with

- `overall`: the `overall` π tables of all fitted gene sets, with columns
  `gene_set` and `genes` (number of genes used) added;
- `pairwise`: the `pairwise` π tables, likewise;
- `gene_sets`: one row per gene set with the genes in the set, the genes
  used, `status` (`done`, `skipped` or `failed`) and a message;
- `results`: the `screen_summarize()` output of each fitted gene set.

With `out`, each gene set is written to its own subfolder with
`screen_write()`, and `overall_scores.csv`, `pairwise_scores.csv` and
`gene_sets.csv` are written to `out`.

---

## reference_profile

```r
reference_profile(counts, cell_type)
```

Pools the raw counts of all cells of each cell type and scales them to
counts per 10,000. The result is the `reference` of `run_spapath()` and
`composition_offset()`.

### Arguments

| Argument | Description |
|----------|-------------|
| `counts` | Cells x genes matrix (dense or sparse) of raw single-cell counts with gene colnames. |
| `cell_type` | Cell type of each cell; cells with `NA` are dropped. |

### Value

Cell types x genes matrix.

---

## build_mesh

```r
build_mesh(xy, max_edge = 0.25, cutoff = 0.08, offset = c(0.1, 0.2))
```

Builds the 2D triangular mesh on which the spatial fields are defined, with
`fmesher::fm_mesh_2d()`. Build it once for a section and pass it to
`screen_joint_fit()` for every gene set. The defaults assume coordinates
rescaled to the unit square, shifted to start at zero and divided by their
largest value.

### Arguments

| Argument | Description |
|----------|-------------|
| `xy` | Numeric matrix of spot coordinates with two columns. |
| `max_edge` | Largest triangle edge. |
| `cutoff` | Minimum distance between mesh nodes; closer spots share a node. |
| `offset` | Inner and outer extension of the mesh around the spots. |

### Value

An `fm_mesh_2d` mesh; `mesh$n` is its number of nodes.

---

## screen_joint_fit

```r
screen_joint_fit(Y, xy, region, mesh, h, mode = c("observed", "adjusted"),
                 X = NULL, offset = NULL,
                 library_size = NULL, dispersion_scale = NULL,
                 composition_offset = NULL, cfg = list())
```

Fits the SpaPath model for one gene set with INLA, integrating over the
hyperparameters with a CCD design. The mode chooses the model:

- `"observed"` fits a Gaussian model to normalized expression and reports the
  total regional difference in pathway activity.
- `"adjusted"` fits a negative binomial model to raw counts with a
  cell-type composition offset and reports the regional difference that
  composition does not explain. For spot i and gene g,
  log E[y<sub>ig</sub>] = log L<sub>i</sub> + a<sub>g</sub> + μ<sub>r(i)</sub>
  + gene effect + spatial field + λ<sub>g</sub> log m<sub>ig</sub>, where
  L<sub>i</sub> is the library size, m<sub>ig</sub> the expression predicted
  from the spot's cell-type proportions, a<sub>g</sub> a gene intercept and
  λ<sub>g</sub> = 1 + v<sub>g</sub> with v<sub>g</sub> ~ N(0, 0.25²). Gene g
  has NB size `nb_size0` × S<sub>g</sub>.

### Arguments

| Argument | Description |
|----------|-------------|
| `Y` | Spots x genes matrix with at least two genes and no missing or infinite values. Normalized expression in observed mode; raw integer counts with unique spot rownames and gene colnames in adjusted mode. Column names are used as gene labels. |
| `xy` | Numeric matrix of spot coordinates with two columns, in the same spot order as `Y`. |
| `region` | Region label of each spot, at least two distinct labels, no `NA`. Labels are converted to character and ordered by first appearance. |
| `mesh` | Mesh for the spatial fields, from [`build_mesh()`](#build_mesh). |
| `h` | Positive spatial length scale in the units of `xy`. The SPDE scale is fixed at `cfg$kappa_h / h`; with the default `kappa_h`, `h` is the practical range of the spatial fields. |
| `mode` | `"observed"` (default) or `"adjusted"`. |
| `X` | Optional numeric covariate matrix with N rows and unique column names. An `intercept` column is added if absent, and no other column may be constant one. |
| `offset` | Optional known additive offset of length N, added to every gene of a spot. |
| `library_size` | Adjusted mode. Library size of each spot, e.g. total UMI over all genes. Only ratios matter; they are rescaled to a geometric mean of one. |
| `dispersion_scale` | Adjusted mode. Relative dispersion scale S<sub>g</sub> of each gene in `Y`, from `estimate_dispersion_scale()`. If named, names must match `colnames(Y)`. |
| `composition_offset` | Adjusted mode. Spots x genes matrix of log m<sub>ig</sub> from `composition_offset()`, finite, with the dimnames of `Y`. |
| `cfg` | List of options; see below. Entries not given keep their defaults. |

### Options in `cfg`

| Option | Default | Description |
|--------|---------|-------------|
| `baseline_region` | first label | Region whose mean is the model intercept. |
| `pi` | `TRUE` | Compute π scores and evidence. |
| `pi_nodes` | `4096` | Initial number of quadrature nodes for the evidence with three or more regions (power of two, ≥ 16). |
| `pi_max_nodes` | `65536` | Maximum number of quadrature nodes (power of two, ≥ 2 × `pi_nodes`). |
| `pi_tolerance` | `0.002` | Convergence tolerance between quadrature refinements. |
| `threads` | `"1:1"` | INLA `num.threads`, as `"A:B"`: A threads evaluate hyperparameter configurations in parallel and each evaluation uses B threads. |
| `strategy` | `"gaussian"`; `"simplified.laplace"` in adjusted mode | INLA approximation strategy. |
| `kappa_h` | `sqrt(8)` | SPDE scale times `h`; must lie in [2.3, 3.2]. |
| `tau_prior` | `"normal"` | Prior on the log precision of the regional gene effects: `"normal"`, `"inla_default"` (log-gamma(1, 5e-5) on the precision), `"pc"` or `"flat"`. |
| `tau_logprecision_mean` | `0` | Mean of the normal prior. |
| `tau_logprecision_sd` | `3` | Standard deviation of the normal prior. |
| `tau_pc_u`, `tau_pc_alpha` | `2`, `0.05` | Parameters of the PC prior. |
| `other_weights` | spot counts | Weights of the other regions in the mean contrast of the baseline region against the rest. Length R − 1, non-negative. |
| `gene_specific_covariates` | none | Columns of `X` that get an independent slope for each gene. Observed mode only. |
| `gene_covariate_prior_sd` | `10` | Standard deviation of the fixed normal prior on the gene-specific slopes. |
| `nb_logsize_mean`, `nb_logsize_sd` | `log(10)`, `1.5` | Normal prior on log `nb_size0`. Adjusted mode only. |
| `low_memory` | `FALSE` | AMD reordering with a small node factor: the same model with less memory, useful for large gene sets. |
| `initial_theta` | `NULL` | Starting hyperparameters (length R + 2), e.g. `fit$fit$mode$theta` from an earlier run. |
| `verbose` | `FALSE` | Print INLA's log. |

### Value

A list with

- `fit`: the INLA fit;
- `joint`: the posterior of the regional means and variances used for π, or
  `NULL` when `cfg$pi = FALSE`;
- `metadata`: dimensions, mode, region labels (baseline first) and spot
  counts, `h`, the options used, the INLA version, the model formula and the
  fit time. In adjusted mode also the dispersion scales and
  `composition_diagnostics`: per gene, the standard deviation of
  log m<sub>ig</sub> within regions (small values mean λ<sub>g</sub> is
  informed mainly by its prior) and its mean.

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
  spatial variance `sigma_s2` and noise variance `sigma_e2`. In adjusted mode
  the regional means are composition-adjusted states named
  `adjusted_state_r1`, …, and `sigma_e2` is replaced by the NB size
  `nb_size0` and its inverse `alpha0`.
- `mean_contrasts`: differences of regional means, each region against the
  rest (`dmu_<r>_vs_rest`) and every pair (`dmu_<r>_vs_<s>`), with
  `p_gt_zero` and `lfsr`.
- `domains`: each region called `"higher"`, `"lower"` or `"undecided"`
  against the rest from `p_gt_zero` ≥ 0.95 or ≤ 0.05 of its mean contrast.
- `gene_covariate_coefficients`: per-gene slopes, or `NULL`.
- `gene_coefficients`: adjusted mode only, the gene intercepts `a_g` and
  composition slopes `lambda_g` of each gene; `NULL` in observed mode.
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
and `complete.txt`. The π files are written only when π was computed. In
adjusted mode `gene_coefficients.csv` and `composition_diagnostics.csv` are
also written.

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

---

## composition_offset

```r
composition_offset(proportions, reference)
```

| Argument | Description |
|----------|-------------|
| `proportions` | Non-negative spots x cell types matrix of cell-type proportions, e.g. RCTD weights, with spot rownames. Rows are rescaled to sum to one. |
| `reference` | Non-negative cell types x genes matrix of reference mean expression on a linear scale, e.g. counts per 10,000. Row names must be the cell types in `colnames(proportions)`. |

Returns the spots x genes matrix log m<sub>ig</sub>, where
m<sub>ig</sub> = Σ<sub>t</sub> p<sub>it</sub> ν<sub>tg</sub> is the
expression predicted from the spot's composition. No pseudocount is added:
a gene with zero predicted expression in some spot gets `-Inf` there and must
be removed before fitting.

---

## estimate_dispersion_scale

```r
estimate_dispersion_scale(counts, library_size, min_total = 10L,
                          min_detected = 3L, fit_type = "parametric",
                          block_size = 256L)
```

| Argument | Description |
|----------|-------------|
| `counts` | Raw counts, spots x genes, dense or sparse, with unique spot and gene names. Use all genes, not only pathway genes. |
| `library_size` | Library size of each spot, e.g. total UMI. |
| `min_total`, `min_detected` | A gene is eligible when its total count is at least `min_total` and it is detected in at least `min_detected` spots. At least 30 genes must be eligible. |
| `fit_type` | DESeq2 dispersion trend type. |
| `block_size` | Genes per block, to bound memory. |

Estimates the NB dispersion α<sub>g</sub> of every eligible gene with DESeq2
(intercept-only design, trend fit and MAP shrinkage) and returns a data frame
with columns `gene`, `total`, `detected`, `eligible`, `dispersion` and
`dispersion_scale` = exp(mean(log α)) / α<sub>g</sub>, which is `NA` for
ineligible genes. Gene-wise and MAP steps run in blocks, while the trend and
its prior variance are fitted once on all genes, so the result equals an
unblocked DESeq2 run. Requires DESeq2.
