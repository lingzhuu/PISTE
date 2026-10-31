# API

```r
library(PathwayLGM)
```

SpaPath answers two questions about a gene set. **Region comparison** asks
whether the pathway is more active in one labelled region than in another
and reports a directional score π with its evidence; it needs region labels
and no coordinates. **Spatial map** asks where in the tissue the pathway is
more or less active; it needs spot coordinates and no region labels.

Each interface can be run in quick mode, which runs the whole workflow for a
panel of gene sets in one call, or step by step with the functions it is
built from.

## Region comparison

### Quick mode

| Function | Description |
|----------|-------------|
| [`compare_regions()`](#compare_regions) | Compare the regions on every gene set of a panel, from raw counts to result tables. |
| [`reference_profile()`](#reference_profile) | Reference expression of each cell type from single-cell counts, for the adjusted mode. |

### Step-by-step mode

| Function | Description |
|----------|-------------|
| [`estimate_dispersion_scale()`](#estimate_dispersion_scale) | Relative negative binomial dispersion of every gene, for the adjusted mode. |
| [`composition_offset()`](#composition_offset) | Log expression predicted from cell-type proportions, for the adjusted mode. |
| [`observed_fit()`](#observed_fit) | Fit the observed-mode model for one gene set, for one section or jointly for adjacent sections. |
| [`screen_joint_fit()`](#screen_joint_fit) | Fit the composition-adjusted model for one gene set with INLA. |
| [`screen_summarize()`](#screen_summarize) | Compute π scores, evidence and posterior summaries from a fit. |
| [`screen_write()`](#screen_write) | Write a summary to CSV and RDS files. |
| [`pathway_pi_integrate()`](#pathway_pi_integrate) | Recompute π scores and evidence from a fit. |
| [`pathway_pi()`](#pathway_pi) | π for given regional means and variances. |

## Spatial map

### Quick mode

| Function | Description |
|----------|-------------|
| [`map_activity()`](#map_activity) | Map the spatial activity of every gene set of a panel, from raw counts to regions of higher and lower activity. |

### Step-by-step mode

| Function | Description |
|----------|-------------|
| [`field_fit()`](#field_fit) | Fit the spatial field of one gene set with INLA and sample it at the spots. |
| [`field_regions()`](#field_regions) | Regions of higher and lower activity from a fit. |

---


## compare_regions

```r
compare_regions(counts, region, gene_sets,
                mode = c("observed", "adjusted"), out = NULL,
                normalized = FALSE, proportions = NULL, reference = NULL,
                baseline = NULL, section = NULL, joint_sections = FALSE,
                min_genes = 5, cfg = list(), verbose = TRUE)
```

Compares the regions on every gene set of a panel in one call:

1. In observed mode, normalizes the counts of each section with SCTransform
   (sctransform v2 Pearson residuals, clipped at ±√(N/30) as in Seurat). In
   adjusted mode, estimates the dispersion scales with
   `estimate_dispersion_scale()` and the composition offset with
   `composition_offset()`, and keeps the genes with a dispersion scale and
   positive predicted expression in every spot.
2. Fits each gene set with `observed_fit()` (observed mode) or
   `screen_joint_fit()` (adjusted mode) and summarizes it with
   `screen_summarize()`. Gene sets with fewer than `min_genes` usable genes
   are skipped; a gene set whose fit fails is reported and the run goes on.

Spot coordinates are not used.

### Arguments

| Argument | Description |
|----------|-------------|
| `counts` | Spots x genes matrix (dense or sparse) of raw counts with spot rownames and gene colnames. With `normalized = TRUE` in observed mode, normalized expression instead. |
| `region` | Region label of each spot, in the spot order of `counts`. Spots without a label are dropped. |
| `gene_sets` | Named list of gene vectors, or the path of a GMT file. |
| `mode` | `"observed"` (default) or `"adjusted"`. |
| `out` | Optional output folder; see Value. |
| `normalized` | Observed mode. `TRUE` when `counts` is already normalized, to skip SCTransform. |
| `proportions` | Adjusted mode. Spots x cell types matrix of cell-type proportions, e.g. normalized RCTD weights, with spot rownames. |
| `reference` | Adjusted mode. Cell types x genes matrix of reference mean expression on a linear scale, e.g. from `reference_profile()`. |
| `baseline` | Baseline region of the model. Default: the first label in `region`. |
| `section` | Optional section label of each spot. Sections are normalized separately and, by default, fitted one by one; the result tables then get a `section` column. |
| `joint_sections` | Observed mode. Fit the sections given by `section` jointly (see [`observed_fit()`](#observed_fit)). Use it only for adjacent sections of the same tissue block with the same region labels. |
| `min_genes` | Minimum number of usable genes for a gene set to be fitted (at least 2). |
| `cfg` | Options; see [Options in `cfg`](#options-in-cfg). |
| `verbose` | Print progress messages. |

### Value

Invisibly, a list with

- `overall`: the `overall` π tables of all fitted gene sets, with columns
  `gene_set` and `genes` (number of genes used) added, and `section` when
  sections are fitted one by one;
- `pairwise`: the `pairwise` π tables, likewise;
- `gene_sets`: one row per gene set with the genes in the set, the genes
  used, `status` (`done`, `skipped` or `failed`) and a message;
- `results`: the `screen_summarize()` output of each fitted gene set.

With `out`, each gene set is written to its own subfolder with
`screen_write()` (`<section>/<gene set>` when sections are fitted one by
one), and `overall_scores.csv`, `pairwise_scores.csv` and
`gene_sets.csv` are written to `out`.

---

## reference_profile

```r
reference_profile(counts, cell_type)
```

Pools the raw counts of all cells of each cell type and scales them to
counts per 10,000. The result is the `reference` of `compare_regions()` and
`composition_offset()`.

### Arguments

| Argument | Description |
|----------|-------------|
| `counts` | Cells x genes matrix (dense or sparse) of raw single-cell counts with gene colnames. |
| `cell_type` | Cell type of each cell; cells with `NA` are dropped. |

### Value

Cell types x genes matrix.

---

## observed_fit

```r
observed_fit(Y, region, section = NULL, cfg = list())
```

Fits the observed-mode model of one gene set: a Gaussian model of normalized
expression with a gene intercept, regional means and gene-level regional
effects,

y<sub>ig</sub> = α + u<sub>g</sub> + μ<sub>r(i)</sub> + γ<sub>g,r(i)</sub> + ε<sub>ig</sub>,
γ<sub>g,r</sub> ~ N(0, τ<sub>r</sub><sup>2</sup>), ε<sub>ig</sub> ~ N(0, σ<sup>2</sup>).

With `section`, adjacent sections of the same tissue block are fitted
jointly. Each section has its own intercept, gene intercepts and noise
variance; the regional means μ<sub>r</sub> and the gene-level regional effects
γ<sub>g,r</sub> are shared, and each section adds its own deviation
η<sub>g,r,k</sub> ~ N(0, ν<sub>r</sub><sup>2</sup>). π is computed from the
shared part, so it measures the regional difference reproduced across
sections. The sections are normalized separately and need not be aligned.

Given the variances, the posterior of the regional means is Gaussian and
exact. The variances are integrated numerically on a grid around their
posterior mode. A fit takes about 0.1 s for one section and a few seconds for
two sections.

### Arguments

| Argument | Description |
|----------|-------------|
| `Y` | Spots x genes matrix of normalized expression. Column names are used as gene labels. |
| `region` | Region label of each spot, at least two regions. |
| `section` | Optional section label of each spot. Every region must be present in every section. |
| `cfg` | Options; see [Options in `cfg`](#options-in-cfg). |

### Options in `cfg`

| Option | Default | Description |
|--------|---------|-------------|
| `pi_alpha` | `0.05` | A region is called higher or lower when the evidence in that direction is at least 1 − `pi_alpha`/2. |
| `pi_report_p_value` | `FALSE` | Add the column `p_value` = 2 min(P(π > 0.5), P(π < 0.5)). It is a posterior probability, not a frequentist p-value. |
| `pi` | `TRUE` | `FALSE` skips the π scores and reports the regional means only. |
| `baseline_region` | first label | Region whose mean is the model intercept. It does not change π. |
| `threads` | `"1:1"` | Adjusted mode only: INLA `num.threads`. |

Other names are refused.

### Value

A list with `joint` (the posterior of the regional means and variances used
for π, or `NULL` when `cfg$pi = FALSE`), `metadata` (dimensions, region and
section labels, spot counts per section and region, the options used and the
fit time) and `grid` (the grid points, their conditional posteriors and
weights). Pass it to [`screen_summarize()`](#screen_summarize).

---

## screen_joint_fit

```r
screen_joint_fit(Y, region, mode = "adjusted", X = NULL, offset = NULL,
                 library_size = NULL, dispersion_scale = NULL,
                 composition_offset = NULL, cfg = list())
```

Fits the composition-adjusted model of one gene set with INLA: a negative
binomial model of raw counts with a cell-type composition offset, reporting
the regional difference that composition does not explain. For spot i and
gene g, log E[y<sub>ig</sub>] = log L<sub>i</sub> + a<sub>g</sub> +
μ<sub>r(i)</sub> + gene effect + λ<sub>g</sub> log m<sub>ig</sub>, where
L<sub>i</sub> is the library size, m<sub>ig</sub> the expression predicted
from the spot's cell-type proportions, a<sub>g</sub> a gene intercept and
λ<sub>g</sub> = 1 + v<sub>g</sub> with v<sub>g</sub> ~ N(0, 0.25²). Gene g
has NB size `nb_size0` × S<sub>g</sub>. The observed mode is fitted by
[`observed_fit()`](#observed_fit).

### Arguments

| Argument | Description |
|----------|-------------|
| `Y` | Spots x genes matrix of raw integer counts with unique spot rownames and gene colnames, at least two genes. |
| `region` | Region label of each spot, at least two distinct labels. |
| `X` | Optional numeric covariate matrix with N rows and unique column names. An `intercept` column is added if absent. |
| `offset` | Optional known additive offset of length N. |
| `library_size` | Library size of each spot, e.g. total UMI over all genes. Only ratios matter. |
| `dispersion_scale` | Relative dispersion scale S<sub>g</sub> of each gene in `Y`, from `estimate_dispersion_scale()`. |
| `composition_offset` | Spots x genes matrix of log m<sub>ig</sub> from `composition_offset()`, with the dimnames of `Y`. |
| `cfg` | Options; see [Options in `cfg`](#options-in-cfg). |

### Value

A list with `fit` (the INLA fit), `joint` and `metadata` as for
[`observed_fit()`](#observed_fit); `metadata` also holds the dispersion scales
and `composition_diagnostics`: per gene, the standard deviation of
log m<sub>ig</sub> within regions (small values mean λ<sub>g</sub> is informed
mainly by its prior) and its mean.

---

## screen_summarize

```r
screen_summarize(result)
```

| Argument | Description |
|----------|-------------|
| `result` | Output of `observed_fit()` or `screen_joint_fit()`. |

### Value

A list with

- `pi$overall`: one row per region, compared with the spot-weighted average of
  all other regions. Columns `region`, `spots`, `score`, `evidence`,
  `evidence_higher`, `evidence_lower`, `call`, `boundary_mass`, and `p_value`
  when `pi_report_p_value = TRUE`.
- `pi$pairwise`: one row per ordered pair of regions. Columns `region`,
  `reference_region`, `score`, `evidence`, `evidence_higher`,
  `evidence_lower`, `call`, `boundary_mass` (and `p_value`).
- `pi$diagnostics`: quadrature convergence, number of nodes and refinement
  history.
- `parameters`: posterior mean, 2.5 %, 50 % and 97.5 % quantiles of
  `baseline_intercept`, `relative_mu_r2`, …, regional means
  `absolute_mu_r1`, …, between-gene variances `tau2_r1`, … and noise variance
  `sigma_e2`. For a joint fit of sections also the section offsets
  `section_offset_2`, …, the section-specific variances `nu2_r1`, … and one
  noise variance per section `sigma_e2_section_1`, …. In adjusted mode the
  regional means are composition-adjusted states named `adjusted_state_r1`,
  …, and `sigma_e2` is replaced by the NB size `nb_size0` and its inverse
  `alpha0`.
- `mean_contrasts`: differences of regional means, each region against the
  rest (`dmu_<r>_vs_rest`) and every pair (`dmu_<r>_vs_<s>`), with
  `p_gt_zero` and `lfsr`.
- `domains`: each region called `"higher"`, `"lower"` or `"undecided"`
  against the rest from `p_gt_zero` ≥ 0.95 or ≤ 0.05 of its mean contrast.
- Adjusted mode only: `gene_coefficients` (gene intercepts `a_g` and
  composition slopes `lambda_g`) and the INLA marginals.
- `metadata`.

`score` is the posterior mean of π, the probability that a gene's expression
state in the region exceeds that in the reference, with the measurement error
of a gene's regional means added to the gene-to-gene variance; 0.5 means no
difference.
`evidence_higher` = P(π > 0.5 | Y) and `evidence_lower` = P(π < 0.5 | Y) are
the posterior probabilities of higher and lower activity; `evidence` equals
`evidence_higher`. `call` is `"higher"` or `"lower"` when that
probability is at least 1 − `pi_alpha`/2 (0.975 by default), otherwise
`"none"`. The optional `p_value` = 2 min(P(π > 0.5), P(π < 0.5)) gives the
same calls at `p_value` ≤ `pi_alpha`; it is a posterior probability, not a
frequentist p-value, and is reported only for comparison with methods that
report p-values. Region indices in
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
`pi_diagnostics.rds`, `pi_refinement.csv` and `complete.txt`
(`joint_marginal_diagnostic.csv` in adjusted mode). The π files are written only when π was computed. In
adjusted mode `gene_coefficients.csv` and `composition_diagnostics.csv` are
also written.

---

## pathway_pi_integrate

```r
pathway_pi_integrate(joint, alpha = 0.05, report_p_value = FALSE,
                     nodes = 4096, max_nodes = 65536, tolerance = 0.002)
```

| Argument | Description |
|----------|-------------|
| `joint` | The `joint` element of `observed_fit()` or `screen_joint_fit()`. |
| `alpha` | A region is called when the evidence in one direction is at least 1 − `alpha`/2. |
| `report_p_value` | Add the `p_value` column. |
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
pathway_pi(mu_r, mu_s, tau2_r, tau2_s, kappa2 = 0)
```

| Argument | Description |
|----------|-------------|
| `mu_r`, `mu_s` | Regional means. |
| `tau2_r`, `tau2_s` | Regional between-gene variances. |
| `kappa2` | Measurement error of a gene's regional means, e.g. σ<sub>ε</sub><sup>2</sup>(1/n<sub>r</sub> + 1/n<sub>s</sub>) in the Gaussian model. |

Returns π = Φ((`mu_r` − `mu_s`) / √(`tau2_r` + `tau2_s` + `kappa2`)).
With `kappa2 = 0` this is P(Z<sub>r</sub> > Z<sub>s</sub>) for independent
Z<sub>r</sub> ~ N(`mu_r`, `tau2_r`) and Z<sub>s</sub> ~ N(`mu_s`, `tau2_s`).
Fits of the Gaussian model use σ<sub>ε</sub><sup>2</sup> at each point of
the hyperparameter grid. In the negative binomial model each count carries the
information w<sub>ig</sub> = 1/(1/m<sub>ig</sub> + 1/size<sub>g</sub>), where
m<sub>ig</sub> is its expected count without regional effects and
size<sub>g</sub> = nb_size0 · S<sub>g</sub> in that configuration; gene g has
κ<sub>g</sub><sup>2</sup> = 1/W<sub>gr</sub> + 1/W<sub>gs</sub>, with
W<sub>gr</sub> the sum of w<sub>ig</sub> over the spots of region r, and
`kappa2` is the harmonic mean of κ<sub>g</sub><sup>2</sup> over genes.
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

---

## map_activity

```r
map_activity(counts, coords, gene_sets, gene_field = TRUE, d = 0.1,
             alpha = 0.05, out = NULL, normalized = FALSE, min_genes = 5,
             threads = "8:1", verbose = TRUE)
```

Normalizes the counts with SCTransform (sctransform v2 Pearson residuals,
clipped at ±√(N/30) as in Seurat), fits each gene set with
[`field_fit()`](#field_fit) and calls its regions with
[`field_regions()`](#field_regions). Gene sets with fewer than `min_genes`
usable genes are skipped; a gene set whose fit fails is reported and the run
goes on. Region labels are not used.

### Arguments

| Argument | Description |
|----------|-------------|
| `counts` | Spots x genes matrix (dense or sparse) of raw counts with spot rownames and gene colnames. With `normalized = TRUE`, normalized expression instead. |
| `coords` | Spots x 2 matrix of spot coordinates, in the spot order of `counts`. |
| `gene_sets` | Named list of gene vectors, or the path of a GMT file. |
| `gene_field` | Add a spatial field of each gene; see [`field_fit()`](#field_fit). `FALSE` is about 50 times faster but less accurate. |
| `d`, `alpha` | Minimum effect and joint posterior level of the regions; see [`field_regions()`](#field_regions). |
| `out` | Optional output folder; see Value. |
| `normalized` | `TRUE` when `counts` is already normalized, to skip SCTransform. |
| `min_genes` | Minimum number of usable genes for a gene set to be fitted (at least 2). |
| `threads` | INLA `num.threads`. |
| `verbose` | Print progress messages. |

### Value

Invisibly, a list with

- `gene_sets`: one row per gene set with the genes in the set, the genes
  used, `status` (`done`, `skipped` or `failed`), a message, the posterior
  mean of the sd of z, the number of spots called higher and lower, and the
  fit time;
- `results`: per fitted gene set, the [`field_regions()`](#field_regions)
  table (`spots`) and the hyperparameters.

With `out`, each gene set is written to its own subfolder (`spots.csv`,
`hyperparameters.csv`) and `gene_sets.csv` is written to `out`.

---

## field_fit

```r
field_fit(Y, coords, gene_field = TRUE, threads = "8:1")
```

Fits one spatial field shared by all genes of the gene set,

y<sub>ig</sub> = a<sub>g</sub> + z(s<sub>i</sub>) + w<sub>g</sub>(s<sub>i</sub>) + ε<sub>ig</sub>,
ε<sub>ig</sub> ~ N(0, σ<sub>ε</sub><sup>2</sup>),

where a<sub>g</sub> is a gene intercept, z is a Matérn field on a triangular
mesh of the spots (the pathway map, constrained to sum to zero) and
w<sub>g</sub> is a Matérn field of gene g on a coarser mesh, with one range
and one sd for all genes, which absorbs the spatial pattern a gene has
outside the pathway. Every gene enters z with weight 1, so genes that change
in opposite directions cancel. Both fields have penalized complexity priors
P(range < 0.1) = 0.5 and P(sd > 1) = 0.05 on coordinates scaled so that the
longer side of the section is 1, and the noise precision has the prior
P(σ<sub>ε</sub> > 2) = 0.05. The model is fitted with INLA and z is sampled
1,000 times jointly at the spots.

### Arguments

| Argument | Description |
|----------|-------------|
| `Y` | Spots x genes matrix of normalized expression with spot rownames, at least two genes. |
| `coords` | Spots x 2 matrix of spot coordinates, in the spot order of `Y`. |
| `gene_field` | Add the gene fields w<sub>g</sub>. `FALSE` is about 50 times faster but less accurate: spatial patterns of single genes that the gene set does not share can then appear as regions of the pathway. |
| `threads` | INLA `num.threads`. |

### Value

A list with `draws` (spots x 1,000 matrix of z), `sigma_e` (posterior mean
of σ<sub>ε</sub> from the noise precision), `hyperparameters` (mean, 2.5 %,
50 % and 97.5 % quantiles of `noise_precision`, `range_z`, `sd_z` and, with
gene fields, `range_w` and `sd_w`; ranges on the scaled coordinates) and
`metadata`. Pass it to [`field_regions()`](#field_regions).

---

## field_regions

```r
field_regions(fit, d = 0.1, alpha = 0.05)
```

Calls the regions where z exceeds the minimum effect u = `d` σ<sub>ε</sub>.
A spot is called `higher` when it belongs to the largest set of spots on
which z > u everywhere with joint posterior probability at least
1 − `alpha`, and `lower` likewise for z < −u. The spots enter the set in
decreasing order of P(z > u). Because the probability is joint, the region as
a whole has posterior probability 1 − `alpha`, not each spot on its own.

| Argument | Description |
|----------|-------------|
| `fit` | Output of [`field_fit()`](#field_fit). |
| `d` | Minimum effect in units of σ<sub>ε</sub>. `0` calls any spot where z is credibly above or below zero, however small. |
| `alpha` | One minus the joint posterior probability of a region. |

Returns a data frame with one row per spot: `spot`, `z_mean`, `z_sd`,
`p_higher` = P(z > u), `p_lower` = P(z < −u) and `call` (`higher`, `lower`
or `none`).

