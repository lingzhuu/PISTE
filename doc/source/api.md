# R API

```r
library(PISTE)
```

The region comparison needs region labels; the spatial map needs spot
coordinates. [`compare()`](#compare) and [`map_activity()`](#map_activity) run
every gene set in one call. The other functions are the steps they are built
from, for when you want to run them yourself.

| Function | Analysis | What it does |
|----------|----------|--------------|
| [`compare()`](#compare) | Region comparison | Compares the regions on every gene set. |
| [`reference_profile()`](#reference_profile) | Adjusted mode | Reference expression of each cell type from single-cell counts. |
| [`estimate_dispersion_scale()`](#estimate_dispersion_scale) | Adjusted mode | Relative negative binomial dispersion of every gene. |
| [`composition_offset()`](#composition_offset) | Adjusted mode | Log expression predicted from cell-type proportions. |
| [`observed_fit()`](#observed_fit) | Region comparison | Fits the observed mode for one gene set, for one section or jointly for adjacent sections. |
| [`screen_joint_fit()`](#screen_joint_fit) | Adjusted mode | Fits the adjusted mode for one gene set with INLA. |
| [`screen_summarize()`](#screen_summarize) | Region comparison | π scores, evidence and posterior summaries of a fit. |
| [`screen_write()`](#screen_write) | Region comparison | Writes a summary to CSV and RDS files. |
| [`pathway_pi_integrate()`](#pathway_pi_integrate) | Region comparison | Recomputes π scores and evidence from a fit. |
| [`pathway_pi()`](#pathway_pi) | Region comparison | π for given regional means and variances. |
| [`map_activity()`](#map_activity) | Spatial map | Maps the activity of every gene set. |
| [`field_fit()`](#field_fit) | Spatial map | Fits the spatial field of one gene set with INLA. |
| [`field_regions()`](#field_regions) | Spatial map | Regions of higher and lower activity from a fit. |

## compare()

```r
compare(counts, region, gene_sets,
        mode = c("observed", "adjusted"), out = NULL,
        proportions = NULL, reference = NULL,
        baseline = NULL, section = NULL, joint_sections = FALSE,
        min_genes = 5, cfg = list(), verbose = TRUE)
```

Compares the regions on every gene set. In the adjusted mode it first
estimates the dispersion scales with `estimate_dispersion_scale()` and the
composition offset with `composition_offset()`, and keeps the genes that have
a dispersion scale and positive predicted expression in every spot. Each gene
set is then fitted with `observed_fit()` or `screen_joint_fit()` and
summarized with `screen_summarize()`. Gene sets with fewer than `min_genes`
usable genes are skipped; a gene set whose fit fails is reported and the run
goes on. Spot coordinates are not used.

| Argument | Description |
|----------|-------------|
| `counts` | Spots × genes matrix (dense or sparse) with spot row names and gene column names: normalized expression in the observed mode, raw counts in the adjusted mode. |
| `region` | Region label of each spot, in the spot order of `counts`. Spots without a label are dropped. |
| `gene_sets` | Named list of gene vectors, or the path of a GMT file. |
| `mode` | `"observed"` (default) or `"adjusted"`. |
| `out` | Optional output folder. |
| `proportions` | Adjusted mode. Spots × cell types matrix of cell-type proportions, for example normalized RCTD weights, with spot row names. |
| `reference` | Adjusted mode. Cell types × genes matrix of reference mean expression on a linear scale, for example from `reference_profile()`. |
| `baseline` | Baseline region of the model. Default: the first label in `region`. |
| `section` | Optional section label of each spot. By default each section is fitted on its own, and the result tables get a `section` column. |
| `joint_sections` | Observed mode. Fits the sections in `section` jointly (see [`observed_fit()`](#observed_fit)). Use it only for adjacent sections of the same tissue block with the same region labels. |
| `min_genes` | Minimum number of usable genes for a gene set to be fitted, at least 2. |
| `cfg` | Options; see [options in `cfg`](#options-in-cfg). |
| `verbose` | Print progress messages. |

**Returns**, invisibly, a list with `overall` (the `overall` π tables of all
fitted gene sets, with columns `gene_set`, `genes` (number of genes used) and,
for sections fitted on their own, `section`), `pairwise` (the `pairwise` π
tables, likewise), `gene_sets` (one row per gene set: genes in the set, genes
used, `status` `done`, `skipped` or `failed`, and a message) and `results`
(the `screen_summarize()` output of each fitted gene set). With `out`, each
gene set is written to its own subfolder with `screen_write()`
(`<section>/<gene set>` for sections fitted on their own), and
`overall_scores.csv`, `pairwise_scores.csv` and `gene_sets.csv` are written
to `out`.

## reference_profile()

```r
reference_profile(counts, cell_type)
```

Pools the raw counts of all cells of each cell type and scales them to counts
per 10,000. The result is the `reference` of `compare()` and
`composition_offset()`.

| Argument | Description |
|----------|-------------|
| `counts` | Cells × genes matrix (dense or sparse) of raw single-cell counts with gene column names. |
| `cell_type` | Cell type of each cell; cells with `NA` are dropped. |

**Returns** a cell types × genes matrix.

## estimate_dispersion_scale()

```r
estimate_dispersion_scale(counts, library_size, min_total = 10L,
                          min_detected = 3L, fit_type = "parametric",
                          block_size = 256L)
```

Estimates the negative binomial dispersion α<sub>g</sub> of every eligible
gene with DESeq2 (intercept-only design, trend fit and MAP shrinkage). The
gene-wise and MAP steps run in blocks, while the trend and its prior variance
are fitted once on all genes, so the result equals an unblocked DESeq2 run.
Requires DESeq2.

| Argument | Description |
|----------|-------------|
| `counts` | Raw counts, spots × genes, dense or sparse, with unique spot and gene names. Use all genes, not only gene set genes. |
| `library_size` | Library size of each spot, for example total UMI. |
| `min_total`, `min_detected` | A gene is eligible when its total count is at least `min_total` and it is detected in at least `min_detected` spots. At least 30 genes must be eligible. |
| `fit_type` | DESeq2 dispersion trend type. |
| `block_size` | Genes per block, to bound memory. |

**Returns** a data frame with columns `gene`, `total`, `detected`,
`eligible`, `dispersion` and `dispersion_scale` = exp(mean(log α)) / α<sub>g</sub>,
which is `NA` for genes that are not eligible.

## composition_offset()

```r
composition_offset(proportions, reference)
```

| Argument | Description |
|----------|-------------|
| `proportions` | Non-negative spots × cell types matrix of cell-type proportions, for example RCTD weights, with spot row names. Rows are rescaled to sum to one. |
| `reference` | Non-negative cell types × genes matrix of reference mean expression on a linear scale, for example counts per 10,000. Row names must be the cell types in `colnames(proportions)`. |

**Returns** the spots × genes matrix log m<sub>ig</sub>, where
m<sub>ig</sub> = Σ<sub>t</sub> p<sub>it</sub> ν<sub>tg</sub> is the expression
predicted from the spot's composition. No pseudocount is added: a gene with
zero predicted expression in some spot gets `-Inf` there and must be removed
before fitting.

## observed_fit()

```r
observed_fit(Y, region, section = NULL, cfg = list())
```

Fits the observed mode for one gene set: a Gaussian model of normalized
expression with a gene intercept, regional means and gene-level regional
effects,

y<sub>ig</sub> = α + u<sub>g</sub> + μ<sub>r(i)</sub> + γ<sub>g,r(i)</sub> + ε<sub>ig</sub>,
γ<sub>g,r</sub> ~ N(0, τ<sub>r</sub><sup>2</sup>), ε<sub>ig</sub> ~ N(0, σ<sup>2</sup>).

With `section`, adjacent sections of the same tissue block are fitted
jointly. Each section has its own intercept, gene intercepts and noise
variance; the regional means μ<sub>r</sub> and the gene-level regional effects
γ<sub>g,r</sub> are shared, and each section adds its own deviation
η<sub>g,r,k</sub> ~ N(0, ν<sub>r</sub><sup>2</sup>). π is computed from the
shared part, so it measures the regional difference the sections reproduce.
The sections need not be aligned.

Given the variances, the posterior of the regional means is Gaussian and
exact; the variances are integrated numerically on a grid around their
posterior mode. A fit takes about 0.1 s for one section and a few seconds for
two.

| Argument | Description |
|----------|-------------|
| `Y` | Spots × genes matrix of normalized expression. Column names are used as gene labels. |
| `region` | Region label of each spot, at least two regions. |
| `section` | Optional section label of each spot. Every region must be present in every section. |
| `cfg` | Options; see below. |

**Options in `cfg`**. Other names are refused.
{ #options-in-cfg }

| Option | Default | Description |
|--------|---------|-------------|
| `pi_alpha` | `0.05` | A region is called higher or lower when the evidence in that direction is at least 1 − `pi_alpha`/2. |
| `pi_report_p_value` | `FALSE` | Adds the column `p_value` = 2 min(P(π > 0.5), P(π < 0.5)). It is a posterior probability, not a frequentist p-value. |
| `pi` | `TRUE` | `FALSE` skips the π scores and reports the regional means only. |
| `baseline_region` | first label | Region whose mean is the model intercept. It does not change π. |
| `threads` | `"1:1"` | Adjusted mode only: INLA `num.threads`. |

**Returns** a list with `joint` (the posterior of the regional means and
variances used for π, or `NULL` when `cfg$pi = FALSE`), `metadata`
(dimensions, region and section labels, spot counts per section and region,
the options used and the fit time) and `grid` (the grid points, their
conditional posteriors and weights). Pass it to
[`screen_summarize()`](#screen_summarize).

## screen_joint_fit()

```r
screen_joint_fit(Y, region, mode = "adjusted", X = NULL, offset = NULL,
                 library_size = NULL, dispersion_scale = NULL,
                 composition_offset = NULL, cfg = list())
```

Fits the adjusted mode for one gene set with INLA: a negative binomial model
of raw counts with a cell-type composition offset, which reports the regional
difference that composition does not explain. For spot i and gene g,
log E[y<sub>ig</sub>] = log L<sub>i</sub> + a<sub>g</sub> +
μ<sub>r(i)</sub> + gene effect + λ<sub>g</sub> log m<sub>ig</sub>, where
L<sub>i</sub> is the library size, m<sub>ig</sub> the expression predicted
from the spot's cell-type proportions, a<sub>g</sub> a gene intercept and
λ<sub>g</sub> = 1 + v<sub>g</sub> with v<sub>g</sub> ~ N(0, 0.25²). Gene g
has negative binomial size `nb_size0` × S<sub>g</sub>.

| Argument | Description |
|----------|-------------|
| `Y` | Spots × genes matrix of raw integer counts with unique spot row names and gene column names, at least two genes. |
| `region` | Region label of each spot, at least two distinct labels. |
| `X` | Optional numeric covariate matrix with N rows and unique column names. An `intercept` column is added if absent. |
| `offset` | Optional known additive offset of length N. |
| `library_size` | Library size of each spot, for example total UMI over all genes. Only ratios matter. |
| `dispersion_scale` | Relative dispersion scale S<sub>g</sub> of each gene in `Y`, from `estimate_dispersion_scale()`. |
| `composition_offset` | Spots × genes matrix of log m<sub>ig</sub> from `composition_offset()`, with the dimnames of `Y`. |
| `cfg` | Options; see [options in `cfg`](#options-in-cfg). |

**Returns** a list with `fit` (the INLA fit), `joint` and `metadata` as for
[`observed_fit()`](#observed_fit). `metadata` also holds the dispersion scales
and `composition_diagnostics`: per gene, the standard deviation of
log m<sub>ig</sub> within regions (small values mean λ<sub>g</sub> is informed
mainly by its prior) and its mean.

## screen_summarize()

```r
screen_summarize(result)
```

`result` is the output of `observed_fit()` or `screen_joint_fit()`.

**Returns** a list with

- `pi$overall`: one row per region, compared with the spot-weighted average of
  all other regions. Columns `region`, `spots`, `score`, `evidence`,
  `evidence_higher`, `evidence_lower`, `call`, `boundary_mass`, and `p_value`
  when `pi_report_p_value = TRUE`.
- `pi$pairwise`: one row per ordered pair of regions. Columns `region`,
  `reference_region`, `score`, `evidence`, `evidence_higher`,
  `evidence_lower`, `call`, `boundary_mass` (and `p_value`).
- `pi$diagnostics`: quadrature convergence, number of nodes and refinement
  history.
- `parameters`: posterior mean and 2.5 %, 50 % and 97.5 % quantiles of
  `baseline_intercept`, `relative_mu_r2`, …, the regional means
  `absolute_mu_r1`, …, the between-gene variances `tau2_r1`, … and the noise
  variance `sigma_e2`. A joint fit of sections adds the section offsets
  `section_offset_2`, …, the section variances `nu2_r1`, … and one noise
  variance per section `sigma_e2_section_1`, …. In the adjusted mode the
  regional means are adjusted states `adjusted_state_r1`, …, and `sigma_e2` is
  replaced by the negative binomial size `nb_size0` and its inverse `alpha0`.
- `mean_contrasts`: differences of regional means, each region against the
  rest (`dmu_<r>_vs_rest`) and every pair (`dmu_<r>_vs_<s>`), with
  `p_gt_zero` and `lfsr`.
- `domains`: each region called `"higher"`, `"lower"` or `"undecided"`
  against the rest, from `p_gt_zero` ≥ 0.95 or ≤ 0.05 of its mean contrast.
- Adjusted mode only: `gene_coefficients` (gene intercepts `a_g` and
  composition slopes `lambda_g`) and the INLA marginals.
- `metadata`.

`score` is the posterior mean of π: the probability that a gene's expression
state in the region exceeds that in the reference, with the measurement error
of a gene's regional means added to the gene-to-gene variance. 0.5 means no
difference. `evidence_higher` = P(π > 0.5 | Y) and `evidence_lower` =
P(π < 0.5 | Y) are the posterior probabilities of higher and lower activity;
`evidence` equals `evidence_higher`. `call` is `"higher"` or `"lower"` when
that probability is at least 1 − `pi_alpha`/2 (0.975 by default), otherwise
`"none"`. The optional `p_value` = 2 min(P(π > 0.5), P(π < 0.5)) gives the
same calls at `p_value` ≤ `pi_alpha`; it is a posterior probability, not a
frequentist p-value, reported only for comparison with methods that report
p-values. Region indices in parameter and contrast names follow
`metadata$region_labels`.

## screen_write()

```r
screen_write(summary, out)
```

Writes the output of `screen_summarize()` (`summary`) to the folder `out`,
created if needed: `screen_summary.rds`, `overall_scores.csv`,
`pairwise_scores.csv`, `parameters.csv`, `mean_contrasts.csv`, `domains.csv`,
`runtime.csv`, `pi_diagnostics.rds`, `pi_refinement.csv` and `complete.txt`.
The π files are written only when π was computed. The adjusted mode adds
`joint_marginal_diagnostic.csv`, `gene_coefficients.csv` and
`composition_diagnostics.csv`.

## pathway_pi_integrate()

```r
pathway_pi_integrate(joint, alpha = 0.05, report_p_value = FALSE,
                     nodes = 4096, max_nodes = 65536, tolerance = 0.002)
```

| Argument | Description |
|----------|-------------|
| `joint` | The `joint` element of `observed_fit()` or `screen_joint_fit()`. |
| `alpha` | A region is called when the evidence in one direction is at least 1 − `alpha`/2. |
| `report_p_value` | Add the `p_value` column. |
| `nodes` | Initial number of quadrature nodes (a power of two, ≥ 16). |
| `max_nodes` | Maximum number of nodes (a power of two, ≥ 2 × `nodes`). |
| `tolerance` | Convergence tolerance between refinements. |

**Returns** a list with `pairwise`, `overall` and `diagnostics`, as in
`screen_summarize()`. With two regions the evidence is computed analytically;
with more it uses deterministic quasi-Monte Carlo quadrature and warns if
`max_nodes` is reached without convergence.

## pathway_pi()

```r
pathway_pi(mu_r, mu_s, tau2_r, tau2_s, kappa2 = 0)
```

| Argument | Description |
|----------|-------------|
| `mu_r`, `mu_s` | Regional means. |
| `tau2_r`, `tau2_s` | Regional between-gene variances. |
| `kappa2` | Measurement error of a gene's regional means, for example σ<sub>ε</sub><sup>2</sup>(1/n<sub>r</sub> + 1/n<sub>s</sub>) in the Gaussian model. |

**Returns** π = Φ((`mu_r` − `mu_s`) / √(`tau2_r` + `tau2_s` + `kappa2`)),
vectorized over its arguments. With `kappa2 = 0` this is
P(Z<sub>r</sub> > Z<sub>s</sub>) for independent
Z<sub>r</sub> ~ N(`mu_r`, `tau2_r`) and Z<sub>s</sub> ~ N(`mu_s`, `tau2_s`).
Fits of the Gaussian model use σ<sub>ε</sub><sup>2</sup> at each point of the
hyperparameter grid. In the negative binomial model each count carries the
information w<sub>ig</sub> = 1/(1/m<sub>ig</sub> + 1/size<sub>g</sub>), where
m<sub>ig</sub> is its expected count without regional effects and
size<sub>g</sub> = nb_size0 · S<sub>g</sub> in that configuration; gene g has
κ<sub>g</sub><sup>2</sup> = 1/W<sub>gr</sub> + 1/W<sub>gs</sub>, with
W<sub>gr</sub> the sum of w<sub>ig</sub> over the spots of region r, and
`kappa2` is the harmonic mean of κ<sub>g</sub><sup>2</sup> over genes.

## map_activity()

```r
map_activity(counts, coords, gene_sets, gene_field = TRUE,
             alpha = 0.05, out = NULL, min_genes = 5,
             threads = "2:4", verbose = TRUE)
```

Fits each gene set with [`field_fit()`](#field_fit) and calls its regions with
[`field_regions()`](#field_regions). Gene sets with fewer than `min_genes`
usable genes are skipped; a gene set whose fit fails is reported and the run
goes on. Region labels are not used.

| Argument | Description |
|----------|-------------|
| `counts` | Spots × genes matrix (dense or sparse) of normalized expression with spot row names and gene column names. |
| `coords` | Spots × 2 matrix of spot coordinates, in the spot order of `counts`. |
| `gene_sets` | Named list of gene vectors, or the path of a GMT file. |
| `gene_field` | Add a spatial field of each gene; see [`field_fit()`](#field_fit). `FALSE` is about 50 times faster but less accurate. |
| `alpha` | Joint posterior level of the regions; see [`field_regions()`](#field_regions). |
| `out` | Optional output folder. |
| `min_genes` | Minimum number of usable genes for a gene set to be fitted, at least 2. |
| `threads` | INLA `num.threads`. |
| `verbose` | Print progress messages. |

**Returns**, invisibly, a list with `gene_sets` (one row per gene set: genes in
the set, genes used, `status` `done`, `skipped` or `failed`, a message, the
posterior mean of the sd of z, the number of spots called higher and lower,
and the fit time) and `results` (for each fitted gene set, the
[`field_regions()`](#field_regions) table `spots` and the hyperparameters).
With `out`, each gene set is written to its own subfolder (`spots.csv`,
`hyperparameters.csv`) and `gene_sets.csv` is written to `out`.

## field_fit()

```r
field_fit(Y, coords, gene_field = TRUE, threads = "2:4")
```

Fits one spatial field shared by all genes of the gene set,

y<sub>ig</sub> = a<sub>g</sub> + z(s<sub>i</sub>) + w<sub>g</sub>(s<sub>i</sub>) + ε<sub>ig</sub>,
ε<sub>ig</sub> ~ N(0, σ<sub>ε</sub><sup>2</sup>),

where a<sub>g</sub> is a gene intercept, z is a Matérn field on a triangular
mesh of the spots (the map of the gene set, constrained to sum to zero) and
w<sub>g</sub> is a Matérn field of gene g on a coarser mesh, with one range
and one sd for all genes, which absorbs the spatial pattern a gene has outside
the gene set. Every gene enters z with weight 1, so genes that change in
opposite directions cancel. Both fields have penalized complexity priors
P(range < 0.1) = 0.5 and P(sd > 1) = 0.05 on coordinates scaled so that the
longer side of the section is 1, and the noise precision has the prior
P(σ<sub>ε</sub> > 2) = 0.05. The model is fitted with INLA.

| Argument | Description |
|----------|-------------|
| `Y` | Spots × genes matrix of normalized expression with spot row names, at least two genes. |
| `coords` | Spots × 2 matrix of spot coordinates, in the spot order of `Y`. |
| `gene_field` | Add the gene fields w<sub>g</sub>. `FALSE` is about 50 times faster but less accurate: spatial patterns of single genes that the gene set does not share can then appear as regions of the gene set. |
| `threads` | INLA `num.threads`. |

**Returns** a list with `inla` (the INLA fit), `node` (the mesh node of each
spot), `spots`, `sigma_e` (posterior mean of σ<sub>ε</sub> from the noise
precision), `hyperparameters` (mean and 2.5 %, 50 % and 97.5 % quantiles of
`noise_precision`, `range_z`, `sd_z` and, with gene fields, `range_w` and
`sd_w`; ranges on the scaled coordinates) and `metadata`. Pass it to
[`field_regions()`](#field_regions).

## field_regions()

```r
field_regions(fit, alpha = 0.05)
```

Calls the regions where z is above or below the baseline of the gene set. A
spot is called `higher` when it belongs to the largest set of spots on which
z > 0 everywhere with joint posterior probability at least 1 − `alpha`, and
`lower` likewise for z < 0. The sets are excursion sets (Bolin & Lindgren
2015), computed with the QC method of the excursions package. Because the
probability is joint, the region as a whole has posterior probability
1 − `alpha`, not each spot on its own.

| Argument | Description |
|----------|-------------|
| `fit` | Output of [`field_fit()`](#field_fit). |
| `alpha` | One minus the joint posterior probability of a region. |

**Returns** a data frame with one row per spot: `spot`, `z_mean`, `z_sd`,
`p_higher` = P(z > 0), `p_lower` = P(z < 0) and `call` (`higher`, `lower` or
`none`).
