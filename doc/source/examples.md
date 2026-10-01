# Example usages

## Screen a collection of gene sets

```r
mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(100, 300), cutoff = 20)

for (set_name in names(gene_sets)) {
  genes <- intersect(gene_sets[[set_name]], colnames(expr))
  if (length(genes) < 2) next
  fit <- screen_joint_fit(expr[, genes], xy, region, mesh, h)
  screen_write(screen_summarize(fit), file.path("results", set_name))
}
```

`gene_sets` is a named list of gene vectors. The mesh is shared by all sets.

## Collect the region-versus-rest results of all sets

```r
overall <- do.call(rbind, lapply(names(gene_sets), function(set_name) {
  f <- file.path("results", set_name, "overall_scores.csv")
  if (file.exists(f)) cbind(gene_set = set_name, read.csv(f))
}))
subset(overall, region == "WM")
```

## Choose the baseline region

```r
fit <- screen_joint_fit(Y, xy, region, mesh, h,
                        cfg = list(baseline_region = "WM"))
```

## Adjust for spot-level covariates

```r
X <- model.matrix(~ batch, data = spot_meta)[, -1, drop = FALSE]
fit <- screen_joint_fit(Y, xy, region, mesh, h, X = X)
```

The intercept is added automatically, so drop the `(Intercept)` column. Column
names must be unique.

To let a covariate have a different slope for each gene:

```r
fit <- screen_joint_fit(Y, xy, region, mesh, h, X = X,
                        cfg = list(gene_specific_covariates = colnames(X)))
res <- screen_summarize(fit)
res$gene_covariate_coefficients
```

## Fit without directional inference

```r
fit <- screen_joint_fit(Y, xy, region, mesh, h, cfg = list(pi = FALSE))
```

The summary then contains the parameters and mean contrasts but no `pi`
results.

## Recompute scores and evidence with more quadrature nodes

With three or more regions the evidence is computed by quadrature. If
`screen_summarize()` warns that the quadrature did not converge, recompute it
from the same fit:

```r
pi <- pathway_pi_integrate(fit$joint, nodes = 8192L, max_nodes = 262144L)
pi$diagnostics$converged
```

## Compute π for given regional moments

```r
pathway_pi(mu_r = 0.2, mu_s = 0, tau2_r = 0.1, tau2_s = 0.1)
```

## Options

Options are passed to `screen_joint_fit()` as `cfg = list(...)`.

| Option | Default | Description |
|--------|---------|-------------|
| `pi` | `TRUE` | Compute directional scores and evidence. |
| `pi_nodes`, `pi_max_nodes`, `pi_tolerance` | `4096`, `65536`, `0.002` | Quadrature settings for three or more regions. |
| `threads` | `"2:1"` | INLA threads. |
| `strategy` | `"gaussian"` | INLA approximation strategy. |
| `kappa_h` | `sqrt(8)` | SPDE scale times `h`, in [2.3, 3.2]. |
| `tau_prior` | `"normal"` | Prior on the log precision of the regional gene effects: `"normal"`, `"inla_default"`, `"pc"` or `"flat"`. |
| `tau_logprecision_mean`, `tau_logprecision_sd` | `0`, `3` | Normal prior parameters. |
| `tau_pc_u`, `tau_pc_alpha` | `2`, `0.05` | PC prior parameters. |
| `baseline_region` | first label | Baseline region. |
| `other_weights` | spot counts | Weights of the other regions in the baseline-versus-rest mean contrast. |
| `gene_specific_covariates` | none | Columns of `X` with a slope per gene. |
| `gene_covariate_prior_sd` | `10` | Prior standard deviation of the gene-specific slopes. |
