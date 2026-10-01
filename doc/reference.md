# Function reference

## `screen_joint_fit()`

```r
screen_joint_fit(Y, xy, region, mesh, h, family = "gaussian",
                 X = NULL, offset = NULL, cfg = list())
```

Fits the joint spatial pathway model for one pathway with INLA.

| Argument | Description |
|----------|-------------|
| `Y` | Spots x genes matrix of normalized expression, K ≥ 2. |
| `xy` | Spots x 2 matrix of coordinates. |
| `region` | Region label of each spot, at least two regions. |
| `mesh` | INLA / fmesher 2D mesh. |
| `h` | Spatial length scale; the SPDE scale is `cfg$kappa_h / h`. |
| `family` | Only `"gaussian"` is supported. |
| `X` | Optional covariate matrix with column names. |
| `offset` | Optional known additive offset of length N. |
| `cfg` | List of [options](options.md). |

Returns a list with `fit` (the INLA object, with configurations removed),
`joint` (the posterior mixture of regional means used for the directional
results, `NULL` when `cfg$pi = FALSE`) and `metadata`.

The model for spot *i* and gene *g* in region *r(i)* is

```
y[i, g] = x_i' beta          shared covariates, including the baseline intercept
        + mu_r(i)            regional offset, zero for the baseline region
        + b_r(i)[g]          gene effect, iid across genes, variance tau2_r
        + x_i' slope_g       optional gene-specific covariate slopes
        + w_g(s_i)           gene-specific Matérn field (SPDE, range fixed by h)
        + e[i, g]            Gaussian noise
```

## `screen_summarize()`

```r
screen_summarize(result)
```

Summarises the output of `screen_joint_fit()`: computes the directional
results with `pathway_pi_integrate()` and builds the parameter, contrast and
region tables. See [Output](output.md).

## `screen_write()`

```r
screen_write(summary, out)
```

Writes the output of `screen_summarize()` to the directory `out`, which is
created if needed. See [Output](output.md#files-written-by-screen_write).

## `pathway_pi_integrate()`

```r
pathway_pi_integrate(joint, nodes = 4096L, max_nodes = 65536L,
                     tolerance = 0.002)
```

Computes scores and evidence from the `joint` element of a fit. Called by
`screen_summarize()`; call it directly to recompute the directional results
with other quadrature settings without refitting:

```r
pi <- pathway_pi_integrate(fit$joint, nodes = 8192L, max_nodes = 131072L)
```

Returns a list with `pairwise`, `overall` and `diagnostics`.

## `pathway_pi()`

```r
pathway_pi(mu_r, mu_s, tau2_r, tau2_s)
```

Returns π = P(Z<sub>r</sub> > Z<sub>s</sub>) for independent
Z<sub>r</sub> ~ N(mu_r, tau2_r) and Z<sub>s</sub> ~ N(mu_s, tau2_s), for given
values of the means and variances. Vectorised over its arguments.

```r
pathway_pi(mu_r = 0.2, mu_s = 0, tau2_r = 0.1, tau2_s = 0.1)
```
