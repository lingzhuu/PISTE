# Troubleshooting

**`Error: ... is not TRUE` at the start of `screen_joint_fit()`.**
An input check failed. Check that `Y`, `xy` and `region` have the same number
of spots, that `Y` has at least two genes and `region` at least two labels,
that there are no missing or infinite values, that `h > 0`, and that
`cfg$kappa_h` lies in [2.3, 3.2].

**`X must have unique, non-empty column names.`**
Give every covariate column a distinct name.

**Two constant columns in `X`.**
The intercept is added automatically. Remove the `(Intercept)` column produced
by `model.matrix()`.

**`Pi inference requires unfixed CCD hyperparameter integration.`**
Directional results need CCD integration. Do not change the internal
integration settings when `cfg$pi = TRUE`.

**`Unsupported configuration in this Gaussian-only package.` /
`gene_cs_rho has been removed` / `remove legacy design_size/non-CCD settings`.**
`cfg` contains settings from an older version of the package. Remove them.

**`Pi quadrature reached max_nodes without two stable refinements.`**
The evidence may be inaccurate for three or more regions. Recompute with more
nodes, e.g. `pathway_pi_integrate(fit$joint, max_nodes = 262144L)`, or a looser
`tolerance`.

**Errors while summarising a variance component.**
Strictly positive marginals are checked by two independent integrations. If
they disagree, the failing marginal is saved under
`logs/positive_marginal_failures/` in the working directory and an error is
raised.

**Slow fits or large memory use.**
The latent model has `mesh$n * K` spatial variables. Use a coarser mesh or more INLA threads
(`cfg$threads`).
