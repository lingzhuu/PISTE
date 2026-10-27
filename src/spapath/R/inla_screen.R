# =============================================================================
# inla_screen.R
#
# Fit, summarise and export the spatial pathway model with INLA, using CCD
# hyperparameter integration and optional deterministic pi integration.
# N spots x K genes; R regions, region 1 = baseline.
#
# Observed mode (Gaussian response, normalised expression)
#
#   y[i, g] = a_g                gene intercept shared by all regions,
#                                a_g = alpha0 + u_g (see below)
#           + x_i' beta          shared covariates other than the intercept
#           + mu_r(i)            relative regional offset, mu_1 = 0
#           + gene_r(i)[g]       region-specific iid gene deviation from a_g,
#                                variance tau2_r
#           + x_i' slope_g       optional gene-specific covariate slopes
#           + spatial_g(s_i)     replicated per-gene SPDE Matern field
#                                (optional: cfg$spatial_field = TRUE adds it;
#                                off by default)
#           + e[i, g]            Gaussian noise, variance sigma_e2
#
#   The shared gene intercept a_g is essential for pi: without it a gene's
#   overall level enters both gene_r[g] and gene_s[g], inflates tau2_r and
#   tau2_s, and pulls pi towards 0.5. With a_g, gene_r[g] and gene_s[g] are
#   conditionally independent and Var(Delta_g) = tau2_r + tau2_s for the
#   paired contrast Delta_g = (mu_r + gene_r[g]) - (mu_s + gene_s[g]), which
#   is the quantity pi is defined on.
#
# Adjusted mode (negative binomial response, raw counts)
#
#   y[i, g] ~ NB(mean L_i exp(eta[i, g]), size nb_size0 * S_g)
#   eta[i, g] = a_g + x_i' beta + mu_r(i) + gene_r(i)[g] + spatial_g(s_i)
#             + lambda_g * log m[i, g]
#
#   L_i      library size; S_g relative dispersion scale of gene g
#   m[i, g]  expression predicted from the spot's cell-type composition
#   a_g      gene intercept, a_g = alpha0 + u_g with alpha0 ~ N(0, 100 (1 - 1/K))
#            and u_g ~ N(0, 100); x_i' beta holds the other covariates.
#            The same a_g is used in both modes.
#   lambda_g composition slope, lambda_g = 1 + v_g with v_g ~ N(0, 0.25^2)
#
# In both modes the baseline state is alpha0 + mean(u), and
# alpha0 + mean(u) + mu_r is the pathway state of region r (adjusted for
# cell-type composition in the adjusted mode).
#
# Free hyperparameters (R + 2): the R gene-effect precisions, one SPDE
# variance parameter (the range is fixed through h) and the Gaussian
# observation precision (observed) or the NB size nb_size0 (adjusted).
#
# Typical workflow
#   mesh    <- build_mesh(xy)
#   fit     <- screen_joint_fit(Y, xy, region, mesh, h)
#   summary <- screen_summarize(fit)
#   screen_write(summary, "results/")
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

#' Default options for screen_joint_fit()
#' @noRd
screen_default_cfg <- function(mode) {
  list(
    # Directional probability (pi) inference
    pi = TRUE,
    pi_delta = .05,
    pi_level = .9,
    pi_nodes = 4096,
    pi_max_nodes = 65536,
    pi_tolerance = .002,

    # Gene-specific spatial fields w_g; off by default, TRUE adds them
    spatial_field = FALSE,

    # INLA threads; SPDE range kappa = kappa_h / h; INLA approximation strategy
    threads = "1:1",
    kappa_h = sqrt(8),
    strategy = if (mode == "observed") "gaussian" else "simplified.laplace",

    # Normal prior on the log NB size (adjusted mode)
    nb_logsize_mean = log(10),
    nb_logsize_sd = 1.5,

    # INLA run controls: AMD reordering with a small node factor to lower
    # memory on large gene sets (same model, permuted precision matrix),
    # starting hyperparameters for a restart, and INLA's verbose log
    low_memory = FALSE,
    initial_theta = NULL,
    verbose = FALSE,

    # Prior on the log precision of the regional gene effects
    tau_prior = "normal",
    tau_logprecision_mean = 0,
    tau_logprecision_sd = 3,
    tau_pc_u = 2,
    tau_pc_alpha = .05,

    # Weights of the other regions in the "region 1 vs rest" contrast
    other_weights = NULL,

    # Columns of X with independent per-gene slopes
    gene_specific_covariates = character(),
    gene_covariate_prior_sd = 10,

    # Baseline region label (default: first label found)
    baseline_region = NULL
  )
}


#' Reject settings from older, incompatible versions of the package
#' @noRd
screen_check_legacy_cfg <- function(cfg) {
  if ("gene_cs_rho" %in% names(cfg)) {
    stop("gene_cs_rho has been removed; region-specific gene effects are conditionally iid. Remove this legacy setting.")
  }
}


#' Reject settings that do not apply to the chosen mode
#' @noRd
screen_check_mode_cfg <- function(cfg, mode) {
  if (mode == "observed" && any(grepl("^nb_", names(cfg)))) {
    stop("nb_* options apply to the adjusted mode only.")
  }
  if (mode == "adjusted" && length(cfg$gene_specific_covariates)) {
    stop("The adjusted mode does not support gene_specific_covariates.")
  }
}


#' @noRd
screen_check_integration_cfg <- function(cfg) {
  if (!is.null(cfg$design_size) ||
      (!is.null(cfg$integration) && cfg$integration != "ccd")) {
    stop("Production screening uses CCD; remove legacy design_size/non-CCD settings.")
  }
}


# -----------------------------------------------------------------------------
# Model building blocks
# -----------------------------------------------------------------------------

#' Region labels with the baseline region first
#' @noRd
screen_region_labels <- function(region, baseline_region) {
  found_labels <- unique(as.character(region))
  baseline <- if (is.null(baseline_region)) found_labels[1] else
    as.character(baseline_region)
  stopifnot(length(baseline) == 1, baseline %in% found_labels)
  c(baseline, setdiff(found_labels, baseline))
}


#' Prepare the covariate matrix and split it into shared and gene-specific parts
#'
#' Adds an `intercept` column when missing; it represents the baseline
#' regional mean and must be the only constant-one column.
#'
#' @return List with the completed matrix `X`, INLA effect names `names`
#'   (`cov_<column>`), column indices `shared` and `gene_specific`, and the
#'   effect / weight names of the gene-specific slopes.
#' @noRd
screen_prepare_covariates <- function(X, N, cfg) {
  if (is.null(X)) X <- matrix(numeric(), N, 0)
  X <- as.matrix(X)
  stopifnot(nrow(X) == N, all(is.finite(X)))

  if (!"intercept" %in% colnames(X)) {
    X <- cbind(intercept = rep(1, N), X)
  }
  stopifnot(all(X[, "intercept"] == 1),
            sum(colSums(X == 1) == N) == 1)

  if (anyDuplicated(colnames(X)) || any(!nzchar(colnames(X)))) {
    stop("X must have unique, non-empty column names.")
  }
  effect_names <- paste0("cov_", make.names(colnames(X), unique = TRUE))

  # gene_x: columns of X with gene-specific slopes.
  gene_x <- match(cfg$gene_specific_covariates, colnames(X))
  stopifnot(!anyNA(gene_x),
            !anyDuplicated(gene_x),
            !(match("intercept", colnames(X)) %in% gene_x),
            length(cfg$gene_covariate_prior_sd) == 1,
            is.finite(cfg$gene_covariate_prior_sd),
            cfg$gene_covariate_prior_sd > 0)

  slope_names <- if (length(gene_x)) {
    paste0("slope_", effect_names[gene_x])
  } else {
    character()
  }

  list(
    X = X,
    names = effect_names,
    shared = setdiff(seq_len(ncol(X)), gene_x),
    gene_specific = gene_x,
    slope_names = slope_names,
    slope_weight_names = if (length(gene_x)) paste0(slope_names, "_weight") else character()
  )
}


#' SPDE Matern model with the range fixed by kappa
#'
#' Only the marginal variance of the field is estimated
#' (one hyperparameter, prior N(0, 1)).
#' @noRd
screen_fixed_range_spde <- function(mesh, kappa) {
  spde <- inla.spde2.matern(
    mesh, alpha = 2, constr = FALSE,
    B.tau = matrix(c(-log(sqrt(4 * pi) * kappa), -1), 1, 2),
    B.kappa = matrix(c(log(kappa), 0), 1, 2),
    theta.prior.mean = 0, theta.prior.prec = 1
  )
  stopifnot(spde$n.theta == 1)
  spde
}


#' Effects table for the long-format (spot, gene) observations
#'
#' One column per shared covariate, per gene-specific slope (index + weight),
#' per regional offset indicator and per regional gene effect. The regional
#' gene index `gene_r<r>` is NA outside region r.
#' @noRd
screen_effects_table <- function(covariates, site, gene, region_long,
                                 offset_names, gene_names) {
  X <- covariates$X
  n_regions <- length(gene_names)

  effects <- data.frame(row = seq_along(site))
  for (j in covariates$shared) {
    effects[[covariates$names[j]]] <- X[site, j]
  }
  for (j in seq_along(covariates$gene_specific)) {
    effects[[covariates$slope_names[j]]] <- gene
    effects[[covariates$slope_weight_names[j]]] <- X[site, covariates$gene_specific[j]]
  }
  for (r in seq_len(n_regions)) {
    if (r > 1) effects[[offset_names[r - 1]]] <- as.numeric(region_long == r)
    effects[[gene_names[r]]] <- ifelse(region_long == r, gene, NA)
  }
  effects$row <- NULL

  # Every observation belongs to exactly one regional gene effect.
  stopifnot(all(rowSums(!is.na(effects[gene_names])) == 1))
  effects
}


#' INLA hyperprior for the regional gene-effect precisions
#'
#' INLA's internal iid hyperparameter is the log precision; for the "normal"
#' prior `param[2]` is the Normal *precision*. The default normal prior is a
#' proper approximation to a log-flat scale prior.
#' @noRd
screen_tau_hyper <- function(cfg) {
  switch(cfg$tau_prior,
    normal = list(prec = list(
      prior = "normal",
      param = c(cfg$tau_logprecision_mean, 1 / cfg$tau_logprecision_sd^2)
    )),
    inla_default = list(prec = list(prior = "loggamma", param = c(1, 5e-5))),
    pc = list(prec = list(prior = "pc.prec",
                          param = c(cfg$tau_pc_u, cfg$tau_pc_alpha))),
    flat = list(prec = list(prior = "flat"))
  )
}


#' Right-hand side of the model formula, as text
#'
#' The `f()` terms refer to `iid_hyper`, `slope_hyper`, `gene_intercept_hyper`,
#' `composition_hyper`, `spde` and `K`, which INLA looks up in the environment
#' of the formula, so the formula itself must be created in
#' `screen_joint_fit()`.
#' @noRd
screen_formula_text <- function(covariates, offset_names, gene_names,
                                adjusted, spatial_field = TRUE) {
  gene_terms <- paste0('f(', gene_names,
                       ',model="iid",n=K,constr=FALSE,hyper=iid_hyper)')
  slope_terms <- if (length(covariates$gene_specific)) {
    paste0('f(', covariates$slope_names, ',', covariates$slope_weight_names,
           ',model="iid",n=K,constr=FALSE,hyper=slope_hyper)')
  } else {
    character()
  }
  # Both modes: auxiliary gene intercepts u_g (a_g = alpha0 + u_g), shared by
  # all regions, so that gene_r[g] are deviations from the gene's own level.
  gene_intercept_term <-
    'f(gene_intercept,model="iid",n=K,constr=FALSE,hyper=gene_intercept_hyper)'
  # Adjusted mode: composition slopes v_g multiplying log m[i, g].
  composition_terms <- if (adjusted) {
    'f(composition_slope,composition_weight,model="iid",n=K,constr=FALSE,hyper=composition_hyper)'
  } else {
    character()
  }
  spatial_term <- if (spatial_field) {
    'f(spatial,model=spde,replicate=spatial.repl,nrep=K,constr=FALSE)'
  } else {
    character()
  }

  paste("y ~ -1 +", paste(c(covariates$names[covariates$shared], offset_names,
                            gene_terms, slope_terms, gene_intercept_term,
                            composition_terms, spatial_term),
                          collapse = " + "))
}


#' Linear combinations of the regional means
#'
#' Builds three groups of INLA linear combinations:
#' * `dmu_<r>_vs_rest`: region r minus a weighted average of all other regions
#'   (weights proportional to spot counts, or `other_weights` for region 1);
#' * `dmu_<r>_vs_<s>`: every pairwise difference;
#' * `absolute_mu_r<r>`: absolute regional means, reported as diagnostics
#'   while mu_r2, ..., mu_rR stay the primary relative coefficients.
#'
#' Contrasts are expressed on the relative offsets only (region 1 has offset
#' zero, so its coefficient is dropped). In adjusted mode (`n_genes` > 0) the
#' absolute means also include mean(u) of the auxiliary gene intercepts, and
#' `gene_intercept_g<g>` = intercept + u_g gives each gene intercept a_g.
#'
#' @param n_genes Number of auxiliary gene intercepts (0 in observed mode).
#' @return List with `lincombs`, the full length-R `coefficients` of each
#'   contrast and the `rest_weights` used for each region.
#' @noRd
screen_region_lincombs <- function(counts, other_weights, intercept_name,
                                   offset_names, n_genes = 0) {
  R <- length(counts)
  lincombs <- list()
  coefficients <- list()
  rest_weights <- vector("list", R)

  # INLA linear combination of the relative offsets for coefficients z.
  offset_lincomb <- function(z, name) {
    setNames(do.call(inla.make.lincomb,
                     as.list(setNames(z[-1], offset_names))), name)
  }

  # Region r versus the rest.
  for (r in seq_len(R)) {
    rest <- setdiff(seq_len(R), r)
    a <- counts[rest] / sum(counts[rest])
    if (r == 1 && !is.null(other_weights)) {
      a <- other_weights
      stopifnot(length(a) == R - 1, all(is.finite(a)), all(a >= 0),
                sum(a) > 0)
      a <- a / sum(a)
    }
    z <- numeric(R)
    z[r] <- 1
    z[rest] <- -a
    rest_weights[[r]] <- a
    name <- paste0("dmu_", r, "_vs_rest")
    lincombs <- c(lincombs, offset_lincomb(z, name))
    coefficients[[name]] <- z
  }

  # All pairwise differences.
  for (pair in combn(seq_len(R), 2, simplify = FALSE)) {
    z <- numeric(R)
    z[pair] <- c(1, -1)
    name <- paste0("dmu_", pair[1], "_vs_", pair[2])
    lincombs <- c(lincombs, offset_lincomb(z, name))
    coefficients[[name]] <- z
  }

  # Absolute regional means: intercept (+ mean(u)) + offset_r.
  for (r in seq_len(R)) {
    z <- c(1, as.numeric(seq.int(2, R) == r))
    names(z) <- c(intercept_name, offset_names)
    z <- as.list(z)
    if (n_genes > 0) z$gene_intercept <- rep(1 / n_genes, n_genes)
    lincombs <- c(lincombs, setNames(do.call(inla.make.lincomb, z),
                                     paste0("absolute_mu_r", r)))
  }

  # Gene intercepts a_g = intercept + u_g.
  for (g in seq_len(n_genes)) {
    z <- setNames(list(1, as.numeric(seq_len(n_genes) == g)),
                  c(intercept_name, "gene_intercept"))
    lincombs <- c(lincombs, setNames(do.call(inla.make.lincomb, z),
                                     paste0("gene_intercept_g", g)))
  }

  list(lincombs = lincombs, coefficients = coefficients,
       rest_weights = rest_weights)
}


# -----------------------------------------------------------------------------
# Adjusted-mode inputs
# -----------------------------------------------------------------------------

#' Validate the dispersion and composition inputs
#' @noRd
screen_check_adjusted_inputs <- function(Y, library_size, dispersion_scale,
                                         composition_offset) {
  if (length(dispersion_scale) != ncol(Y) || any(!is.finite(dispersion_scale)) ||
      any(dispersion_scale <= 0)) {
    stop("dispersion_scale must be positive and finite, one value per gene.")
  }
  if (!is.null(names(dispersion_scale)) &&
      !identical(names(dispersion_scale), colnames(Y))) {
    stop("names(dispersion_scale) must match colnames(Y) in order.")
  }
  if (!is.matrix(composition_offset) || !is.numeric(composition_offset) ||
      !identical(dimnames(composition_offset), dimnames(Y)) ||
      any(!is.finite(composition_offset))) {
    stop("composition_offset must be a finite matrix with the dimnames of Y.")
  }
}


#' Expected counts without regional effects, for the NB measurement error
#'
#' m[i, g] = c_g L_i exp(offset_i) m_composition[i, g], where the gene scale
#' c_g matches the expected to the observed total count of gene g. The
#' composition slope is taken at its prior mean lambda_g = 1.
#' @noRd
screen_expected_counts <- function(Y, library_size, offset, composition_offset,
                                   dispersion_scale, region_index) {
  expected <- library_size * exp(offset + composition_offset)
  gene_scale <- colSums(Y) / colSums(expected)
  list(expected = sweep(expected, 2, gene_scale, "*"),
       dispersion_scale = as.numeric(dispersion_scale),
       region_index = region_index)
}


#' Per-gene diagnostics of the composition offset
#'
#' `within_region_sd` is the standard deviation of log m[i, g] after removing
#' the covariates and region means. Small values mean lambda_g is identified
#' mainly by its prior.
#' @noRd
screen_composition_diagnostics <- function(composition_offset, X, region_index) {
  region_indicators <- outer(region_index, seq.int(2, max(region_index)), "==") * 1
  design <- qr(cbind(X, region_indicators))
  if (design$rank < ncol(design$qr)) {
    stop("X is rank deficient with the intercept/region indicators; remove duplicate region or constant columns.")
  }
  data.frame(gene = colnames(composition_offset),
             within_region_sd = apply(qr.resid(design, composition_offset), 2,
                                      stats::sd),
             mean_log_offset = colMeans(composition_offset),
             row.names = NULL)
}


# -----------------------------------------------------------------------------
# Model fitting
# -----------------------------------------------------------------------------

#' Fit the joint spatial pathway model with INLA
#'
#' @param Y N x K matrix (spots x genes), K >= 2. Normalised expression in
#'   observed mode; raw integer counts with unique spot and gene names in
#'   adjusted mode.
#' @param xy N x 2 matrix of spot coordinates.
#' @param region Length-N vector of region labels (at least two regions).
#' @param mesh Mesh for the SPDE fields, from [build_mesh()]. Needed only
#'   when `cfg$spatial_field` is `TRUE`; otherwise it may be `NULL`.
#' @param h Spatial length scale; the SPDE parameter kappa is fixed at
#'   `cfg$kappa_h / h`. Needed only when `cfg$spatial_field` is `TRUE`;
#'   otherwise it may be `NULL`.
#' @param mode `"observed"` (Gaussian model of normalised expression) or
#'   `"adjusted"` (negative binomial model of raw counts, adjusted for
#'   cell-type composition).
#' @param X Optional N x p covariate matrix with unique column names. An
#'   `intercept` column is added if absent; it represents the baseline
#'   regional mean.
#' @param offset Optional known additive offset of length N.
#' @param library_size Adjusted mode: length-N library sizes, e.g. total UMI
#'   counts per spot. Only their ratios matter; they are rescaled to a
#'   geometric mean of one.
#' @param dispersion_scale Adjusted mode: length-K relative dispersion scales
#'   S_g from [estimate_dispersion_scale()]; gene g has NB size
#'   `nb_size0 * S_g`.
#' @param composition_offset Adjusted mode: N x K matrix of log predicted
#'   expression log m[i, g] from [composition_offset()], with the dimnames of
#'   `Y`.
#' @param cfg List of options overriding the defaults:
#'   * `pi` (TRUE), `pi_delta` (0.05), `pi_level` (0.95), `pi_nodes`,
#'     `pi_max_nodes`, `pi_tolerance`: directional probability inference,
#'     see [pathway_pi_integrate()].
#'   * `threads` ("1:1"), `strategy` ("gaussian" in observed mode,
#'     "simplified.laplace" in adjusted mode): INLA settings.
#'   * `spatial_field` (FALSE): add the gene-specific spatial fields
#'     spatial_g. By default the model is fitted without them, so the
#'     regional gene effects alone carry each gene's departure from the
#'     regional means; TRUE adds the fields (one more hyperparameter, slower).
#'   * `kappa_h` (`sqrt(8)`): kappa * h, must lie in \[2.3, 3.2\].
#'   * `tau_prior`: prior on the gene-effect log precisions, one of
#'     `"normal"` (default; `tau_logprecision_mean`, `tau_logprecision_sd`),
#'     `"inla_default"`, `"pc"` (`tau_pc_u`, `tau_pc_alpha`) or `"flat"`.
#'   * `nb_logsize_mean` (`log(10)`), `nb_logsize_sd` (1.5): Normal prior on
#'     log nb_size0, adjusted mode only.
#'   * `other_weights`: weights of the non-baseline regions in the
#'     "region 1 vs rest" contrast (default: proportional to spot counts).
#'   * `gene_specific_covariates`, `gene_covariate_prior_sd` (10): columns of
#'     `X` that get independent per-gene slopes with a fixed Normal prior
#'     (observed mode only).
#'   * `baseline_region`: label of the baseline region (default: the first
#'     label found in `region`).
#'   * `low_memory` (FALSE): AMD reordering with a small global node factor;
#'     the same model with lower memory use for large gene sets.
#'   * `initial_theta` (NULL): starting hyperparameters (length R + 1, or
#'     R + 2 with the spatial fields), e.g. `fit$fit$mode$theta` of an
#'     earlier run.
#'   * `verbose` (FALSE): print INLA's log.
#' @param .fixed_theta Internal structural test hook: fix the hyperparameters.
#' @param .int_strategy INLA integration strategy, `"ccd"` or `"eb"`.
#' @return A list with the INLA `fit` (configurations removed), the extracted
#'   `joint` mixture used by [pathway_pi_integrate()] (NULL when `cfg$pi` is
#'   FALSE) and a `metadata` list.
#' @export
screen_joint_fit <- function(Y, xy, region, mesh, h,
                             mode = c("observed", "adjusted"),
                             X = NULL, offset = NULL,
                             library_size = NULL, dispersion_scale = NULL,
                             composition_offset = NULL,
                             cfg = list(), .fixed_theta = NULL,
                             .int_strategy = "ccd") {
  # --- Configuration ---------------------------------------------------------
  mode <- match.arg(mode)
  adjusted <- mode == "adjusted"
  family <- if (adjusted) "nbinomial" else "gaussian"

  screen_check_legacy_cfg(cfg)
  screen_check_mode_cfg(cfg, mode)
  stopifnot(.int_strategy %in% c("ccd", "eb"))
  screen_check_integration_cfg(cfg)

  cfg <- modifyList(screen_default_cfg(mode), cfg)
  if (isTRUE(cfg$pi) && (.int_strategy != "ccd" || !is.null(.fixed_theta))) {
    stop("Pi inference requires unfixed CCD hyperparameter integration.")
  }
  stopifnot(is.logical(cfg$spatial_field), length(cfg$spatial_field) == 1,
            !is.na(cfg$spatial_field))
  spatial_field <- cfg$spatial_field

  # --- Dimensions and regions ------------------------------------------------
  N <- nrow(Y)
  K <- ncol(Y)
  labels <- screen_region_labels(region, cfg$baseline_region)
  region_index <- match(as.character(region), labels)
  R <- length(labels)
  # Hyperparameters: R regional log precisions, the observation precision or
  # NB size, and the spatial log standard deviation when fields are fitted.
  n_hyper <- R + 1 + spatial_field

  if (spatial_field) {
    stopifnot(!is.null(mesh), is.numeric(h), length(h) == 1,
              is.finite(h), h > 0,
              cfg$kappa_h >= 2.3, cfg$kappa_h <= 3.2)
  }
  stopifnot(N == nrow(xy), length(region) == N,
            R >= 2, K >= 2,
            all(is.finite(Y)), !anyNA(region),
            is.finite(cfg$tau_logprecision_mean),
            is.finite(cfg$tau_logprecision_sd),
            cfg$tau_logprecision_sd > 0,
            cfg$tau_prior %in% c("normal", "inla_default", "pc", "flat"),
            is.null(cfg$initial_theta) ||
              (length(cfg$initial_theta) == n_hyper &&
                 all(is.finite(cfg$initial_theta))))

  covariates <- screen_prepare_covariates(X, N, cfg)
  X <- covariates$X

  if (is.null(offset)) offset <- rep(0, N)
  stopifnot(length(offset) == N, all(is.finite(offset)))

  composition <- NULL
  if (adjusted) {
    screen_check_adjusted_inputs(Y, library_size, dispersion_scale,
                                 composition_offset)
    stopifnot(cfg$nb_logsize_sd > 0)
    library_size <- library_size / exp(mean(log(library_size)))
    composition <- screen_composition_diagnostics(composition_offset, X,
                                                  region_index)
  }

  # Long format: one observation per (spot, gene); one field replicate per gene.
  site <- rep(seq_len(N), K)
  gene <- rep(seq_len(K), each = N)
  region_long <- region_index[site]

  # --- Spatial fields (optional) ---------------------------------------------
  kappa <- NA_real_
  spde <- NULL
  if (spatial_field) {
    kappa <- cfg$kappa_h / h
    spde <- screen_fixed_range_spde(mesh, kappa)

    spatial_index <- inla.spde.make.index("spatial", mesh$n, n.repl = K)
    A <- inla.spde.make.A(mesh, loc = xy, index = site, repl = gene, n.repl = K)

    # Each row of A must only touch its own gene's replicate and sum to one.
    nz <- summary(A)  # non-zero entries (i, j, x)
    stopifnot(all((nz$j - 1) %/% mesh$n + 1 == gene[nz$i]),
              max(abs(rowSums(A) - 1)) < 1e-9)
  }

  # --- Effects and data stack ------------------------------------------------
  # The intercept in X is the baseline regional mean. Other mu coefficients
  # are offsets from that baseline, with delta_1 fixed at zero.
  offset_names <- paste0("mu_r", seq.int(2, R))
  intercept_name <- covariates$names[match("intercept", colnames(X))]
  gene_names <- paste0("gene_r", seq_len(R))

  effects <- screen_effects_table(covariates, site, gene, region_long,
                                  offset_names, gene_names)
  if (adjusted) {
    # Composition slope lambda_g = 1 + v_g: the fixed part enters through the
    # INLA offset below, v_g multiplies log m[i, g] here.
    effects$composition_slope <- gene
    effects$composition_weight <- as.vector(composition_offset)
  }

  # Auxiliary gene intercepts u_g, one per observation's gene (both modes).
  gene_map <- Matrix::sparseMatrix(i = seq_len(N * K), j = gene, x = 1,
                                   dims = c(N * K, K))
  stack_A <- list(1, gene_map)
  stack_effects <- list(effects, list(gene_intercept = seq_len(K)))
  if (spatial_field) {
    stack_A <- c(stack_A, list(A))
    stack_effects <- c(stack_effects, list(spatial_index))
  }
  stack <- inla.stack(list(y = as.vector(Y)), A = stack_A,
                      effects = stack_effects,
                      compress = FALSE, remove.unused = FALSE)

  # --- Priors and formula ----------------------------------------------------
  # `iid_hyper`, `slope_hyper`, `gene_intercept_hyper`, `composition_hyper`,
  # `spde` and `K` are referenced by name inside the formula's f() terms.
  iid_hyper <- screen_tau_hyper(cfg)

  # Independent per-gene slopes with a fixed prior precision: no extra
  # hyperparameters and no dense gene-by-covariate dummy matrix.
  slope_hyper <- list(prec = list(initial = log(1 / cfg$gene_covariate_prior_sd^2),
                                  fixed = TRUE))

  # u_g ~ N(0, 100) and v_g ~ N(0, 0.25^2), both with fixed precision.
  gene_intercept_hyper <- list(prec = list(initial = log(.01), fixed = TRUE))
  composition_hyper <- list(prec = list(initial = log(1 / .25^2), fixed = TRUE))

  formula <- as.formula(screen_formula_text(covariates, offset_names,
                                            gene_names, adjusted,
                                            spatial_field))

  # --- Regional contrasts ----------------------------------------------------
  counts <- as.numeric(table(factor(region_index, levels = seq_len(R))))
  contrasts <- screen_region_lincombs(counts, cfg$other_weights,
                                      intercept_name, offset_names,
                                      n_genes = K)

  # --- Run INLA --------------------------------------------------------------
  observation_family <- if (adjusted) {
    list(variant = 2, hyper = list(theta = list(
      prior = "normal", param = c(cfg$nb_logsize_mean, 1 / cfg$nb_logsize_sd^2))))
  } else {
    list(hyper = list(prec = list(prior = "pc.prec", param = c(2, .05))))
  }

  inla_args <- list(
    formula = formula,
    family = family,
    data = inla.stack.data(stack),
    lincomb = contrasts$lincombs,
    control.family = observation_family,
    control.fixed = list(mean = 0, prec = .01),
    control.predictor = list(A = inla.stack.A(stack), compute = FALSE),
    control.compute = list(config = isTRUE(cfg$pi), dic = FALSE,
                           waic = FALSE, cpo = FALSE),
    control.inla = list(strategy = cfg$strategy, int.strategy = .int_strategy),
    inla.mode = "experimental",
    num.threads = cfg$threads,
    verbose = isTRUE(cfg$verbose)
  )
  # alpha0 ~ N(0, 100 (1 - 1/K)), so that alpha0 + mean(u) ~ N(0, 100).
  inla_args$control.fixed$prec <- c(list(default = .01),
                                    setNames(list(1 / (100 * (1 - 1 / K))),
                                             intercept_name))
  if (adjusted) {
    inla_args$E <- library_size[site]
    inla_args$scale <- as.numeric(dispersion_scale[gene])
    inla_args$offset <- offset[site] + as.vector(composition_offset)
  } else if (any(offset != 0)) {
    inla_args$offset <- offset[site]
  }
  if (isTRUE(cfg$low_memory)) {
    # A permutation of the same precision matrix; no approximation or prior change.
    inla_args$control.inla$reordering <- "amd"
    inla_args$control.inla$global.node.factor <- .01
  }
  if (!is.null(.fixed_theta)) {
    # Structural test hook only: conditions known hyperparameters, no sampling.
    inla_args$control.mode <- list(theta = .fixed_theta, fixed = TRUE,
                                   restart = FALSE)
  } else if (!is.null(cfg$initial_theta)) {
    inla_args$control.mode <- list(theta = cfg$initial_theta, restart = TRUE,
                                   fixed = FALSE)
  }

  # The data and predictor matrix are now in inla_args; free the intermediates
  # before INLA starts.
  rm(stack, stack_A, stack_effects, effects, gene_map)
  if (spatial_field) rm(A, nz, spatial_index)
  gc(FALSE)

  started <- proc.time()[["elapsed"]]
  fit <- do.call(inla, inla_args)
  fit_seconds <- proc.time()[["elapsed"]] - started
  gc(FALSE)

  # --- Extract the low-dimensional joint posterior, then drop configs --------
  gene_intercept_name <- "gene_intercept"
  joint <- NULL
  config_bytes <- 0
  if (isTRUE(cfg$pi)) {
    config_bytes <- as.numeric(object.size(fit$misc$configs))
    nb_counts <- if (adjusted) {
      screen_expected_counts(Y, library_size, offset, composition_offset,
                             dispersion_scale, region_index)
    }
    joint <- pi_extract_joint(fit, labels, counts, intercept_name, offset_names,
                              gene_intercept_name, nb_counts)
    rm(nb_counts)
    fit$misc$configs <- NULL
  }

  stopifnot(is.null(fit$misc$configs),
            all(vapply(fit$summary.random[gene_names], nrow,
                       integer(1)) == K))
  if (spatial_field) {
    stopifnot(nrow(fit$summary.random$spatial) == mesh$n * K)
  } else {
    stopifnot(is.null(fit$summary.random$spatial))
  }
  if (is.null(.fixed_theta)) {
    stopifnot(length(fit$mode$theta) == n_hyper,
              nrow(fit$summary.hyperpar) == n_hyper)
  }

  # --- Result ----------------------------------------------------------------
  is_normal_prior <- cfg$tau_prior == "normal"
  offset_convention <- if (adjusted) {
    "log library size via INLA exposure E plus lambda_g * log m[i, g]"
  } else if (any(offset != 0)) {
    "known Gaussian additive offset"
  } else {
    "zero (SCT/normalized Gaussian response)"
  }

  metadata <- list(
    N = N, K = K, R = R,
    mode = mode,
    region_labels = labels,
    region_spots = counts,
    other_weights = contrasts$rest_weights[[1]],
    rest_weights = contrasts$rest_weights,
    spatial_field = spatial_field,
    h = if (spatial_field) h else NA_real_,
    kappa0 = kappa,
    rho0 = if (spatial_field) sqrt(8) / kappa else NA_real_,
    gene_dependence = "conditional_iid",
    model_version = as.character(utils::packageVersion("PathwayLGM")),
    cfg = cfg,
    family = family,
    fit_seconds = fit_seconds,
    hyperparameters = n_hyper,
    tau_prior = list(
      internal_parameter = "log precision",
      prior = cfg$tau_prior,
      mean = if (is_normal_prior) cfg$tau_logprecision_mean else NA,
      sd = if (is_normal_prior) cfg$tau_logprecision_sd else NA,
      parameters = iid_hyper$prec$param
    ),
    nb_prior = if (adjusted) observation_family$hyper else NULL,
    dispersion_scale = if (adjusted) dispersion_scale else NULL,
    composition_diagnostics = composition,
    composition_slope_prior = if (adjusted) list(mean = 1, sd = .25) else NULL,
    baseline_region = labels[1],
    intercept_name = intercept_name,
    gene_intercept_name = gene_intercept_name,
    regional_offset_names = offset_names,
    covariate_names = covariates$names[covariates$shared],
    gene_covariate_names = covariates$slope_names,
    gene_covariate_columns = colnames(X)[covariates$gene_specific],
    gene_labels = colnames(Y),
    offset_convention = offset_convention,
    integration = .int_strategy,
    config = FALSE,
    fit_config = isTRUE(cfg$pi),
    pi_enabled = isTRUE(cfg$pi),
    config_bytes_before_extraction = config_bytes,
    joint_bytes = as.numeric(object.size(joint)),
    posterior_sampling = FALSE,
    contrast_coefficients = contrasts$coefficients,
    INLA_version = as.character(packageVersion("INLA")),
    formula = deparse(formula)
  )

  list(fit = fit, joint = joint, metadata = metadata)
}




# -----------------------------------------------------------------------------
# Posterior summaries
# -----------------------------------------------------------------------------

#' Mean and 95% interval of `transform(x)` for an INLA marginal
#'
#' Set `decreasing = TRUE` for a decreasing `transform` (e.g. log precision
#' -> variance) so the interval ends are swapped.
#'
#' @return Named vector `c(mean, lower, median, upper)`.
#' @noRd
screen_marginal_summary <- function(m, transform = identity, decreasing = FALSE) {
  quantiles <- inla.qmarginal(c(.025, .5, .975), m)
  if (decreasing) quantiles <- rev(quantiles)
  c(mean = inla.emarginal(transform, m),
    lower = transform(quantiles[1]),
    median = transform(quantiles[2]),
    upper = transform(quantiles[3]))
}


#' Convert an INLA summary table to (parameter, mean, lower, median, upper)
#' @noRd
screen_summary_frame <- function(summary_table) {
  data.frame(
    parameter = rownames(summary_table),
    mean = summary_table$mean,
    lower = summary_table$`0.025quant`,
    median = summary_table$`0.5quant`,
    upper = summary_table$`0.975quant`,
    row.names = NULL
  )
}


#' Fixed effects: baseline intercept, relative offsets, covariates, then the
#' absolute regional means
#'
#' In adjusted mode the baseline intercept is alpha0 + mean(u), taken from
#' the `absolute_mu_r1` linear combination, and the absolute means are named
#' `adjusted_state_r<r>`: composition-adjusted states, not expression means.
#' @noRd
screen_fixed_effect_table <- function(fit, metadata) {
  intercept_name <- metadata$intercept_name
  offset_names <- metadata$regional_offset_names
  adjusted <- metadata$mode == "adjusted"

  fixed <- fit$summary.fixed[c(intercept_name, offset_names,
                               setdiff(metadata$covariate_names, intercept_name)), ,
                             drop = FALSE]
  table <- screen_summary_frame(fixed)
  absolute <- screen_summary_frame(fit$summary.lincomb.derived[
    paste0("absolute_mu_r", seq_len(metadata$R)), , drop = FALSE])

  # The baseline intercept is alpha0 + mean(u) in both modes.
  is_intercept <- table$parameter == intercept_name
  table[is_intercept, -1] <- absolute[1, -1]

  is_offset <- table$parameter %in% offset_names
  table$parameter[is_intercept] <- "baseline_intercept"
  table$parameter[is_offset] <- paste0("relative_", table$parameter[is_offset])

  if (adjusted) {
    absolute$parameter <- sub("^absolute_mu_", "adjusted_state_", absolute$parameter)
  }
  rbind(table, absolute)
}


#' Internal-scale marginal of a hyperparameter
#'
#' INLA integrates over the hyperparameters on their internal scale (log
#' precision, log size, Theta1), where the marginals are much less skewed
#' and well covered by INLA's grid. On the user scale a precision marginal
#' can be extremely right-skewed, so its quantiles are computed here from the
#' internal marginal and transformed afterwards. `internal.marginals.hyperpar`
#' lists the hyperparameters in the same order as `marginals.hyperpar`, so
#' the internal marginal is looked up by the user-scale name.
#' @noRd
screen_internal_marginal <- function(fit, name) {
  internal <- fit$internal.marginals.hyperpar
  j <- match(name, names(fit$marginals.hyperpar))
  if (is.na(j) || length(internal) != length(fit$marginals.hyperpar)) {
    stop("INLA did not return the internal marginal of '", name, "'.")
  }
  internal[[j]]
}


#' Variance components on the variance scale
#'
#' Summarised from the internal-scale marginals (see
#' `screen_internal_marginal()`): tau2_r = exp(-log precision of gene_r) and
#' sigma_s2 = exp(2 * Theta1) for the spatial fields (when fitted); then
#' sigma_e2 = exp(-log Gaussian observation precision) (observed mode) or the
#' NB size nb_size0 = exp(log size) and its inverse alpha0 (adjusted mode).
#' @noRd
screen_variance_table <- function(fit, n_regions, family, spatial_field = TRUE) {
  hyper_names <- paste0("Precision for gene_r", seq_len(n_regions))
  names(hyper_names) <- paste0("tau2_r", seq_len(n_regions))
  if (spatial_field) hyper_names <- c(hyper_names, sigma_s2 = "Theta1 for spatial")

  inverse_exp <- function(theta) exp(-theta)
  summarise <- function(parameter, name, transform, decreasing) {
    summary <- screen_marginal_summary(screen_internal_marginal(fit, name),
                                       transform, decreasing)
    data.frame(parameter = parameter, t(summary), row.names = NULL)
  }

  table <- NULL
  for (j in seq_along(hyper_names)) {
    is_precision <- j <= n_regions
    transform <- if (is_precision) inverse_exp else function(theta) exp(2 * theta)
    table <- rbind(table, summarise(names(hyper_names)[j], hyper_names[[j]],
                                    transform, is_precision))
  }

  if (family == "nbinomial") {
    size <- "size for the nbinomial observations (1/overdispersion)"
    return(rbind(table,
                 summarise("nb_size0", size, exp, FALSE),
                 summarise("alpha0", size, inverse_exp, TRUE)))
  }

  rbind(table, summarise("sigma_e2", "Precision for the Gaussian observations",
                         inverse_exp, TRUE))
}


#' Regional mean contrasts with P(> 0) and local false sign rate (lfsr)
#' @noRd
screen_contrast_table <- function(fit, contrast_names) {
  do.call(rbind, lapply(contrast_names, function(name) {
    m <- fit$marginals.lincomb.derived[[name]]
    stopifnot(!is.null(m))
    p_positive <- 1 - inla.pmarginal(0, m)
    data.frame(quantity = name, t(screen_marginal_summary(m)),
               p_gt_zero = p_positive, lfsr = min(p_positive, 1 - p_positive))
  }))
}


#' Call each region "higher", "lower" or "undecided" relative to the rest
#'
#' Uses the `dmu_<r>_vs_rest` contrasts: "higher" when P(> 0) >= 0.95,
#' "lower" when P(> 0) <= 0.05.
#' @noRd
screen_domain_calls <- function(contrasts, region_labels) {
  n_regions <- length(region_labels)
  domains <- contrasts[match(paste0("dmu_", seq_len(n_regions), "_vs_rest"),
                             contrasts$quantity), ]
  domains$region <- region_labels
  domains$direction <- ifelse(domains$p_gt_zero >= .95, "higher",
                              ifelse(domains$p_gt_zero <= .05, "lower",
                                     "undecided"))
  domains
}


#' Per-gene slopes of the gene-specific covariates (NULL if there are none)
#' @noRd
screen_gene_slope_table <- function(result) {
  fit <- result$fit
  slope_names <- result$metadata$gene_covariate_names
  if (!length(slope_names)) return(NULL)

  do.call(rbind, lapply(seq_along(slope_names), function(j) {
    z <- fit$summary.random[[slope_names[j]]]
    stopifnot(nrow(z) == result$metadata$K)
    genes <- if (is.null(result$metadata$gene_labels)) seq_len(nrow(z)) else
      result$metadata$gene_labels
    data.frame(gene = genes,
               covariate = result$metadata$gene_covariate_columns[j],
               mean = z$mean,
               lower = z$`0.025quant`,
               median = z$`0.5quant`,
               upper = z$`0.975quant`)
  }))
}


#' Gene intercepts a_g (both modes) and composition slopes lambda_g
#' (adjusted mode only)
#' @noRd
screen_gene_coefficient_table <- function(result) {
  fit <- result$fit
  K <- result$metadata$K

  intercepts <- screen_summary_frame(
    fit$summary.lincomb.derived[paste0("gene_intercept_g", seq_len(K)), , drop = FALSE])
  table <- transform(intercepts, parameter = "a_g")
  genes <- result$metadata$gene_labels

  if (result$metadata$mode == "adjusted") {
    slopes <- screen_summary_frame(fit$summary.random$composition_slope)
    stopifnot(nrow(slopes) == K)
    slopes[, -1] <- 1 + slopes[, -1]  # lambda_g = 1 + v_g
    table <- rbind(table, transform(slopes, parameter = "lambda_g"))
    genes <- c(genes, result$metadata$gene_labels)
  }

  cbind(gene = genes, table)
}


#' Summarise a fit returned by screen_joint_fit()
#'
#' @param result Output of [screen_joint_fit()].
#' @return A list with
#'   * `pi`: directional scores and evidence from [pathway_pi_integrate()]
#'     (NULL when pi inference is disabled);
#'   * `parameters`: fixed effects, absolute regional means (named
#'     `adjusted_state_r<r>` in adjusted mode) and variance components
#'     (tau2_r, sigma_s2, and sigma_e2 or nb_size0 and alpha0);
#'   * `mean_contrasts`: all regional contrasts with `p_gt_zero` and `lfsr`;
#'   * `domains`: per-region calls from the region-vs-rest contrasts;
#'   * `gene_covariate_coefficients`: gene-specific slopes (or NULL);
#'   * `gene_coefficients`: adjusted mode, gene intercepts `a_g` and
#'     composition slopes `lambda_g` (NULL in observed mode);
#'   * the raw INLA marginals, `metadata`, INLA `warnings` and timing.
#' @export
screen_summarize <- function(result) {
  started <- proc.time()[["elapsed"]]
  fit <- result$fit
  metadata <- result$metadata

  parameters <- rbind(screen_fixed_effect_table(fit, metadata),
                      screen_variance_table(fit, metadata$R, metadata$family,
                                            isTRUE(metadata$spatial_field)))
  contrasts <- screen_contrast_table(fit, names(metadata$contrast_coefficients))
  domains <- screen_domain_calls(contrasts, metadata$region_labels)
  gene_slopes <- screen_gene_slope_table(result)

  # The intercept's marginal is that of alpha0 + mean(u) in both modes.
  fixed_marginals <- fit$marginals.fixed
  fixed_marginals[[metadata$intercept_name]] <-
    fit$marginals.lincomb.derived[["absolute_mu_r1"]]

  pi_result <- if (!is.null(result$joint)) {
    pathway_pi_integrate(result$joint,
                         delta = metadata$cfg$pi_delta,
                         level = metadata$cfg$pi_level,
                         nodes = metadata$cfg$pi_nodes,
                         max_nodes = metadata$cfg$pi_max_nodes,
                         tolerance = metadata$cfg$pi_tolerance)
  } else {
    NULL
  }

  list(
    pi = pi_result,
    joint = result$joint,
    parameters = parameters,
    gene_covariate_coefficients = gene_slopes,
    gene_coefficients = screen_gene_coefficient_table(result),
    mean_contrasts = contrasts,
    domains = domains,
    mean_contrast = contrasts[contrasts$quantity == "dmu_1_vs_rest", , drop = FALSE],
    # Separate hyperparameter marginals cannot identify a joint variance ratio
    # posterior. Do not multiply marginal densities or substitute modal ratios.
    variance_contrast = data.frame(
      quantity = "log_tau2_1_over_weighted_within_rest",
      mean = NA, lower = NA, median = NA, upper = NA,
      p_gt_zero = NA,
      status = "not_computed_from_separate_marginals"
    ),
    fixed_marginals = fixed_marginals,
    hyper_marginals = fit$marginals.hyperpar,
    internal_hyper_marginals = fit$internal.marginals.hyperpar,
    contrast_marginals = fit$marginals.lincomb.derived,
    metadata = metadata,
    warnings = fit$misc$warnings,
    postprocess_seconds = proc.time()[["elapsed"]] - started
  )
}


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

#' Write a screening summary to disk
#'
#' Files written to `out`:
#' * `screen_summary.rds`: the full summary object;
#' * `parameters.csv`, `mean_contrasts.csv`, `domains.csv`, `runtime.csv`;
#' * when pi inference is enabled: `pairwise_scores.csv`,
#'   `overall_scores.csv`, `pi_diagnostics.rds`, `pi_refinement.csv` and
#'   `joint_marginal_diagnostic.csv`;
#' * in adjusted mode: `gene_coefficients.csv` and
#'   `composition_diagnostics.csv`;
#' * `complete.txt`: a short completion marker.
#'
#' @param summary Output of [screen_summarize()].
#' @param out Output directory (created if needed).
#' @export
screen_write <- function(summary, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out_file <- function(name) file.path(out, name)
  write_table <- function(x, name) write.csv(x, out_file(name), row.names = FALSE)

  saveRDS(summary, out_file("screen_summary.rds"))

  if (!is.null(summary$pi)) {
    write_table(summary$pi$pairwise, "pairwise_scores.csv")
    write_table(summary$pi$overall, "overall_scores.csv")
    saveRDS(summary$pi$diagnostics, out_file("pi_diagnostics.rds"))
    write_table(summary$pi$diagnostics$refinement, "pi_refinement.csv")
    write_table(summary$joint$marginal_diagnostic, "joint_marginal_diagnostic.csv")
  }

  write_table(summary$parameters, "parameters.csv")
  if (!is.null(summary$gene_coefficients)) {
    write_table(summary$gene_coefficients, "gene_coefficients.csv")
    write_table(summary$metadata$composition_diagnostics,
                "composition_diagnostics.csv")
  }
  write_table(summary$mean_contrasts, "mean_contrasts.csv")
  write_table(summary$domains, "domains.csv")
  write_table(data.frame(step = c("fit", "marginal_summary"),
                         seconds = c(summary$metadata$fit_seconds,
                                     summary$postprocess_seconds)),
              "runtime.csv")

  pi_status <- if (is.null(summary$pi)) "Pi inference disabled." else
    paste("Pi numerical refinement converged:", summary$pi$diagnostics$converged)
  writeLines(c("CCD screening completed; config=FALSE; no posterior sampling.",
               paste("Gene effects:", summary$metadata$gene_dependence),
               paste("Free hyperparameters:", summary$metadata$hyperparameters),
               pi_status),
             out_file("complete.txt"))
}
