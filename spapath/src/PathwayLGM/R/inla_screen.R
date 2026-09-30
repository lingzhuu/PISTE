# CCD inference with optional low-dimensional deterministic pi integration.
#
# Model (Gaussian response, N spots x K genes, R regions):
#
#   y[i, g] = x_i' beta                 shared covariates, incl. the baseline
#                                       intercept (mean of region 1)
#           + mu_r(i)                   relative regional offset (mu_1 = 0)
#           + gene_r(i)[g]              region-specific iid gene effect,
#                                       variance tau2_r
#           + x_i' slope_g              optional gene-specific covariate slopes
#           + spatial_g(s_i)            replicated per-gene SPDE Matern field
#           + e[i, g]                   Gaussian noise, variance sigma_e2
#
# Free hyperparameters: R gene-effect precisions, one SPDE parameter and the
# Gaussian observation precision (R + 2 in total).


# -----------------------------------------------------------------------------
# Model fitting
# -----------------------------------------------------------------------------

#' Fit the joint Gaussian spatial pathway model with INLA
#'
#' @param Y Numeric N x K matrix of normalised expression (spots x genes).
#' @param xy N x 2 matrix of spot coordinates.
#' @param region Length-N vector of region labels (at least two regions).
#' @param mesh `fmesher`/INLA mesh for the SPDE field.
#' @param h Spatial length scale; the SPDE range parameter is fixed at
#'   kappa = `cfg$kappa_h / h`.
#' @param family Only `"gaussian"` is supported.
#' @param X Optional N x p covariate matrix. An `intercept` column is added if
#'   absent; it represents the baseline regional mean.
#' @param offset Optional known additive offset of length N.
#' @param cfg List of options overriding the defaults below:
#'   * `pi`, `pi_nodes`, `pi_max_nodes`, `pi_tolerance`: directional
#'     probability inference (see [pathway_pi_integrate()]).
#'   * `threads`: INLA `num.threads`.
#'   * `kappa_h`: kappa * h, must lie in \[2.3, 3.2\] (default `sqrt(8)`).
#'   * `strategy`: INLA approximation strategy.
#'   * `tau_prior`: prior on the gene-effect log precisions, one of
#'     `"normal"`, `"inla_default"`, `"pc"`, `"flat"`, with parameters
#'     `tau_logprecision_mean`, `tau_logprecision_sd`, `tau_pc_u`,
#'     `tau_pc_alpha`.
#'   * `other_weights`: weights of the non-baseline regions in the
#'     "region 1 vs rest" contrast (default: proportional to spot counts).
#'   * `gene_specific_covariates`, `gene_covariate_prior_sd`: columns of `X`
#'     that get independent per-gene slopes with a fixed Normal prior.
#'   * `baseline_region`: label of the baseline region (default: the first
#'     label found).
#' @param .fixed_theta Internal test hook: fix the hyperparameters.
#' @param .int_strategy INLA integration strategy, `"ccd"` or `"eb"`.
#' @return A list with the INLA `fit`, the extracted `joint` mixture (or NULL
#'   when `cfg$pi` is FALSE) and a `metadata` list.
#' @export
screen_joint_fit <- function(Y, xy, region, mesh, h, family = "gaussian",
                             X = NULL, offset = NULL,
                             cfg = list(), .fixed_theta = NULL,
                             .int_strategy = "ccd") {
  # --- Reject legacy / unsupported settings ----------------------------------
  if ("gene_cs_rho" %in% names(cfg))
    stop("gene_cs_rho has been removed; region-specific gene effects are conditionally iid. Remove this legacy setting.")
  family <- match.arg(family, "gaussian")
  if (any(grepl("^nb_", names(cfg)))) stop("Unsupported configuration in this Gaussian-only package.")
  stopifnot(.int_strategy %in% c("ccd", "eb"))
  if (!is.null(cfg$design_size) || (!is.null(cfg$integration) && cfg$integration != "ccd"))
    stop("Production screening uses CCD; remove legacy design_size/non-CCD settings.")

  # --- Default configuration -------------------------------------------------
  cfg <- modifyList(list(pi = TRUE, pi_nodes = 4096L, pi_max_nodes = 65536L,
                         pi_tolerance = .002,
                         threads = "2:1", kappa_h = sqrt(8),
                         strategy = "gaussian",
                         tau_prior = "normal",
                         tau_logprecision_mean = 0, tau_logprecision_sd = 3,
                         tau_pc_u = 2, tau_pc_alpha = .05, other_weights = NULL,
                         gene_specific_covariates = character(),
                         gene_covariate_prior_sd = 10,
                         baseline_region = NULL),
                    cfg)
  if (isTRUE(cfg$pi) && (.int_strategy != "ccd" || !is.null(.fixed_theta)))
    stop("Pi inference requires unfixed CCD hyperparameter integration.")

  # --- Dimensions and regions (baseline region first) ------------------------
  N <- nrow(Y)
  K <- ncol(Y)
  found_labels <- unique(as.character(region))
  baseline <- if (is.null(cfg$baseline_region)) found_labels[1L] else
    as.character(cfg$baseline_region)
  stopifnot(length(baseline) == 1L, baseline %in% found_labels)
  labels <- c(baseline, setdiff(found_labels, baseline))
  rr <- match(as.character(region), labels)
  R <- length(labels)

  stopifnot(N == nrow(xy), ncol(xy) == 2L, length(region) == N, R >= 2L, K >= 2L,
            all(is.finite(Y)), all(is.finite(xy)), !anyNA(region),
            is.finite(h), h > 0,
            cfg$kappa_h >= 2.3, cfg$kappa_h <= 3.2,
            is.finite(cfg$tau_logprecision_mean),
            is.finite(cfg$tau_logprecision_sd),
            cfg$tau_logprecision_sd > 0,
            cfg$tau_prior %in% c("normal", "inla_default", "pc", "flat"))

  # --- Covariates ------------------------------------------------------------
  if (is.null(X)) X <- matrix(numeric(), N, 0L)
  X <- as.matrix(X)
  stopifnot(nrow(X) == N, all(is.finite(X)))
  if (!"intercept" %in% colnames(X))
    X <- cbind(intercept = rep(1, N), X)
  # Exactly one constant-one column, and it must be the intercept.
  stopifnot(all(X[, "intercept"] == 1),
            sum(colSums(X == 1) == N) == 1L)
  if (ncol(X) > 0L) {
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    if (anyDuplicated(colnames(X)) || any(!nzchar(colnames(X))))
      stop("X must have unique, non-empty column names.")
    x_names <- paste0("cov_", make.names(colnames(X), unique = TRUE))
  } else x_names <- character()

  # Split covariates into shared fixed effects and gene-specific slopes.
  gene_x <- match(cfg$gene_specific_covariates, colnames(X))
  stopifnot(!anyNA(gene_x), !anyDuplicated(gene_x),
            !(match("intercept", colnames(X)) %in% gene_x),
            length(cfg$gene_covariate_prior_sd) == 1L,
            is.finite(cfg$gene_covariate_prior_sd),
            cfg$gene_covariate_prior_sd > 0)
  shared_x <- setdiff(seq_len(ncol(X)), gene_x)
  gene_x_names <- if (length(gene_x)) paste0("slope_", x_names[gene_x]) else character()
  gene_x_weights <- if (length(gene_x)) paste0(gene_x_names, "_weight") else character()

  if (is.null(offset)) offset <- rep(0, N)
  stopifnot(length(offset) == N, all(is.finite(offset)))

  # --- SPDE spatial field ----------------------------------------------------
  # kappa is fixed by h; the single free parameter controls the marginal
  # variance of the field.
  kappa <- cfg$kappa_h / h
  spde <- inla.spde2.matern(mesh, alpha = 2, constr = FALSE,
                            B.tau = matrix(c(-log(sqrt(4 * pi) * kappa), -1), 1, 2),
                            B.kappa = matrix(c(log(kappa), 0), 1, 2),
                            theta.prior.mean = 0, theta.prior.prec = 1)
  stopifnot(spde$n.theta == 1L)

  # Long format: observation (site, gene), one field replicate per gene.
  site <- rep(seq_len(N), K)
  gene <- rep(seq_len(K), each = N)
  rlong <- rr[site]
  idx <- inla.spde.make.index("spatial", mesh$n, n.repl = K)
  A <- inla.spde.make.A(mesh, loc = xy, index = site, repl = gene, n.repl = K)
  nz <- summary(A)
  stopifnot(all((nz$j - 1L) %/% mesh$n + 1L == gene[nz$i]),
            max(abs(rowSums(A) - 1)) < 1e-9)

  # --- Effects table ---------------------------------------------------------
  # The intercept in X is the baseline regional mean. Other mu coefficients
  # are offsets from that baseline, with delta_1 fixed at zero.
  mu_names <- paste0("mu_r", seq.int(2L, R))
  intercept_name <- x_names[match("intercept", colnames(X))]
  gene_names <- paste0("gene_r", seq_len(R))

  effects <- data.frame(row = seq_len(N * K))
  for (j in shared_x) effects[[x_names[j]]] <- X[site, j]
  for (j in seq_along(gene_x)) {
    effects[[gene_x_names[j]]] <- gene
    effects[[gene_x_weights[j]]] <- X[site, gene_x[j]]
  }
  for (r in seq_len(R)) {
    # Regional offset indicator, and gene index active only in region r.
    if (r > 1L) effects[[mu_names[r - 1L]]] <- as.numeric(rlong == r)
    effects[[gene_names[r]]] <- ifelse(rlong == r, gene, NA_integer_)
  }
  effects$row <- NULL
  stopifnot(all(rowSums(!is.na(effects[gene_names])) == 1L))

  stk <- inla.stack(list(y = as.vector(Y)), A = list(1, A),
                    effects = list(effects, idx),
                    compress = FALSE, remove.unused = FALSE)

  # --- Priors and formula ----------------------------------------------------
  # Proper approximation to a log-flat scale prior. INLA's internal iid
  # hyperparameter is log precision; param[2] is Normal precision.
  iid_hyper <- switch(cfg$tau_prior,
    normal = list(prec = list(prior = "normal",
                              param = c(cfg$tau_logprecision_mean,
                                        1 / cfg$tau_logprecision_sd^2))),
    inla_default = list(prec = list(prior = "loggamma", param = c(1, 5e-5))),
    pc = list(prec = list(prior = "pc.prec", param = c(cfg$tau_pc_u, cfg$tau_pc_alpha))),
    flat = list(prec = list(prior = "flat")))
  terms <- paste0('f(', gene_names, ',model="iid",n=K,constr=FALSE,hyper=iid_hyper)')

  # Proper independent per-gene slopes with fixed prior precision: no extra
  # hyperparameters and no dense gene-by-covariate dummy matrix.
  slope_hyper <- list(prec = list(initial = log(1 / cfg$gene_covariate_prior_sd^2),
                                  fixed = TRUE))
  slope_terms <- if (length(gene_x)) paste0('f(', gene_x_names, ',', gene_x_weights,
    ',model="iid",n=K,constr=FALSE,hyper=slope_hyper)') else character()

  formula <- as.formula(paste("y ~ -1 +", paste(c(
    x_names[shared_x], mu_names, terms, slope_terms,
    'f(spatial,model=spde,replicate=spatial.repl,nrep=K,constr=FALSE)'),
    collapse = " + ")))

  # --- Linear combinations (regional contrasts) ------------------------------
  counts <- as.numeric(table(factor(rr, levels = seq_len(R))))
  weights <- vector("list", R)
  lc <- list()
  contrast_info <- list()

  # Region r vs a weighted average of all other regions.
  for (r in seq_len(R)) {
    rest <- setdiff(seq_len(R), r)
    a <- counts[rest] / sum(counts[rest])
    if (r == 1L && !is.null(cfg$other_weights)) {
      a <- cfg$other_weights
      stopifnot(length(a) == R - 1L, all(is.finite(a)), all(a >= 0), sum(a) > 0)
      a <- a / sum(a)
    }
    z <- numeric(R)
    z[r] <- 1
    z[rest] <- -a
    weights[[r]] <- a
    name <- paste0("dmu_", r, "_vs_rest")
    lc <- c(lc, setNames(do.call(inla.make.lincomb,
                                 as.list(setNames(z[-1L], mu_names))), name))
    contrast_info[[name]] <- z
  }

  # All pairwise region differences.
  for (pair in combn(seq_len(R), 2, simplify = FALSE)) {
    z <- numeric(R)
    z[pair] <- c(1, -1)
    name <- paste0("dmu_", pair[1], "_vs_", pair[2])
    lc <- c(lc, setNames(do.call(inla.make.lincomb,
                                 as.list(setNames(z[-1L], mu_names))), name))
    contrast_info[[name]] <- z
  }

  # Report absolute fitted means as diagnostics while keeping mu_r2,...,mu_rR
  # as the primary relative regional coefficients.
  for (r in seq_len(R)) {
    z <- c(1, as.numeric(seq.int(2L, R) == r))
    names(z) <- c(intercept_name, mu_names)
    lc <- c(lc, setNames(do.call(inla.make.lincomb, as.list(z)),
                         paste0("absolute_mu_r", r)))
  }

  # --- Run INLA --------------------------------------------------------------
  fam <- list(hyper = list(prec = list(prior = "pc.prec", param = c(2, .05))))
  call <- list(formula = formula, family = family,
               data = inla.stack.data(stk), lincomb = lc,
               control.family = fam,
               control.fixed = list(mean = 0, prec = .01),
               control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
               control.compute = list(config = isTRUE(cfg$pi), dic = FALSE,
                                      waic = FALSE, cpo = FALSE),
               control.inla = list(strategy = cfg$strategy,
                                   int.strategy = .int_strategy),
               inla.mode = "experimental", num.threads = cfg$threads,
               verbose = FALSE)
  if (any(offset != 0)) call$offset <- offset[site]
  # Structural test hook only: conditions known hyperparameters, no sampling.
  if (!is.null(.fixed_theta))
    call$control.mode <- list(theta = .fixed_theta, fixed = TRUE, restart = FALSE)

  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(inla, call)
  elapsed <- proc.time()[["elapsed"]] - t0

  # --- Extract the low-dimensional joint posterior, then drop configs --------
  joint <- NULL
  config_bytes <- 0
  if (isTRUE(cfg$pi)) {
    config_bytes <- as.numeric(object.size(fit$misc$configs))
    joint <- pi_extract_joint(fit, labels, counts, intercept_name, mu_names)
    fit$misc$configs <- NULL
  }
  stopifnot(is.null(fit$misc$configs),
            nrow(fit$summary.random$spatial) == mesh$n * K,
            all(vapply(fit$summary.random[gene_names], nrow, integer(1)) == K))
  if (is.null(.fixed_theta))
    stopifnot(length(fit$mode$theta) == R + 2L,
              nrow(fit$summary.hyperpar) == R + 2L)

  # --- Result ----------------------------------------------------------------
  list(fit = fit, joint = joint,
       metadata = list(
         N = N, K = K, R = R,
         region_labels = labels, region_spots = counts,
         other_weights = weights[[1]], rest_weights = weights,
         h = h, kappa0 = kappa, rho0 = sqrt(8) / kappa,
         gene_dependence = "conditional_iid", model_version = "0.5.0",
         cfg = cfg, family = family, fit_seconds = elapsed,
         hyperparameters = R + 2L,
         tau_prior = list(
           internal_parameter = "log precision", prior = cfg$tau_prior,
           mean = if (cfg$tau_prior == "normal") cfg$tau_logprecision_mean else NA_real_,
           sd = if (cfg$tau_prior == "normal") cfg$tau_logprecision_sd else NA_real_,
           parameters = iid_hyper$prec$param),
         baseline_region = baseline, intercept_name = intercept_name,
         regional_offset_names = mu_names,
         covariate_names = x_names[shared_x],
         gene_covariate_names = gene_x_names,
         gene_covariate_columns = colnames(X)[gene_x],
         gene_labels = colnames(Y),
         offset_convention = if (any(offset != 0))
           "known Gaussian additive offset" else "zero (SCT/normalized Gaussian response)",
         integration = .int_strategy, config = FALSE,
         fit_config = isTRUE(cfg$pi),
         pi_enabled = isTRUE(cfg$pi),
         config_bytes_before_extraction = config_bytes,
         joint_bytes = as.numeric(object.size(joint)),
         posterior_sampling = FALSE,
         contrast_coefficients = contrast_info,
         INLA_version = as.character(packageVersion("INLA")),
         formula = deparse(formula)))
}


# -----------------------------------------------------------------------------
# Posterior summaries
# -----------------------------------------------------------------------------

#' Mean and 95% interval of `transform(x)` for an INLA marginal
#'
#' Strictly positive marginals are routed to the checked log-scale integrator
#' `positive_marginal_summary()`.
#' @noRd
screen_marginal_summary <- function(m, transform = identity, decreasing = FALSE) {
  if (all(m[, 1] > 0)) return(positive_marginal_summary(m, transform, decreasing))
  qs <- inla.qmarginal(c(.025, .5, .975), m)
  if (decreasing) qs <- rev(qs)
  c(mean = inla.emarginal(transform, m), lower = transform(qs[1]),
    median = transform(qs[2]), upper = transform(qs[3]))
}


#' Summarise a fit returned by screen_joint_fit()
#'
#' Produces parameter summaries (baseline intercept, relative and absolute
#' regional means, covariates, variance components), regional mean contrasts
#' with P(> 0) and local false sign rate, domain calls, optional gene-specific
#' slopes and, when enabled, directional pi scores and evidence.
#'
#' Domains are called "higher" when P(contrast > 0) >= 0.95, "lower" when
#' <= 0.05, and "undecided" otherwise.
#'
#' @param result Output of [screen_joint_fit()].
#' @return A list of summary tables and marginals.
#' @export
screen_summarize <- function(result) {
  t0 <- proc.time()[["elapsed"]]
  fit <- result$fit
  R <- result$metadata$R
  stopifnot(is.null(fit$misc$configs))

  # --- Fixed effects ---------------------------------------------------------
  fixed <- fit$summary.fixed[c(result$metadata$intercept_name,
                               result$metadata$regional_offset_names,
                               setdiff(result$metadata$covariate_names,
                                       result$metadata$intercept_name)), ,
                             drop = FALSE]
  pars <- data.frame(parameter = rownames(fixed), mean = fixed$mean,
                     lower = fixed$`0.025quant`, median = fixed$`0.5quant`,
                     upper = fixed$`0.975quant`, row.names = NULL)
  pars$parameter[pars$parameter == result$metadata$intercept_name] <-
    "baseline_intercept"
  pars$parameter[pars$parameter %in% result$metadata$regional_offset_names] <-
    paste0("relative_", pars$parameter[pars$parameter %in%
                                         result$metadata$regional_offset_names])

  # Absolute regional means (from the linear combinations).
  absolute <- fit$summary.lincomb.derived[paste0("absolute_mu_r", seq_len(R)), ,
                                          drop = FALSE]
  pars <- rbind(pars, data.frame(parameter = rownames(absolute),
                                 mean = absolute$mean,
                                 lower = absolute$`0.025quant`,
                                 median = absolute$`0.5quant`,
                                 upper = absolute$`0.975quant`,
                                 row.names = NULL))

  # --- Variance components ---------------------------------------------------
  # tau2_r = 1 / precision of gene_r; sigma_s2 = exp(2 * Theta1) for the field.
  keys <- c(paste0("Precision for gene_r", seq_len(R)), "Theta1 for spatial")
  names(keys) <- c(paste0("tau2_r", seq_len(R)), "sigma_s2")
  for (j in seq_along(keys)) {
    fun <- if (j <= R) function(x) 1 / x else function(x) exp(2 * x)
    pars <- rbind(pars, data.frame(
      parameter = names(keys)[j],
      t(screen_marginal_summary(fit$marginals.hyperpar[[keys[j]]], fun, j <= R)),
      row.names = NULL))
  }
  pars <- rbind(pars, data.frame(parameter = "sigma_e2", t(screen_marginal_summary(
    fit$marginals.hyperpar[["Precision for the Gaussian observations"]],
    function(x) 1 / x, TRUE))))

  # --- Regional mean contrasts -----------------------------------------------
  contrasts <- do.call(rbind, lapply(names(result$metadata$contrast_coefficients), function(name) {
    m <- fit$marginals.lincomb.derived[[name]]
    stopifnot(!is.null(m))
    p <- 1 - inla.pmarginal(0, m)
    data.frame(quantity = name, t(screen_marginal_summary(m)),
               p_gt_zero = p, lfsr = min(p, 1 - p))
  }))

  # Domain calls from the "region vs rest" contrasts.
  domain <- contrasts[match(paste0("dmu_", seq_len(R), "_vs_rest"), contrasts$quantity), ]
  domain$region <- result$metadata$region_labels
  domain$direction <- ifelse(domain$p_gt_zero >= .95, "higher",
                             ifelse(domain$p_gt_zero <= .05, "lower", "undecided"))

  # --- Gene-specific covariate slopes ----------------------------------------
  gene_slopes <- if (length(result$metadata$gene_covariate_names)) do.call(rbind,
    lapply(seq_along(result$metadata$gene_covariate_names), function(j) {
      z <- fit$summary.random[[result$metadata$gene_covariate_names[j]]]
      stopifnot(nrow(z) == result$metadata$K)
      data.frame(gene = if (is.null(result$metadata$gene_labels)) seq_len(nrow(z)) else
                   result$metadata$gene_labels,
                 covariate = result$metadata$gene_covariate_columns[j],
                 mean = z$mean,
                 lower = z$`0.025quant`, median = z$`0.5quant`,
                 upper = z$`0.975quant`)
    })) else NULL

  # --- Directional pi scores -------------------------------------------------
  pi_result <- if (!is.null(result$joint)) pathway_pi_integrate(result$joint,
    nodes = result$metadata$cfg$pi_nodes,
    max_nodes = result$metadata$cfg$pi_max_nodes,
    tolerance = result$metadata$cfg$pi_tolerance) else NULL

  list(pi = pi_result,
       joint = result$joint,
       parameters = pars,
       gene_covariate_coefficients = gene_slopes,
       mean_contrasts = contrasts,
       domains = domain,
       mean_contrast = contrasts[contrasts$quantity == "dmu_1_vs_rest", , drop = FALSE],
       # Separate hyperparameter marginals cannot identify a joint variance ratio
       # posterior. Do not multiply marginal densities or substitute modal ratios.
       variance_contrast = data.frame(
         quantity = "log_tau2_1_over_weighted_within_rest",
         mean = NA_real_, lower = NA_real_, median = NA_real_, upper = NA_real_,
         p_gt_zero = NA_real_,
         status = "not_computed_from_separate_marginals"),
       fixed_marginals = fit$marginals.fixed,
       hyper_marginals = fit$marginals.hyperpar,
       contrast_marginals = fit$marginals.lincomb.derived,
       metadata = result$metadata,
       warnings = fit$misc$warnings,
       postprocess_seconds = proc.time()[["elapsed"]] - t0)
}


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

#' Write a screening summary to disk
#'
#' Writes `screen_summary.rds`, CSV tables (parameters, mean contrasts,
#' domains, runtime and, when pi inference is enabled, pairwise/overall scores
#' and diagnostics) and a short `complete.txt` marker into `out`.
#'
#' @param summary Output of [screen_summarize()].
#' @param out Output directory (created if needed).
#' @export
screen_write <- function(summary, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(summary, file.path(out, "screen_summary.rds"))

  if (!is.null(summary$pi)) {
    write.csv(summary$pi$pairwise, file.path(out, "pairwise_scores.csv"), row.names = FALSE)
    write.csv(summary$pi$overall, file.path(out, "overall_scores.csv"), row.names = FALSE)
    saveRDS(summary$pi$diagnostics, file.path(out, "pi_diagnostics.rds"))
    write.csv(summary$pi$diagnostics$refinement,
              file.path(out, "pi_refinement.csv"), row.names = FALSE)
    write.csv(summary$joint$marginal_diagnostic,
              file.path(out, "joint_marginal_diagnostic.csv"), row.names = FALSE)
  }

  write.csv(summary$parameters, file.path(out, "parameters.csv"), row.names = FALSE)
  write.csv(summary$mean_contrasts, file.path(out, "mean_contrasts.csv"), row.names = FALSE)
  write.csv(summary$domains, file.path(out, "domains.csv"), row.names = FALSE)
  write.csv(data.frame(step = c("fit", "marginal_summary"),
                       seconds = c(summary$metadata$fit_seconds,
                                   summary$postprocess_seconds)),
            file.path(out, "runtime.csv"), row.names = FALSE)

  writeLines(c("CCD screening completed; config=FALSE; no posterior sampling.",
               paste("Gene effects:", summary$metadata$gene_dependence),
               paste("Free hyperparameters:", summary$metadata$hyperparameters),
               if (is.null(summary$pi)) "Pi inference disabled." else
                 paste("Pi numerical refinement converged:",
                       summary$pi$diagnostics$converged)),
             file.path(out, "complete.txt"))
}
