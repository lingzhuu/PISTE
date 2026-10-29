# =============================================================================
# closed_form_screen.R
#
# Closed-form fit of the observed-mode model for one section. The conditional
# posterior given the hyperparameters is Gaussian in closed form and the
# hyperparameters are integrated on a grid; no INLA. The result is summarised
# into the fields of screen_summarize(), so screen_write() writes the same
# files. Adjacent sections are fitted jointly by joint_sections.R.
#
# Model (N labelled spots, K genes, R regions, region 1 = baseline)
#
#   y[i, g] = alpha0 + u_g + mu_r(i) + gene_r(i)[g] + e[i, g],   mu_1 = 0
#
#   alpha0       ~ N(0, 100 (1 - 1/K))     so that alpha0 + mean(u) ~ N(0, 100)
#   mu_r, r > 1  ~ N(0, 100)               (INLA control.fixed prec .01)
#   u_g          ~ N(0, 100) iid           (gene intercept, fixed precision)
#   gene_r[g]    ~ N(0, tau2_r) iid
#   e[i, g]      ~ N(0, sigma_e2)
#   log(1 / tau2_r)   ~ N(tau_logprecision_mean, tau_logprecision_sd^2)
#   1 / sigma_e2      ~ pc.prec(u = 2, alpha = .05)
#
# Hyperparameters theta = (log 1/sigma_e2, log 1/tau2_1, ..., log 1/tau2_R).
#
# Exact conditional posterior. Given theta, integrating out u_g and gene_r[g]
# leaves every gene's vector of regional means ybar_g (length R) Gaussian,
#
#   ybar_g | alpha0, mu ~ N(X beta, S),  beta = (alpha0, mu_2, ..., mu_R),
#   X = [1 | 0; I_(R-1)],  S = 100 * 11' + diag(tau2_r + sigma_e2 / n_r),
#
# with the same S for every gene (every gene is observed at every labelled
# spot). The data enter only through the per-gene, per-region means, the spot
# counts n_r and the total within-region sum of squares, so each theta costs
# O(K R^2). beta has a Gaussian posterior in closed form; the canonical
# regional states (alpha0 + mean(u), mu_2, ..., mu_R), the gene intercepts
# a_g = alpha0 + u_g and all regional contrasts are linear in (beta, u) and
# Gaussian given theta. The marginal likelihood p(y | theta) is exact.
#
# Integration over theta. The log posterior of theta is maximised (BFGS) and
# its Hessian taken numerically. theta is then integrated on a regular grid in
# the standardised coordinates z of that Gaussian approximation,
# theta = mode + V Lambda^(-1/2) z, with cfg$closed_form_grid_points points per
# axis over +- cfg$closed_form_grid_halfwidth standard deviations. Grid weights
# are proportional to exp(log posterior); points with relative weight below
# cfg$closed_form_weight_threshold are dropped. INLA instead uses a central
# composite design and stops its mode search earlier. On 240 DLPFC random null
# sets the closed form evaluated at INLA's own design points reproduces INLA's
# pi and evidence to < 5e-6; on the grid it differs from INLA by at most 0.003
# in pi and 0.03 in evidence (calls identical on 240/240), on 200 simulated sets
# by at most 0.009 in pi and 0.014 in the evidences (calls identical on
# 200/200). Refining the grid up to 31 points per axis over +-7 standard
# deviations moves the closed form by < 0.003, so the grid is the more accurate
# of the two integrations.
#
# Posterior summaries. Latent quantities: mixtures of Gaussians over the grid
# (mean exact; quantiles by inverting the mixture CDF). Hyperparameters:
# weighted grid distribution of theta_j (mean of the transformed value; 2.5, 50
# and 97.5% quantiles = first grid value where the weighted CDF reaches p).
# Pi and evidence: pathway_pi_integrate() on the same mixture format as
# pi_extract_joint() (mu, covariance, tau2, kappa2 per grid point;
# kappa2_rs = sigma_e2 (1/n_r + 1/n_s)).
#
# Scope: observed mode, intercept-only design, no offset, tau_prior = "normal".
# The grid has grid_points^(R+1) points, so closed_form_joint_fit() stops when
# it would exceed cfg$closed_form_max_points.
# =============================================================================


#' Sufficient statistics: spot counts, per-gene regional means, within SS
#' @noRd
cf_sufficient_statistics <- function(Y, region_index, R) {
  n <- tabulate(region_index, R)
  if (any(n == 0)) stop("Every region needs at least one spot.")
  sums <- rowsum(Y, region_index, reorder = TRUE)            # R x K
  means <- sums / n
  within <- sum(Y^2) - sum(sums * means)                     # sum over g, r of SS - n ybar^2
  list(n = n, ybar = t(means), within = within, N = nrow(Y), K = ncol(Y), R = R)
}


#' Everything at one theta: log posterior and the conditional Gaussians
#' @noRd
cf_at_theta <- function(theta, stats, prior) {
  K <- stats$K; R <- stats$R; n <- stats$n
  sigma_e2 <- exp(-theta[1]); tau2 <- exp(-theta[-1])
  X <- diag(R); X[, 1] <- 1                                   # beta -> regional means
  S <- prior$u_var + diag(tau2 + sigma_e2 / n, R)             # Cov(ybar_g | beta)
  S_chol <- chol(S)
  S_inv <- chol2inv(S_chol)
  P_beta <- diag(c(1 / prior$alpha0_var, rep(1 / prior$mu_var, R - 1)), R)
  sum_ybar <- colSums(stats$ybar)
  XtSinv <- crossprod(X, S_inv)
  Q_beta <- P_beta + K * XtSinv %*% X
  Q_chol <- chol(Q_beta)
  b <- XtSinv %*% sum_ybar
  m_beta <- backsolve(Q_chol, forwardsolve(t(Q_chol), b))
  C_beta <- chol2inv(Q_chol)

  # log p(y | theta): gene-level means part + within-region part
  quad <- sum((stats$ybar %*% S_inv) * stats$ybar)
  log_lik <- -K * R / 2 * log(2 * pi) - K * sum(log(diag(S_chol))) - quad / 2 +
    sum(log(diag(P_beta))) / 2 - sum(log(diag(Q_chol))) + sum(b * m_beta) / 2 -
    K * sum(n - 1) / 2 * log(2 * pi * sigma_e2) - K / 2 * sum(log(n)) -
    stats$within / (2 * sigma_e2)
  lambda <- -log(prior$pc_alpha) / prior$pc_u                 # pc.prec on 1/sigma_e2, log scale
  log_prior <- log(lambda / 2) - lambda * exp(-theta[1] / 2) - theta[1] / 2 +
    sum(dnorm(theta[-1], prior$tau_mean, prior$tau_sd, log = TRUE))

  # u_g | beta, ybar_g: E = h'(ybar_g - X beta), Var = u_var - u_var^2 1'S^-1 1
  h <- prior$u_var * rowSums(S_inv)
  var_u <- prior$u_var - prior$u_var * sum(h)
  # canonical states s = (alpha0 + mean(u), mu_2..mu_R) = L beta + c (+ noise in s_1)
  L <- diag(R); L[1, ] <- L[1, ] - as.numeric(crossprod(h, X))
  c0 <- c(sum(h * sum_ybar) / K, rep(0, R - 1))
  s_mean <- as.numeric(L %*% m_beta) + c0
  s_cov <- L %*% C_beta %*% t(L); s_cov[1, 1] <- s_cov[1, 1] + var_u / K
  # gene intercepts a_g = alpha0 + u_g
  a_vec <- c(1, rep(0, R - 1)) - as.numeric(crossprod(X, h))
  a_mean <- m_beta[1] + as.numeric(stats$ybar %*% h) - sum(h * (X %*% m_beta))
  a_var <- as.numeric(crossprod(a_vec, C_beta %*% a_vec)) + var_u
  list(log_post = log_lik + log_prior, s_mean = s_mean, s_cov = (s_cov + t(s_cov)) / 2,
       a_mean = a_mean, a_var = a_var, sigma_e2 = sigma_e2, tau2 = tau2)
}


#' Mixture-of-Gaussians summary of a scalar: mean, 2.5/50/97.5% quantiles, P(> 0)
#' @noRd
cf_mixture_summary <- function(m, s, w) {
  cdf <- function(q) sum(w * pnorm((q - m) / s))
  lo <- min(m - 8 * s); hi <- max(m + 8 * s)
  qf <- function(p) uniroot(function(q) cdf(q) - p, c(lo, hi), tol = 1e-10)$root
  c(mean = sum(w * m), lower = qf(.025), median = qf(.5), upper = qf(.975),
    p_gt_zero = 1 - cdf(0))
}


#' Weighted-grid summary of a hyperparameter on the variance scale
#' @noRd
cf_hyper_summary <- function(theta_j, w, transform, decreasing) {
  o <- order(theta_j); x <- theta_j[o]; cw <- cumsum(w[o])
  qf <- function(p) x[min(which(cw >= p))]
  q <- c(qf(.025), qf(.5), qf(.975)); if (decreasing) q <- rev(q)
  c(mean = sum(w * transform(theta_j)), lower = transform(q[1]),
    median = transform(q[2]), upper = transform(q[3]))
}


#' Fit the observed-mode model of one section in closed form
#'
#' The observed mode: a Gaussian model of normalised expression with a gene
#' intercept, regional means and gene-level regional effects. Given the
#' hyperparameters the posterior is Gaussian, and the R + 1 hyperparameters
#' (the error variance and the gene-effect variance of every region) are
#' integrated on a regular grid around their posterior mode. The fit takes a
#' fraction of a second. Adjacent sections are fitted jointly by
#' [joint_sections_fit()].
#'
#' @param Y Numeric matrix of normalized expression, spots x genes.
#' @param region Factor or character vector of region labels per spot.
#' @param X Optional design matrix; only an intercept column is supported.
#' @param offset Optional offset; must be absent or zero.
#' @param mode Must be `"observed"`.
#' @param cfg Options as for [screen_joint_fit()]. The grid is set by
#'   `closed_form_grid_points` (points per hyperparameter axis),
#'   `closed_form_grid_halfwidth` (half width in posterior standard
#'   deviations), `closed_form_weight_threshold` (relative weight below which
#'   a grid point is dropped), `closed_form_max_points` (largest grid
#'   allowed) and `closed_form_hessian_floor` (smallest Hessian eigenvalue
#'   used to scale the grid; `metadata$hessian_floored` records whether it
#'   was applied).
#' @return A list with `fit = NULL`, `joint` (format of the `joint` element
#'   of [screen_joint_fit()], or `NULL` when `cfg$pi = FALSE`), `metadata`
#'   (with `engine = "closed_form"` and the grid used) and `closed_form` (the
#'   conditional Gaussians, grid points and weights). Pass it to
#'   [screen_summarize()].
#' @export
closed_form_joint_fit <- function(Y, region, X = NULL, offset = NULL,
                                  mode = "observed", cfg = list()) {
  started <- proc.time()[["elapsed"]]
  mode <- match.arg(mode, "observed")
  screen_check_legacy_cfg(cfg)
  screen_check_mode_cfg(cfg, mode)
  screen_check_integration_cfg(cfg)
  cfg <- modifyList(screen_default_cfg(mode), cfg)
  if (cfg$tau_prior != "normal") stop("The closed form implements tau_prior = 'normal' only.")
  if (length(cfg$gene_specific_covariates)) stop("The closed form has no gene-specific covariates.")
  if (!is.null(X) && !(ncol(X) == 1 && all(X == 1))) stop("The closed form supports X = intercept only.")
  if (!is.null(offset) && any(offset != 0)) stop("The closed form supports no offset.")
  stopifnot(is.numeric(cfg$closed_form_grid_points), cfg$closed_form_grid_points >= 3,
            cfg$closed_form_grid_halfwidth > 0, cfg$closed_form_weight_threshold >= 0)

  Y <- as.matrix(Y)
  N <- nrow(Y); K <- ncol(Y)
  labels <- screen_region_labels(region, cfg$baseline_region)
  region_index <- match(as.character(region), labels)
  R <- length(labels)
  stopifnot(length(region) == N, R >= 2, K >= 2, all(is.finite(Y)), !anyNA(region_index))
  d <- R + 1
  n_points <- cfg$closed_form_grid_points^d
  if (n_points > cfg$closed_form_max_points) {
    stop("The hyperparameter grid would have ", n_points, " points (R = ", R,
         "); reduce closed_form_grid_points or raise closed_form_max_points.")
  }

  stats <- cf_sufficient_statistics(Y, region_index, R)
  prior <- list(alpha0_var = 100 * (1 - 1 / K), mu_var = 100, u_var = 100,
                pc_u = 2, pc_alpha = .05,
                tau_mean = cfg$tau_logprecision_mean, tau_sd = cfg$tau_logprecision_sd)

  # --- Mode and curvature of the hyperparameter posterior ------------------------
  neg_log_post <- function(theta) {
    v <- tryCatch(-cf_at_theta(theta, stats, prior)$log_post, error = function(e) Inf)
    if (is.finite(v)) v else 1e100
  }
  start <- c(-log(stats$within / (N * K - R * K)), rep(cfg$tau_logprecision_mean, R))
  opt <- optim(start, neg_log_post, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 1000))
  eig <- eigen(optimHess(opt$par, neg_log_post), symmetric = TRUE)
  # A variance at the boundary (tau2_r -> 0) can leave a flat or slightly
  # indefinite curvature; floor the eigenvalues so that the grid spans that
  # axis widely. The prior keeps the posterior proper.
  hessian_floored <- any(eig$values < cfg$closed_form_hessian_floor)
  eig$values <- pmax(eig$values, cfg$closed_form_hessian_floor)
  to_theta <- eig$vectors %*% diag(1 / sqrt(eig$values), d)

  # --- Grid in standardised coordinates ------------------------------------------
  axis <- seq(-cfg$closed_form_grid_halfwidth, cfg$closed_form_grid_halfwidth,
              length.out = cfg$closed_form_grid_points)
  Z <- as.matrix(expand.grid(rep(list(axis), d)))
  thetas <- t(opt$par + to_theta %*% t(Z))
  points <- lapply(seq_len(nrow(thetas)), function(j) cf_at_theta(thetas[j, ], stats, prior))
  log_post <- vapply(points, `[[`, numeric(1), "log_post")
  w <- exp(log_post - max(log_post)); w <- w / sum(w)
  keep <- w >= cfg$closed_form_weight_threshold
  points <- points[keep]; thetas <- thetas[keep, , drop = FALSE]; w <- w[keep] / sum(w[keep])

  # --- Mixture in the format of pi_extract_joint() --------------------------------
  to_absolute <- pi_absolute_mean_transform(R)
  spot_error <- outer(1 / stats$n, 1 / stats$n, "+"); diag(spot_error) <- 0
  theta_names <- c("Log precision for the Gaussian observations",
                   paste0("Log precision for gene_r", seq_len(R)))
  configurations <- lapply(seq_along(points), function(k) {
    p <- points[[k]]
    list(mu = as.numeric(to_absolute %*% p$s_mean),
         covariance = to_absolute %*% p$s_cov %*% t(to_absolute),
         tau2 = p$tau2, kappa2 = p$sigma_e2 * spot_error,
         theta = setNames(thetas[k, ], theta_names))
  })
  mix_mean <- Reduce("+", Map(function(cf, wk) cf$mu * wk, configurations, w))
  mix_second <- Reduce("+", Map(function(cf, wk) (cf$covariance + tcrossprod(cf$mu)) * wk,
                                configurations, w))
  joint <- list(
    configurations = configurations, weights = w, region_labels = labels,
    region_spots = stats$n, estimand = PI_ESTIMAND,
    baseline_region = labels[1],
    approximation = sprintf("closed-form conditional Gaussians on a %d^%d hyperparameter grid (+-%g sd)",
                            cfg$closed_form_grid_points, d, cfg$closed_form_grid_halfwidth),
    weight_source = "normalized exp(log posterior of theta) on the grid",
    max_reported_mean_skewness = 0,
    marginal_diagnostic = data.frame(region = labels, mixture_mean = mix_mean,
                                     INLA_marginal_mean = NA,
                                     mixture_sd = sqrt(pmax(0, diag(mix_second) - mix_mean^2)),
                                     INLA_marginal_sd = NA),
    extraction_seconds = 0)

  counts <- stats$n
  contrasts <- screen_region_contrasts(counts, cfg$other_weights)
  metadata <- list(
    N = N, K = K, R = R, mode = mode, engine = "closed_form",
    region_labels = labels, region_spots = counts,
    other_weights = contrasts$rest_weights[[1]], rest_weights = contrasts$rest_weights,
    gene_dependence = "conditional_iid",
    model_version = as.character(utils::packageVersion("PathwayLGM")),
    cfg = cfg, family = "gaussian",
    fit_seconds = proc.time()[["elapsed"]] - started,
    hyperparameters = d, hessian_floored = hessian_floored,
    tau_prior = list(internal_parameter = "log precision", prior = cfg$tau_prior,
                     mean = cfg$tau_logprecision_mean, sd = cfg$tau_logprecision_sd,
                     parameters = c(cfg$tau_logprecision_mean, 1 / cfg$tau_logprecision_sd^2)),
    grid = list(points_per_axis = cfg$closed_form_grid_points,
                halfwidth_sd = cfg$closed_form_grid_halfwidth,
                points_kept = length(w), theta_mode = setNames(opt$par, theta_names)),
    intercept_name = "intercept", regional_offset_names = paste0("mu_r", seq.int(2, R)),
    contrast_coefficients = contrasts$coefficients,
    gene_labels = colnames(Y), composition_diagnostics = NULL,
    offset_convention = "zero (SCT/normalized Gaussian response)")
  list(fit = NULL, joint = if (isTRUE(cfg$pi)) joint else NULL, metadata = metadata,
       closed_form = list(points = points, thetas = thetas, weights = w))
}


#' Summaries of a closed-form fit with the fields of screen_summarize()
#' @noRd
closed_form_summarize <- function(result) {
  started <- proc.time()[["elapsed"]]
  md <- result$metadata; cf <- result$closed_form; w <- cf$weights; R <- md$R; K <- md$K
  to_absolute <- pi_absolute_mean_transform(R)
  state_mean <- sapply(cf$points, `[[`, "s_mean")              # R x points
  summ_lin <- function(z) {   # linear combination z' s of the canonical states
    m <- as.numeric(crossprod(z, state_mean))
    s <- sqrt(vapply(cf$points, function(p) as.numeric(crossprod(z, p$s_cov %*% z)), numeric(1)))
    cf_mixture_summary(m, s, w)
  }
  e <- function(r) as.numeric(seq_len(R) == r)

  # fixed effects: baseline intercept (alpha0 + mean(u)), relative offsets, absolute means
  fixed <- rbind(
    data.frame(parameter = "baseline_intercept", t(summ_lin(e(1))[1:4])),
    do.call(rbind, lapply(seq.int(2, R), function(r)
      data.frame(parameter = paste0("relative_mu_r", r), t(summ_lin(e(r))[1:4])))),
    do.call(rbind, lapply(seq_len(R), function(r)
      data.frame(parameter = paste0("absolute_mu_r", r), t(summ_lin(to_absolute[r, ])[1:4])))))
  inverse_exp <- function(x) exp(-x)
  hyper <- rbind(
    do.call(rbind, lapply(seq_len(R), function(r)
      data.frame(parameter = paste0("tau2_r", r),
                 t(cf_hyper_summary(cf$thetas[, r + 1], w, inverse_exp, TRUE))))),
    data.frame(parameter = "sigma_e2", t(cf_hyper_summary(cf$thetas[, 1], w, inverse_exp, TRUE))))
  parameters <- rbind(fixed, hyper); rownames(parameters) <- NULL

  # regional contrasts (same coefficients as the INLA linear combinations)
  contrasts <- do.call(rbind, lapply(names(md$contrast_coefficients), function(name) {
    z_abs <- md$contrast_coefficients[[name]]              # coefficients on absolute means
    s <- summ_lin(as.numeric(crossprod(to_absolute, z_abs)))
    data.frame(quantity = name, mean = s[["mean"]], lower = s[["lower"]], median = s[["median"]],
               upper = s[["upper"]], p_gt_zero = s[["p_gt_zero"]],
               lfsr = min(s[["p_gt_zero"]], 1 - s[["p_gt_zero"]]))
  }))
  domains <- screen_domain_calls(contrasts, md$region_labels)

  gene_coefficients <- do.call(rbind, lapply(seq_len(K), function(g) {
    m <- vapply(cf$points, function(p) p$a_mean[g], numeric(1))
    s <- sqrt(vapply(cf$points, `[[`, numeric(1), "a_var"))
    data.frame(gene = if (is.null(md$gene_labels)) g else md$gene_labels[g],
               parameter = "a_g", t(cf_mixture_summary(m, s, w)[1:4]))
  }))

  pi_result <- if (!is.null(result$joint)) {
    pathway_pi_integrate(result$joint, alpha = md$cfg$pi_alpha,
                         report_p_value = md$cfg$pi_report_p_value,
                         nodes = md$cfg$pi_nodes, max_nodes = md$cfg$pi_max_nodes,
                         tolerance = md$cfg$pi_tolerance)
  }
  list(
    pi = pi_result, joint = result$joint, parameters = parameters,
    gene_covariate_coefficients = NULL, gene_coefficients = gene_coefficients,
    mean_contrasts = contrasts, domains = domains,
    mean_contrast = contrasts[contrasts$quantity == "dmu_1_vs_rest", , drop = FALSE],
    variance_contrast = data.frame(quantity = "log_tau2_1_over_weighted_within_rest",
                                   mean = NA, lower = NA, median = NA, upper = NA,
                                   p_gt_zero = NA, status = "not_computed_from_separate_marginals"),
    fixed_marginals = NULL, hyper_marginals = NULL, internal_hyper_marginals = NULL,
    contrast_marginals = NULL,
    theta_grid = data.frame(cf$thetas, weight = w),
    metadata = md, warnings = NULL,
    postprocess_seconds = proc.time()[["elapsed"]] - started)
}
