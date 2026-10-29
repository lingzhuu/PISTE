# =============================================================================
# joint_sections.R
#
# Closed-form joint fit of adjacent sections (observed mode). The regional
# comparison of a single section cannot separate a gene's regional effect from
# gene-specific spatial variation on the scale of a region. Sections of the
# same tissue block share the regional effects, whereas the spatial variation
# specific to each section is independent, so a joint fit of several sections
# splits the gene-level regional effect into a shared and a section-specific
# component and defines pi on the shared one.
#
# Model (S sections, K genes, R regions, region 1 = baseline, section 1 =
# reference section)
#
#   y[i, g, k] = alpha0 + delta_k + u_{g,k} + mu_r(i) + gamma_{g,r(i)}
#                + eta_{g,r(i),k} + e[i, g, k],      mu_1 = 0, delta_1 = 0
#
#   alpha0        ~ N(0, 100 (1 - 1/K))
#   delta_k, k>1  ~ N(0, 100)             section offset
#   mu_r, r > 1   ~ N(0, 100)             regional means shared by the sections
#   u_{g,k}       ~ N(0, 100) iid         gene intercept of each section
#   gamma_{g,r}   ~ N(0, tau2_r) iid      gene x region effect shared by the sections
#   eta_{g,r,k}   ~ N(0, nu2_r) iid       gene x region x section deviation
#   e[i, g, k]    ~ N(0, sigma_e2_k)
#   log(1 / tau2_r), log(1 / nu2_r) ~ N(tau_logprecision_mean, tau_logprecision_sd^2)
#   1 / sigma_e2_k ~ pc.prec(u = 2, alpha = .05)
#
# Hyperparameters theta = (log 1/sigma_e2_1..S, log 1/tau2_1..R, log 1/nu2_1..R).
#
# Given theta, integrating out u, gamma and eta leaves every gene's vector of
# the S x R section-by-region means Gaussian with the same covariance for every
# gene,
#
#   V[c, c'] = 100 [k = k'] + tau2_r [r = r'] + (nu2_r + sigma_e2_k / n_c) [c = c'],
#
# so the closed form of closed_form_screen.R carries over: the data enter
# through the per-gene cell means, the cell sizes and the within-cell sums of
# squares of each section, and theta is integrated on a regular grid in the
# standardised coordinates of its Laplace approximation. When a variance sits
# at the boundary (nu2_r -> 0) the Hessian can be flat or slightly indefinite;
# its eigenvalues are floored so that the grid spans that axis widely, and the
# prior keeps the posterior proper (metadata$hessian_floored records it).
#
# pi is defined on the shared component: tau2_r only, with the regional means
# of the reference section and kappa2_rs = sigma_e2_1 (1 / n_r1 + 1 / n_s1).
# The section-specific variances nu2_r are reported in the parameter table.
# Spot positions play no role, so the sections need not be aligned.
# =============================================================================


#' Sufficient statistics of every section: cell sizes, per-gene cell means,
#' within-cell sums of squares
#' @noRd
js_sufficient_statistics <- function(Y_list, region_index_list, R) {
  S <- length(Y_list); K <- ncol(Y_list[[1]])
  ybar <- NULL; n <- NULL; within <- numeric(S); N <- numeric(S)
  for (k in seq_len(S)) {
    nk <- tabulate(region_index_list[[k]], R)
    if (any(nk == 0)) stop("Every region needs at least one spot in every section.")
    sums <- rowsum(Y_list[[k]], region_index_list[[k]], reorder = TRUE)    # R x K
    ybar <- cbind(ybar, t(sums / nk))                                       # K x (S R)
    n <- c(n, nk)
    within[k] <- sum(Y_list[[k]]^2) - sum(sums^2 / nk)
    N[k] <- nrow(Y_list[[k]])
  }
  sec <- rep(seq_len(S), each = R); reg <- rep(seq_len(R), S)
  list(ybar = ybar, n = n, within = within, N = N, K = K, R = R, S = S,
       sec = sec, reg = reg, same_sec = outer(sec, sec, "=="),
       same_reg = outer(reg, reg, "=="))
}


#' Everything at one theta: log posterior and the conditional Gaussian of the
#' regional means of the reference section
#' @noRd
js_at_theta <- function(theta, stats, prior) {
  K <- stats$K; R <- stats$R; S <- stats$S; C <- length(stats$n)
  sigma_e2 <- exp(-theta[seq_len(S)])
  tau2 <- exp(-theta[S + seq_len(R)])
  nu2 <- exp(-theta[S + R + seq_len(R)])

  # beta = (alpha0, mu_2..mu_R, delta_2..delta_S) -> cell means X beta
  X <- cbind(1, outer(stats$reg, seq.int(2, R), "==") * 1,
             if (S > 1) outer(stats$sec, seq.int(2, S), "==") * 1)
  V <- prior$u_var * stats$same_sec + stats$same_reg * tau2[stats$reg] +
    diag(nu2[stats$reg] + sigma_e2[stats$sec] / stats$n, C)
  V_chol <- chol(V)
  V_inv <- chol2inv(V_chol)
  P_beta <- diag(c(1 / prior$alpha0_var, rep(1 / prior$mu_var, ncol(X) - 1)), ncol(X))
  XtVinv <- crossprod(X, V_inv)
  Q_chol <- chol(P_beta + K * XtVinv %*% X)
  b <- XtVinv %*% colSums(stats$ybar)
  m_beta <- backsolve(Q_chol, forwardsolve(t(Q_chol), b))
  C_beta <- chol2inv(Q_chol)

  quad <- sum((stats$ybar %*% V_inv) * stats$ybar)
  cells_per_section <- tabulate(stats$sec, S)
  log_lik <- -K * C / 2 * log(2 * pi) - K * sum(log(diag(V_chol))) - quad / 2 +
    sum(log(diag(P_beta))) / 2 - sum(log(diag(Q_chol))) + sum(b * m_beta) / 2 -
    sum(K * (stats$N - cells_per_section) / 2 * log(2 * pi * sigma_e2)) -
    K / 2 * sum(log(stats$n)) - sum(stats$within / (2 * sigma_e2))
  lambda <- -log(prior$pc_alpha) / prior$pc_u
  log_prior <- sum(log(lambda / 2) - lambda * exp(-theta[seq_len(S)] / 2) - theta[seq_len(S)] / 2) +
    sum(dnorm(theta[-seq_len(S)], prior$tau_mean, prior$tau_sd, log = TRUE))

  # absolute regional means of the reference section: alpha0 + mu_r
  L <- cbind(1, rbind(0, diag(R - 1)), matrix(0, R, max(S - 1, 0)))
  list(log_post = log_lik + log_prior,
       beta_mean = as.numeric(m_beta), beta_cov = C_beta,
       s_mean = as.numeric(L %*% m_beta), s_cov = L %*% C_beta %*% t(L),
       sigma_e2 = sigma_e2, tau2 = tau2, nu2 = nu2)
}


#' Fit the observed-mode model jointly to several sections
#'
#' Sections of the same tissue block share the regional means and the
#' gene-level regional effects, whereas gene-specific variation on the scale
#' of a region is independent between sections. The joint fit splits the
#' gene-level regional effect of every region into a component shared by the
#' sections, with variance tau2_r, and a component specific to each section,
#' with variance nu2_r, and defines pi on the shared component. Each section
#' has its own intercept, gene intercepts and error variance, so the sections
#' may be normalised separately and need not be aligned. The fit is in closed
#' form, like [closed_form_joint_fit()], with the S + 2R hyperparameters
#' integrated on a grid.
#'
#' @param Y_list List of spots x genes matrices of normalized expression, one
#'   per section, with the same genes in the same order.
#' @param region_list List of region labels per spot, one vector per section.
#'   Every region must be present in every section.
#' @param cfg Options as for [screen_joint_fit()]. The grid is set by
#'   `joint_grid_points` (points per hyperparameter axis),
#'   `closed_form_grid_halfwidth`, `closed_form_weight_threshold` and
#'   `closed_form_max_points`.
#' @return A list with `fit = NULL`, `joint` (format of the `joint` element
#'   of [screen_joint_fit()], defined on the shared component and the regional
#'   means of the first section), `metadata` (with
#'   `engine = "closed_form_sections"`, `sections`, `section_spots` and
#'   `hessian_floored`) and `closed_form` (the conditional Gaussians, grid
#'   points and weights). Pass it to [screen_summarize()].
#' @export
joint_sections_fit <- function(Y_list, region_list, cfg = list()) {
  started <- proc.time()[["elapsed"]]
  screen_check_legacy_cfg(cfg)
  screen_check_mode_cfg(cfg, "observed")
  screen_check_integration_cfg(cfg)
  cfg <- modifyList(screen_default_cfg("observed"), cfg)
  if (cfg$tau_prior != "normal") stop("The joint fit implements tau_prior = 'normal' only.")
  if (length(cfg$gene_specific_covariates)) stop("The joint fit has no gene-specific covariates.")
  stopifnot(is.list(Y_list), is.list(region_list), length(Y_list) == length(region_list),
            length(Y_list) >= 1, is.numeric(cfg$joint_grid_points), cfg$joint_grid_points >= 3)

  Y_list <- lapply(Y_list, as.matrix)
  S <- length(Y_list); K <- ncol(Y_list[[1]])
  genes <- colnames(Y_list[[1]])
  for (k in seq_len(S)) {
    stopifnot(ncol(Y_list[[k]]) == K, nrow(Y_list[[k]]) == length(region_list[[k]]),
              all(is.finite(Y_list[[k]])))
    if (!is.null(genes) && !identical(colnames(Y_list[[k]]), genes)) {
      stop("All sections need the same genes in the same order.")
    }
  }
  labels <- screen_region_labels(unlist(lapply(region_list, as.character)), cfg$baseline_region)
  R <- length(labels)
  stopifnot(R >= 2, K >= 2)
  region_index_list <- lapply(region_list, function(r) {
    idx <- match(as.character(r), labels)
    if (anyNA(idx)) stop("Region labels contain NA.")
    idx
  })
  d <- S + 2 * R
  n_points <- cfg$joint_grid_points^d
  if (n_points > cfg$closed_form_max_points) {
    stop("The hyperparameter grid would have ", n_points, " points (S = ", S, ", R = ", R,
         "); reduce joint_grid_points or raise closed_form_max_points.")
  }

  stats <- js_sufficient_statistics(Y_list, region_index_list, R)
  prior <- list(alpha0_var = 100 * (1 - 1 / K), mu_var = 100, u_var = 100,
                pc_u = 2, pc_alpha = .05,
                tau_mean = cfg$tau_logprecision_mean, tau_sd = cfg$tau_logprecision_sd)

  # --- Mode and curvature of the hyperparameter posterior ------------------------
  neg_log_post <- function(theta) {
    v <- tryCatch(-js_at_theta(theta, stats, prior)$log_post, error = function(e) Inf)
    if (is.finite(v)) v else 1e100
  }
  start <- c(-log(stats$within / (stats$N * K)), rep(cfg$tau_logprecision_mean, 2 * R))
  opt <- optim(start, neg_log_post, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 2000))
  eig <- eigen(optimHess(opt$par, neg_log_post), symmetric = TRUE)
  hessian_floored <- any(eig$values < cfg$closed_form_hessian_floor)
  eig$values <- pmax(eig$values, cfg$closed_form_hessian_floor)
  to_theta <- eig$vectors %*% diag(1 / sqrt(eig$values), d)

  # --- Grid in standardised coordinates ------------------------------------------
  axis <- seq(-cfg$closed_form_grid_halfwidth, cfg$closed_form_grid_halfwidth,
              length.out = cfg$joint_grid_points)
  Z <- as.matrix(expand.grid(rep(list(axis), d)))
  thetas <- t(opt$par + to_theta %*% t(Z))
  points <- lapply(seq_len(nrow(thetas)), function(j) js_at_theta(thetas[j, ], stats, prior))
  log_post <- vapply(points, `[[`, numeric(1), "log_post")
  w <- exp(log_post - max(log_post)); w <- w / sum(w)
  keep <- w >= cfg$closed_form_weight_threshold
  points <- points[keep]; thetas <- thetas[keep, , drop = FALSE]; w <- w[keep] / sum(w[keep])

  # --- Mixture in the format of pi_extract_joint(), on the shared component --------
  n1 <- stats$n[stats$sec == 1]
  spot_error <- outer(1 / n1, 1 / n1, "+"); diag(spot_error) <- 0
  theta_names <- c(paste0("Log precision for the Gaussian observations, section ", seq_len(S)),
                   paste0("Log precision for gene_r", seq_len(R)),
                   paste0("Log precision for gene_r", seq_len(R), " by section"))
  configurations <- lapply(seq_along(points), function(k) {
    p <- points[[k]]
    list(mu = p$s_mean, covariance = (p$s_cov + t(p$s_cov)) / 2,
         tau2 = p$tau2, kappa2 = p$sigma_e2[1] * spot_error,
         theta = setNames(thetas[k, ], theta_names))
  })
  mix_mean <- Reduce("+", Map(function(cf, wk) cf$mu * wk, configurations, w))
  mix_second <- Reduce("+", Map(function(cf, wk) (cf$covariance + tcrossprod(cf$mu)) * wk,
                                configurations, w))
  joint <- list(
    configurations = configurations, weights = w, region_labels = labels,
    region_spots = n1, estimand = PI_ESTIMAND, baseline_region = labels[1],
    approximation = sprintf("closed-form conditional Gaussians on a %d^%d hyperparameter grid (+-%g sd), %d sections",
                            cfg$joint_grid_points, d, cfg$closed_form_grid_halfwidth, S),
    weight_source = "normalized exp(log posterior of theta) on the grid",
    max_reported_mean_skewness = 0,
    marginal_diagnostic = data.frame(region = labels, mixture_mean = mix_mean,
                                     INLA_marginal_mean = NA,
                                     mixture_sd = sqrt(pmax(0, diag(mix_second) - mix_mean^2)),
                                     INLA_marginal_sd = NA),
    extraction_seconds = 0)

  contrasts <- screen_region_contrasts(n1, cfg$other_weights)
  section_spots <- matrix(stats$n, nrow = S, byrow = TRUE, dimnames = list(NULL, labels))
  metadata <- list(
    N = sum(stats$N), K = K, R = R, S = S, mode = "observed", engine = "closed_form_sections",
    region_labels = labels, region_spots = n1, sections = S, section_spots = section_spots,
    other_weights = contrasts$rest_weights[[1]], rest_weights = contrasts$rest_weights,
    gene_dependence = "conditional_iid",
    model_version = as.character(utils::packageVersion("PathwayLGM")),
    cfg = cfg, family = "gaussian",
    fit_seconds = proc.time()[["elapsed"]] - started,
    hyperparameters = d, hessian_floored = hessian_floored,
    tau_prior = list(internal_parameter = "log precision", prior = cfg$tau_prior,
                     mean = cfg$tau_logprecision_mean, sd = cfg$tau_logprecision_sd,
                     parameters = c(cfg$tau_logprecision_mean, 1 / cfg$tau_logprecision_sd^2)),
    grid = list(points_per_axis = cfg$joint_grid_points,
                halfwidth_sd = cfg$closed_form_grid_halfwidth,
                points_kept = length(w), theta_mode = setNames(opt$par, theta_names)),
    intercept_name = "intercept", regional_offset_names = paste0("mu_r", seq.int(2, R)),
    contrast_coefficients = contrasts$coefficients,
    gene_labels = genes, composition_diagnostics = NULL,
    offset_convention = "zero (SCT/normalized Gaussian response)")
  list(fit = NULL, joint = if (isTRUE(cfg$pi)) joint else NULL, metadata = metadata,
       closed_form = list(points = points, thetas = thetas, weights = w))
}


#' Summaries of a joint-section fit with the fields of screen_summarize()
#' @noRd
joint_sections_summarize <- function(result) {
  started <- proc.time()[["elapsed"]]
  md <- result$metadata; cf <- result$closed_form; w <- cf$weights
  R <- md$R; S <- md$S
  state_mean <- sapply(cf$points, `[[`, "s_mean")              # R x points
  summ_lin <- function(z) {   # linear combination z' s of the absolute regional means
    m <- as.numeric(crossprod(z, state_mean))
    s <- sqrt(vapply(cf$points, function(p) as.numeric(crossprod(z, p$s_cov %*% z)), numeric(1)))
    cf_mixture_summary(m, s, w)
  }
  e <- function(r) as.numeric(seq_len(R) == r)
  beta_summary <- function(j) {
    m <- vapply(cf$points, function(p) p$beta_mean[j], numeric(1))
    s <- sqrt(vapply(cf$points, function(p) p$beta_cov[j, j], numeric(1)))
    cf_mixture_summary(m, s, w)
  }

  fixed <- rbind(
    data.frame(parameter = "baseline_intercept", t(beta_summary(1)[1:4])),
    do.call(rbind, lapply(seq.int(2, R), function(r)
      data.frame(parameter = paste0("relative_mu_r", r), t(beta_summary(r)[1:4])))),
    do.call(rbind, lapply(seq_len(R), function(r)
      data.frame(parameter = paste0("absolute_mu_r", r), t(summ_lin(e(r))[1:4])))),
    if (S > 1) do.call(rbind, lapply(seq.int(2, S), function(k)
      data.frame(parameter = paste0("section_offset_", k), t(beta_summary(R + k - 1)[1:4])))))
  inverse_exp <- function(x) exp(-x)
  hyper <- rbind(
    do.call(rbind, lapply(seq_len(R), function(r)
      data.frame(parameter = paste0("tau2_r", r),
                 t(cf_hyper_summary(cf$thetas[, S + r], w, inverse_exp, TRUE))))),
    do.call(rbind, lapply(seq_len(R), function(r)
      data.frame(parameter = paste0("nu2_r", r),
                 t(cf_hyper_summary(cf$thetas[, S + R + r], w, inverse_exp, TRUE))))),
    do.call(rbind, lapply(seq_len(S), function(k)
      data.frame(parameter = paste0("sigma_e2_section_", k),
                 t(cf_hyper_summary(cf$thetas[, k], w, inverse_exp, TRUE))))))
  parameters <- rbind(fixed, hyper); rownames(parameters) <- NULL

  contrasts <- do.call(rbind, lapply(names(md$contrast_coefficients), function(name) {
    s <- summ_lin(md$contrast_coefficients[[name]])
    data.frame(quantity = name, mean = s[["mean"]], lower = s[["lower"]], median = s[["median"]],
               upper = s[["upper"]], p_gt_zero = s[["p_gt_zero"]],
               lfsr = min(s[["p_gt_zero"]], 1 - s[["p_gt_zero"]]))
  }))
  domains <- screen_domain_calls(contrasts, md$region_labels)

  pi_result <- if (!is.null(result$joint)) {
    pathway_pi_integrate(result$joint, alpha = md$cfg$pi_alpha,
                         report_p_value = md$cfg$pi_report_p_value,
                         nodes = md$cfg$pi_nodes, max_nodes = md$cfg$pi_max_nodes,
                         tolerance = md$cfg$pi_tolerance)
  }
  list(
    pi = pi_result, joint = result$joint, parameters = parameters,
    gene_covariate_coefficients = NULL, gene_coefficients = NULL,
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
