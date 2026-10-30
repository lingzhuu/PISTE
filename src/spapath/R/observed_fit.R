# Observed mode: Gaussian model of normalised expression, for one section or
# jointly for adjacent sections of the same tissue block.
#
#   y[i, g, k] = alpha0 + delta_k + u[g, k] + mu_r(i) + gamma[g, r(i)]
#                + eta[g, r(i), k] + e[i, g, k],          mu_1 = 0, delta_1 = 0
#
#   alpha0 ~ N(0, 100 (1 - 1/K)), delta_k ~ N(0, 100), mu_r ~ N(0, 100)
#   u[g, k]        ~ N(0, 100)       gene intercept of each section
#   gamma[g, r]    ~ N(0, tau2_r)    gene x region effect, shared by the sections
#   eta[g, r, k]   ~ N(0, nu2_r)     section-specific deviation (two or more sections)
#   e[i, g, k]     ~ N(0, sigma2_k)
#   log(1/tau2_r), log(1/nu2_r) ~ N(0, 3^2);  1/sigma2_k ~ pc.prec(2, 0.05)
#
# Given the variances theta, u, gamma and eta integrate out and every gene's
# vector of section x region means is Gaussian with the same covariance
#
#   V[c, c'] = 100 [k = k'] + tau2_r [r = r'] + (nu2_r + sigma2_k / n_c) [c = c'],
#
# so the posterior given theta is exact and needs only the cell means, cell
# sizes and within-cell sums of squares. theta is integrated on a grid around
# its posterior mode, as in INLA's grid strategy. pi uses the shared tau2_r,
# the regional means and the noise of a gene's regional difference combined
# over the sections, kappa2 = 1 / sum_k 1 / (sigma2_k (1/n_kr + 1/n_ks) + nu2_r + nu2_s).

# Grid: points per axis (one section / several sections), half width in
# posterior sd, smallest relative weight kept, and the smallest curvature used
# when a variance sits at zero and the Hessian is flat.
GRID <- list(points = 7, points_sections = 5, halfwidth = 3, min_weight = 1e-10,
             min_curvature = .05)


observed_statistics <- function(Y, region, section, R) {
  S <- max(section)
  cells <- (section - 1) * R + region
  n <- tabulate(cells, S * R)
  if (any(n == 0)) stop("Every region needs spots in every section.")
  sums <- rowsum(Y, cells, reorder = TRUE)
  list(ybar = t(sums / n), n = n, K = ncol(Y), R = R, S = S,
       N = tabulate(section, S),
       within = vapply(seq_len(S), function(k)
         sum(Y[section == k, ]^2) - sum(sums[(k - 1) * R + seq_len(R), ]^2 / n[(k - 1) * R + seq_len(R)]),
         numeric(1)),
       sec = rep(seq_len(S), each = R), reg = rep(seq_len(R), S))
}


observed_at_theta <- function(theta, st) {
  K <- st$K; R <- st$R; S <- st$S; C <- S * R
  sigma2 <- exp(-theta[seq_len(S)])
  tau2 <- exp(-theta[S + seq_len(R)])
  nu2 <- if (S > 1) exp(-theta[S + R + seq_len(R)]) else rep(0, R)

  # beta = (alpha0, mu_2..mu_R, delta_2..delta_S)
  X <- cbind(1, outer(st$reg, seq.int(2, R), "==") * 1,
             if (S > 1) outer(st$sec, seq.int(2, S), "==") * 1)
  V <- 100 * outer(st$sec, st$sec, "==") + outer(st$reg, st$reg, "==") * tau2[st$reg] +
    diag(nu2[st$reg] + sigma2[st$sec] / st$n, C)
  V_chol <- chol(V)
  V_inv <- chol2inv(V_chol)
  P <- diag(c(1 / (100 * (1 - 1 / K)), rep(1 / 100, ncol(X) - 1)), ncol(X))
  XtVinv <- crossprod(X, V_inv)
  Q_chol <- chol(P + K * XtVinv %*% X)
  b <- XtVinv %*% colSums(st$ybar)
  m <- backsolve(Q_chol, forwardsolve(t(Q_chol), b))

  log_lik <- -K * C / 2 * log(2 * pi) - K * sum(log(diag(V_chol))) -
    sum((st$ybar %*% V_inv) * st$ybar) / 2 +
    sum(log(diag(P))) / 2 - sum(log(diag(Q_chol))) + sum(b * m) / 2 -
    sum(K * (st$N - R) / 2 * log(2 * pi * sigma2)) - K / 2 * sum(log(st$n)) -
    sum(st$within / (2 * sigma2))
  lambda <- -log(.05) / 2
  log_prior <- sum(log(lambda / 2) - lambda * exp(-theta[seq_len(S)] / 2) - theta[seq_len(S)] / 2) +
    sum(dnorm(theta[-seq_len(S)], 0, 3, log = TRUE))

  # regional means of section 1, alpha0 + mu_r
  L <- cbind(1, rbind(0, diag(R - 1)), matrix(0, R, S - 1))
  C_beta <- chol2inv(Q_chol)
  list(log_post = log_lik + log_prior, beta_mean = as.numeric(m), beta_cov = C_beta,
       mean = as.numeric(L %*% m), cov = L %*% C_beta %*% t(L),
       sigma2 = sigma2, tau2 = tau2, nu2 = nu2)
}


#' Fit the observed-mode model of one gene set
#'
#' Fits the Gaussian model of normalised expression with a gene intercept,
#' regional means and gene-level regional effects. With `section`, adjacent
#' sections of the same tissue block are fitted jointly: the gene-level
#' regional effect is split into a part shared by the sections (variance
#' tau2_r) and a part specific to each section (variance nu2_r), and pi uses
#' the shared part. Each section keeps its own intercept, gene intercepts and
#' error variance.
#'
#' @param Y Spots x genes matrix of normalised expression.
#' @param region Region label of each spot.
#' @param section Optional section label of each spot; two or more sections
#'   are fitted jointly. Every region must be present in every section.
#' @param cfg Options as for [run_spapath()]: `pi`, `pi_alpha`,
#'   `pi_report_p_value`, `baseline_region`.
#' @return A fit to pass to [screen_summarize()].
#' @export
observed_fit <- function(Y, region, section = NULL, cfg = list()) {
  started <- proc.time()[["elapsed"]]
  cfg <- screen_merge_cfg(cfg, "observed")
  Y <- as.matrix(Y)
  labels <- screen_region_labels(region, cfg$baseline_region)
  if (is.null(section)) section <- rep(1, nrow(Y))
  sections <- unique(section)
  st <- observed_statistics(Y, match(region, labels), match(section, sections), length(labels))
  K <- st$K; R <- st$R; S <- st$S
  d <- S + R + if (S > 1) R else 0
  points_per_axis <- if (S > 1) GRID$points_sections else GRID$points
  if (points_per_axis^d > 1e6) stop("Too many regions or sections for the grid.")

  # Mode and curvature of the posterior of theta
  neg_log_post <- function(theta) {
    v <- tryCatch(-observed_at_theta(theta, st)$log_post, error = function(e) Inf)
    if (is.finite(v)) v else 1e100
  }
  start <- c(-log(st$within / (st$N * K)), rep(0, d - S))
  opt <- optim(start, neg_log_post, method = "BFGS", control = list(reltol = 1e-12, maxit = 2000))
  eig <- eigen(optimHess(opt$par, neg_log_post), symmetric = TRUE)
  curvature_floored <- any(eig$values < GRID$min_curvature)
  eig$values <- pmax(eig$values, GRID$min_curvature)

  # Grid in the standardised coordinates of the Laplace approximation
  axis <- seq(-GRID$halfwidth, GRID$halfwidth, length.out = points_per_axis)
  z <- as.matrix(expand.grid(rep(list(axis), d)))
  thetas <- t(opt$par + eig$vectors %*% diag(1 / sqrt(eig$values), d) %*% t(z))
  points <- lapply(seq_len(nrow(thetas)), function(j) observed_at_theta(thetas[j, ], st))
  log_post <- vapply(points, `[[`, numeric(1), "log_post")
  w <- exp(log_post - max(log_post))
  keep <- w / sum(w) >= GRID$min_weight
  points <- points[keep]; thetas <- thetas[keep, , drop = FALSE]; w <- w[keep] / sum(w[keep])

  # Mixture over the grid in the format read by pathway_pi_integrate()
  n_cell <- matrix(st$n, S, R, byrow = TRUE)
  n1 <- n_cell[1, ]
  combined_kappa2 <- function(p) {
    k2 <- 1 / Reduce(`+`, lapply(seq_len(S), function(k)
      1 / (p$sigma2[k] * outer(1 / n_cell[k, ], 1 / n_cell[k, ], "+") + outer(p$nu2, p$nu2, "+"))))
    diag(k2) <- 0
    k2
  }
  configurations <- lapply(points, function(p)
    list(mu = p$mean, covariance = (p$cov + t(p$cov)) / 2, tau2 = p$tau2,
         kappa2 = combined_kappa2(p)))
  joint <- list(configurations = configurations, weights = w, region_labels = labels,
                region_spots = n1, estimand = PI_ESTIMAND, baseline_region = labels[1],
                approximation = sprintf("exact given the variances, %d^%d grid", points_per_axis, d))

  contrasts <- screen_region_contrasts(n1, NULL)
  metadata <- list(
    N = sum(st$N), K = K, R = R, S = S, mode = "observed", engine = "observed",
    sections = as.character(sections), region_labels = labels, region_spots = n1,
    section_spots = matrix(st$n, S, R, byrow = TRUE, dimnames = list(sections, labels)),
    rest_weights = contrasts$rest_weights, contrast_coefficients = contrasts$coefficients,
    gene_labels = colnames(Y), gene_dependence = "conditional_iid", cfg = cfg,
    model_version = as.character(utils::packageVersion("PathwayLGM")),
    hyperparameters = d, curvature_floored = curvature_floored,
    fit_seconds = proc.time()[["elapsed"]] - started)
  list(fit = NULL, joint = if (isTRUE(cfg$pi)) joint else NULL, metadata = metadata,
       grid = list(points = points, thetas = thetas, weights = w))
}


# Mean and 2.5 / 50 / 97.5% quantiles of a mixture of normals, and P(> 0)
mixture_summary <- function(m, s, w) {
  cdf <- function(q) sum(w * pnorm((q - m) / s))
  q <- function(p) uniroot(function(x) cdf(x) - p, range(m - 8 * s, m + 8 * s), tol = 1e-10)$root
  c(mean = sum(w * m), lower = q(.025), median = q(.5), upper = q(.975), p_gt_zero = 1 - cdf(0))
}


# Summaries of observed_fit() with the fields of screen_summarize()
observed_summarize <- function(result) {
  started <- proc.time()[["elapsed"]]
  md <- result$metadata; g <- result$grid; w <- g$weights; R <- md$R; S <- md$S
  linear <- function(z) {   # z' (alpha0 + mu_r) of section 1
    mixture_summary(vapply(g$points, function(p) sum(z * p$mean), numeric(1)),
                    sqrt(vapply(g$points, function(p) as.numeric(crossprod(z, p$cov %*% z)), numeric(1))), w)
  }
  coefficient <- function(j) {
    mixture_summary(vapply(g$points, function(p) p$beta_mean[j], numeric(1)),
                    sqrt(vapply(g$points, function(p) p$beta_cov[j, j], numeric(1))), w)
  }
  variance <- function(j) {   # grid distribution of exp(-theta_j)
    o <- order(g$thetas[, j], decreasing = TRUE)
    x <- exp(-g$thetas[o, j]); cw <- cumsum(w[o])
    q <- function(p) x[min(which(cw >= p))]
    c(mean = sum(w * exp(-g$thetas[, j])), lower = q(.025), median = q(.5), upper = q(.975))
  }
  row <- function(name, v) data.frame(parameter = name, t(v[1:4]))
  e <- function(r) as.numeric(seq_len(R) == r)
  parameters <- rbind(
    row("baseline_intercept", coefficient(1)),
    do.call(rbind, lapply(seq.int(2, R), function(r) row(paste0("relative_mu_r", r), coefficient(r)))),
    do.call(rbind, lapply(seq_len(R), function(r) row(paste0("absolute_mu_r", r), linear(e(r))))),
    if (S > 1) do.call(rbind, lapply(seq.int(2, S), function(k)
      row(paste0("section_offset_", k), coefficient(R + k - 1)))),
    do.call(rbind, lapply(seq_len(R), function(r) row(paste0("tau2_r", r), variance(S + r)))),
    if (S > 1) do.call(rbind, lapply(seq_len(R), function(r) row(paste0("nu2_r", r), variance(S + R + r)))),
    do.call(rbind, lapply(seq_len(S), function(k)
      row(if (S > 1) paste0("sigma_e2_section_", k) else "sigma_e2", variance(k)))))
  rownames(parameters) <- NULL

  contrasts <- do.call(rbind, lapply(names(md$contrast_coefficients), function(name) {
    s <- linear(md$contrast_coefficients[[name]])
    data.frame(quantity = name, t(s), lfsr = min(s[["p_gt_zero"]], 1 - s[["p_gt_zero"]]))
  }))
  pi_result <- if (!is.null(result$joint)) {
    pathway_pi_integrate(result$joint, alpha = md$cfg$pi_alpha,
                         report_p_value = md$cfg$pi_report_p_value,
                         nodes = md$cfg$pi_nodes, max_nodes = md$cfg$pi_max_nodes,
                         tolerance = md$cfg$pi_tolerance)
  }
  list(pi = pi_result, joint = result$joint, parameters = parameters,
       mean_contrasts = contrasts, domains = screen_domain_calls(contrasts, md$region_labels),
       theta_grid = data.frame(g$thetas, weight = w), metadata = md,
       postprocess_seconds = proc.time()[["elapsed"]] - started)
}
