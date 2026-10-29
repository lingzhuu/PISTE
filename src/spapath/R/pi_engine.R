# =============================================================================
# pi_engine.R
#
# Directional probability engine for the conditional-iid pathway model.
#
# For a gene in region r its pathway state is Z_r ~ N(mu_r, tau2_r), where
# mu_r is the regional mean and tau2_r the between-gene variance in region r.
# For two regions r and s we report
#
#   score    E[ pi_rs ]          with  pi_rs = P(Z_r > Z_s | mu, tau2)
#                                     = Phi((mu_r - mu_s) /
#                                           sqrt(tau2_r + tau2_s + kappa2_rs))
#   evidence P( pi_rs > 0.5 )    i.e. the posterior probability that r is
#                                enriched relative to s
#   evidence_higher  P( pi_rs > 0.5 )    the same probability, and
#   evidence_lower   P( pi_rs < 0.5 )    its complement, reported as a pair
#
# A region pair is called "higher" or "lower" when the posterior probability
# of that direction is at least 1 - alpha / 2 (alpha = 0.05 by default), i.e.
# when the equal-tailed (1 - alpha) credible interval of pi_rs excludes 0.5.
# Optionally (report_p_value = TRUE) the table also carries
# p_value = 2 min(P(pi_rs > 0.5), P(pi_rs < 0.5)), the posterior probability
# of the less likely direction, doubled, so that the call is p_value <= alpha.
# It is a posterior probability, not a frequentist p-value: it says nothing
# about the frequency of false calls under repeated sampling.
#
# kappa2_rs is the measurement error of a gene's regional means. It keeps
# pi_rs near 0.5 when both the regional difference and tau2 are at the noise
# level.
#
#   Gaussian model           kappa2_rs = sigma_e2 * (1 / n_r + 1 / n_s),
#                            with n_r spots in region r.
#   Negative binomial model  each count carries the information
#                            w_ig = 1 / (1 / m_ig + 1 / size_g) on the log
#                            scale, where m_ig is its expected count without
#                            regional effects and size_g = nb_size0 * S_g.
#                            Gene g has kappa2_g = 1 / W_gr + 1 / W_gs, with
#                            W_gr the sum of w_ig over the spots of region r,
#                            and kappa2_rs is the harmonic mean of kappa2_g
#                            over genes. This weights genes by their
#                            information, as the model does when it estimates
#                            the regional means, so that nearly empty genes do
#                            not dominate. With equal noise in every spot and
#                            gene it reduces to the Gaussian formula.
#
# The expectation / probability is over the posterior of the regional
# means. That posterior is approximated by a mixture of multivariate Gaussians,
# one per INLA CCD hyperparameter configuration (see `pi_extract_joint()`).
#
# "Overall" results compare each region with a spot-count-weighted average of
# all other regions ("rest").
#
# Contents
#   pathway_pi()                  exported: conditional pi for given moments
#   pathway_pi_integrate()        exported: scores + evidence over the mixture
#   pi_extract_joint()            INLA fit -> low-dimensional Gaussian mixture
#   pi_* helpers                  internal building blocks
# =============================================================================


# Values of pi within this distance of 0.5 are snapped to exactly 0.5, so that
# ties (pi == 0.5) are detected reliably despite floating point round-off.
PI_TIE_TOLERANCE <- 32 * .Machine$double.eps

PI_ESTIMAND <- "Phi((mu_r - mu_s) / sqrt(tau2_r + tau2_s + kappa2_rs))"


# -----------------------------------------------------------------------------
# Conditional directional probability
# -----------------------------------------------------------------------------

#' P(Z_r > Z_s) for (possibly correlated) Gaussian states
#'
#' Vectorised over all arguments; each argument must be a scalar or share a
#' common length.
#'
#' @param mr,ms Means of Z_r and Z_s.
#' @param vr,vs Variances of Z_r and Z_s (must be > 0).
#' @param covariance Cov(Z_r, Z_s).
#' @param extra_variance Nonnegative variance added to Var(Z_r - Z_s), such as
#'   the measurement error kappa2_rs.
#' @return Probabilities in \[0, 1\].
#' @noRd
pi_gaussian_directional <- function(mr, ms, vr, vs, covariance = 0,
                                    extra_variance = 0) {
  sizes <- lengths(list(mr, ms, vr, vs, covariance, extra_variance))
  n <- max(sizes)
  if (n == 0 || any(!sizes %in% c(1, n))) {
    stop("Probability inputs must be scalar or have a common positive length.")
  }

  mr <- rep_len(mr, n)
  ms <- rep_len(ms, n)
  vr <- rep_len(vr, n)
  vs <- rep_len(vs, n)
  covariance <- rep_len(covariance, n)
  extra_variance <- rep_len(extra_variance, n)

  if (any(!is.finite(c(mr, ms, vr, vs, covariance, extra_variance))) ||
      any(vr <= 0 | vs <= 0) || any(extra_variance < 0)) {
    stop("Means must be finite and variances positive.")
  }
  if (any(abs(covariance) > sqrt(vr * vs) + 1e-8)) {
    stop("Invalid covariance.")
  }

  # Z_r - Z_s is Gaussian; a positive difference is directional enrichment.
  difference_variance <- vr + vs - 2 * covariance + extra_variance
  if (any(difference_variance <= 0)) stop("Degenerate difference variance.")

  prob <- stats::pnorm((mr - ms) / sqrt(difference_variance))
  prob[abs(prob - 0.5) <= PI_TIE_TOLERANCE] <- 0.5
  pmax(0, pmin(1, prob))
}


#' Conditional probability that the first pathway state exceeds the second
#'
#' Computes pi = Phi((mu_r - mu_s) / sqrt(tau2_r + tau2_s + kappa2)) for
#' independent Gaussian pathway states Z_r ~ N(mu_r, tau2_r) and
#' Z_s ~ N(mu_s, tau2_s), where kappa2 is the measurement error of a gene's
#' regional means.
#'
#' @param mu_r,mu_s Regional means.
#' @param tau2_r,tau2_s Regional (between-gene) variances.
#' @param kappa2 Measurement-error variance of a gene's regional means, e.g.
#'   sigma_e2 * (1 / n_r + 1 / n_s) in the Gaussian model; 0 gives
#'   P(Z_r > Z_s).
#' @return Numeric vector of probabilities in \[0, 1\].
#' @export
pathway_pi <- function(mu_r, mu_s, tau2_r, tau2_s, kappa2 = 0) {
  pi_gaussian_directional(mu_r, mu_s, tau2_r, tau2_s, extra_variance = kappa2)
}


# -----------------------------------------------------------------------------
# Extracting the low-dimensional joint posterior from an INLA fit
# -----------------------------------------------------------------------------

#' Covariance of linear combinations of the latent field
#'
#' Returns Cov(t(W) x) = t(W) Q^{-1} W for a sparse weight matrix W, using
#' sparse Cholesky solves instead of inverting the full precision matrix.
#' Accepts INLA's full, upper- or lower-triangular storage of a symmetric
#' matrix.
#'
#' @param Q Precision matrix (dense or sparse).
#' @param weights Sparse matrix with one column per linear combination.
#' @return Dense symmetric covariance matrix.
#' @noRd
pi_selected_covariance <- function(Q, weights) {
  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)

  if (isSymmetric(Q)) {
    Q <- Matrix::forceSymmetric(Q)
  } else if (Matrix::nnzero(Matrix::tril(Q, -1)) == 0) {
    Q <- Matrix::forceSymmetric(Q, uplo = "U")
  } else if (Matrix::nnzero(Matrix::triu(Q, 1)) == 0) {
    Q <- Matrix::forceSymmetric(Q, uplo = "L")
  } else {
    stop("Unsupported nonsymmetric INLA precision storage.")
  }

  solution <- Matrix::solve(Matrix::Cholesky(Q, LDL = FALSE), weights,
                            system = "A")
  covariance <- as.matrix(Matrix::crossprod(weights, solution))

  # Symmetrise to remove round-off asymmetry.
  (covariance + t(covariance)) / 2
}


#' Sparse map from the latent field to the canonical regional states
#'
#' Column 1 is the baseline intercept and columns 2..R the relative regional
#' offsets. In adjusted mode the gene intercepts are an intercept plus K free
#' auxiliary effects u_g, and the canonical baseline intercept is
#' intercept + mean(u), so column 1 also puts weight 1/K on every u_g.
#'
#' @param configs `fit$misc$configs` from an INLA fit.
#' @param n_latent Length of each configuration's latent mean vector.
#' @param effect_names Names of the intercept and the relative offsets.
#' @param gene_intercept_name Name of the auxiliary gene-intercept effect, or
#'   NULL in observed mode.
#' @return List with the sparse `n_latent` x R matrix `map` and `state_rows`,
#'   the latent positions of the intercept and the offsets.
#' @noRd
pi_regional_state_map <- function(configs, n_latent, effect_names,
                                  gene_intercept_name = NULL) {
  contents <- configs$contents
  starts <- contents$start

  if (isTRUE(configs$.preopt)) {
    # With .preopt the saved Q and mean exclude the (leading) predictor
    # blocks, so shift the start offsets accordingly.
    is_predictor <- contents$tag %in% c("APredictor", "Predictor")
    if (any(is_predictor) &&
        !identical(which(is_predictor), seq_len(sum(is_predictor)))) {
      stop("Unrecognized INLA predictor layout.")
    }
    starts <- starts - sum(contents$length[is_predictor])
  }

  selected <- match(effect_names, contents$tag)
  if (anyNA(selected) || any(contents$length[selected] != 1)) {
    stop("Cannot identify intercept and regional offsets in INLA configuration.")
  }
  rows <- starts[selected]
  cols <- seq_along(effect_names)
  values <- rep(1, length(effect_names))

  if (!is.null(gene_intercept_name)) {
    block <- match(gene_intercept_name, contents$tag)
    if (is.na(block)) stop("Cannot identify the gene intercepts in INLA configuration.")
    n_genes <- contents$length[block]
    rows <- c(rows, starts[block] + seq_len(n_genes) - 1)
    cols <- c(cols, rep(1, n_genes))
    values <- c(values, rep(1 / n_genes, n_genes))
  }

  if (any(rows < 1 | rows > n_latent)) {
    stop("INLA compact precision and mean layout do not agree.")
  }
  list(map = Matrix::sparseMatrix(i = rows, j = cols, x = values,
                                  dims = c(n_latent, length(effect_names))),
       state_rows = starts[selected])
}


#' Normalised CCD configuration weights
#'
#' Weights are proportional to `exp(log.posterior)`, as in INLA's own
#' configuration mixture.
#' @noRd
pi_configuration_weights <- function(configs) {
  log_weights <- vapply(configs$config, function(config) config$log.posterior,
                        numeric(1))
  if (any(!is.finite(log_weights))) stop("Nonfinite CCD configuration weights.")

  weights <- exp(log_weights - max(log_weights))
  weights / sum(weights)
}


#' Matrix mapping (intercept, offset_2, ..., offset_R) to absolute means
#'
#' Region 1 is the baseline: mu_1 = intercept, mu_r = intercept + offset_r.
#' @noRd
pi_absolute_mean_transform <- function(n_regions) {
  transform <- matrix(0, n_regions, n_regions)
  transform[, 1] <- 1
  for (r in seq.int(2, n_regions)) transform[r, r] <- 1
  transform
}


#' Negative binomial measurement error of a gene's regional means
#'
#' Each count carries the log-scale information w = 1 / (1 / m + 1 / size).
#' For every gene and region pair (r, s) the measurement error is
#' kappa2_g = 1 / W_gr + 1 / W_gs, where W_gr sums w over the spots of region
#' r; the pathway value is the harmonic mean of kappa2_g over genes.
#'
#' @param expected N x K matrix of expected counts m without regional effects.
#' @param size Length-K NB sizes, nb_size0 * S_g.
#' @param region_index Length-N region index in 1..n_regions.
#' @param n_regions Number of regions.
#' @return n_regions x n_regions matrix with a zero diagonal.
#' @noRd
pi_nb_kappa2 <- function(expected, size, region_index, n_regions) {
  information <- expected * rep(size, each = nrow(expected)) /
    (expected + rep(size, each = nrow(expected)))
  region_information <- rowsum(information, region_index, reorder = TRUE)
  if (nrow(region_information) != n_regions) {
    stop("Every region needs at least one spot.")
  }

  kappa2 <- matrix(0, n_regions, n_regions)
  for (r in seq_len(n_regions - 1)) {
    for (s in seq.int(r + 1, n_regions)) {
      # 1 / kappa2_g = W_gr W_gs / (W_gr + W_gs); zero for a gene without counts.
      total <- region_information[r, ] + region_information[s, ]
      gene_information <- ifelse(total > 0,
                                 region_information[r, ] *
                                   region_information[s, ] / total, 0)
      if (mean(gene_information) <= 0) {
        stop("No gene carries counts in both regions.")
      }
      kappa2[r, s] <- kappa2[s, r] <- 1 / mean(gene_information)
    }
  }
  kappa2
}


#' Extract the CCD mixture of regional means from an INLA fit
#'
#' For every saved CCD hyperparameter configuration, takes the conditional
#' Gaussian of (baseline intercept, relative regional offsets), maps it to the
#' absolute regional means and records the regional variances tau2_r and the
#' measurement-error variances kappa2_rs. The fit
#' must have been run with `control.compute = list(config = TRUE)` and an
#' unconstrained latent model.
#'
#' @param fit INLA fit object.
#' @param region_labels Region labels, baseline first.
#' @param region_spots Number of spots per region.
#' @param intercept_name Name of the baseline intercept latent effect.
#' @param regional_offset_names Names of the relative offsets mu_r2..mu_rR.
#' @param gene_intercept_name Name of the auxiliary gene-intercept effect in
#'   adjusted mode, NULL in observed mode.
#' @param nb_counts Adjusted mode: list with `expected` (N x K expected counts
#'   without regional effects), `dispersion_scale` (S_g) and `region_index`,
#'   used for the negative binomial kappa2_rs. NULL in observed mode.
#' @return A list with `configurations` (each with `mu`, `covariance`, `tau2`,
#'   `kappa2`, `theta`), normalised `weights`, region information and a
#'   `marginal_diagnostic` table comparing mixture moments with INLA's own
#'   marginals.
#' @noRd
pi_extract_joint <- function(fit, region_labels, region_spots,
                             intercept_name, regional_offset_names,
                             gene_intercept_name = NULL, nb_counts = NULL) {
  started <- proc.time()[["elapsed"]]

  configs <- fit$misc$configs
  if (is.null(configs) || length(configs$config) < 2) {
    stop("Pi inference requires multiple saved CCD configurations.")
  }
  if (!is.null(configs$constr) && length(configs$constr$A) > 0) {
    stop("Pi extraction currently requires the package's unconstrained latent model.")
  }

  n_regions <- length(region_labels)
  states <- pi_regional_state_map(configs, nrow(configs$config[[1]]$Q),
                                  c(intercept_name, regional_offset_names),
                                  gene_intercept_name)
  to_absolute <- pi_absolute_mean_transform(n_regions)
  weights <- pi_configuration_weights(configs)
  precision_tags <- paste0("Log precision for gene_r", seq_len(n_regions))
  inverse_spots <- 1 / region_spots
  spot_error <- outer(inverse_spots, inverse_spots, "+")
  diag(spot_error) <- 0

  mixture <- vector("list", length(weights))
  max_skewness <- 0

  for (k in seq_along(weights)) {
    # Release the previous configuration's sparse factorisation first.
    gc(FALSE)
    config <- configs$config[[k]]
    if (length(config$improved.mean) != nrow(config$Q) ||
        nrow(states$map) != nrow(config$Q)) {
      stop("INLA compact precision and mean layout do not agree.")
    }

    # Regional variances: tau2_r = 1 / precision of gene_r.
    precision_index <- match(precision_tags, names(config$theta))
    if (anyNA(precision_index)) {
      stop("Cannot identify regional log precisions by name.")
    }
    tau2 <- exp(-config$theta[precision_index])

    # Measurement error of a gene's regional means.
    if (is.null(nb_counts)) {
      log_noise_precision <- config$theta["Log precision for the Gaussian observations"]
      if (is.na(log_noise_precision)) {
        stop("Cannot identify the Gaussian observation precision.")
      }
      kappa2 <- exp(-log_noise_precision) * spot_error
    } else {
      size_index <- grep("size for the nbinomial observations",
                         names(config$theta), fixed = TRUE)
      if (length(size_index) != 1) {
        stop("Cannot identify the negative binomial size.")
      }
      kappa2 <- pi_nb_kappa2(nb_counts$expected,
                             exp(config$theta[[size_index]]) *
                               nb_counts$dispersion_scale,
                             nb_counts$region_index, n_regions)
    }

    covariance <- pi_selected_covariance(config$Q, states$map)
    if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
      stop("Invalid selected covariance.")
    }

    # Transform (baseline intercept, relative offsets) to absolute regional
    # states only for the directional probability integration.
    mixture[[k]] <- list(
      mu = as.numeric(to_absolute %*%
                        as.numeric(Matrix::crossprod(states$map, config$improved.mean))),
      covariance = to_absolute %*% covariance %*% t(to_absolute),
      tau2 = as.numeric(tau2),
      kappa2 = kappa2,
      theta = config$theta
    )

    skewness <- config$skewness[states$state_rows]
    if (any(is.finite(skewness))) {
      max_skewness <- max(max_skewness, abs(skewness[is.finite(skewness)]))
    }
  }

  list(
    configurations = mixture,
    weights = weights,
    region_labels = region_labels,
    region_spots = region_spots,
    estimand = PI_ESTIMAND,
    baseline_region = region_labels[1],
    approximation = "CCD mixture of improved-mean conditional multivariate Gaussians; no skew transform",
    weight_source = "normalized exp(config$log.posterior), as used by INLA's configuration mixture",
    max_reported_mean_skewness = max_skewness,
    marginal_diagnostic = pi_mixture_diagnostic(fit, mixture, weights, region_labels),
    extraction_seconds = proc.time()[["elapsed"]] - started
  )
}


#' Compare the mixture's mean and sd with INLA's marginals of absolute_mu_r
#' @noRd
pi_mixture_diagnostic <- function(fit, mixture, weights, region_labels) {
  n_regions <- length(region_labels)

  mixture_mean <- Reduce("+", Map(function(component, w) component$mu * w,
                                  mixture, weights))
  mixture_second_moment <- Reduce("+", Map(function(component, w)
    (component$covariance + tcrossprod(component$mu)) * w, mixture, weights))
  mixture_sd <- sqrt(pmax(0, diag(mixture_second_moment) - mixture_mean^2))

  reference <- fit$summary.lincomb.derived[
    paste0("absolute_mu_r", seq_len(n_regions)), , drop = FALSE]

  data.frame(
    region = region_labels,
    mixture_mean = mixture_mean,
    INLA_marginal_mean = reference$mean,
    mixture_sd = mixture_sd,
    INLA_marginal_sd = reference$sd
  )
}


# -----------------------------------------------------------------------------
# Integrating scores and evidence over the mixture
# -----------------------------------------------------------------------------

#' Integrate pi scores and evidence from a low-dimensional CCD Gaussian mixture
#'
#' For every ordered region pair (r, s) computes
#' * **score**: the posterior mean of pi_rs = P(Z_r > Z_s), obtained
#'   analytically for each mixture component;
#' * **evidence**: the posterior probability that pi_rs > 0.5;
#' * **evidence_higher**, **evidence_lower**: P(pi_rs > 0.5) and
#'   P(pi_rs < 0.5), the posterior probabilities of each direction;
#' * **call**: `"higher"` or `"lower"` when the posterior probability of that
#'   direction is at least 1 - `alpha` / 2, otherwise `"none"`;
#' * **p_value** (only with `report_p_value = TRUE`): 2 min(P(pi_rs > 0.5),
#'   P(pi_rs < 0.5)), the doubled posterior probability of the less likely
#'   direction, so that the call is `p_value <= alpha`. It is a posterior
#'   probability, not a frequentist p-value.
#'
#' The "overall" table compares each region with a spot-count-weighted average
#' of all other regions. With two regions the evidence is computed
#' analytically. Otherwise it uses antithetic, nested, unrandomised Sobol
#' quadrature whose node count is doubled until two consecutive refinements are
#' stable within `tolerance`. The quadrature is deterministic and never touches
#' the RNG state.
#'
#' @param joint The `joint` element returned by [screen_joint_fit()].
#' @param alpha Level of the call, in (0, 1): a direction is called when its
#'   posterior probability is at least 1 - alpha / 2.
#' @param report_p_value Add the `p_value` column.
#' @param nodes Initial number of quadrature nodes (power of two, >= 16).
#' @param max_nodes Maximum number of nodes (power of two, >= 2 * nodes).
#' @param tolerance Convergence tolerance for refinement differences.
#' @return A list with data frames `pairwise` (one row per ordered pair) and
#'   `overall` (one row per region), and a list of `diagnostics`.
#' @export
pathway_pi_integrate <- function(joint, alpha = 0.05, report_p_value = FALSE,
                                 nodes = 4096, max_nodes = 65536,
                                 tolerance = 0.002) {
  started <- proc.time()[["elapsed"]]

  pi_validate_mixture(joint, nodes, max_nodes, tolerance)
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Use 0 < alpha < 1.")
  }

  n_regions <- length(joint$region_labels)
  weights <- joint$weights
  mixture <- joint$configurations
  counts <- joint$region_spots
  rest_weights <- pi_rest_weights(counts)
  pairs <- t(combn(seq_len(n_regions), 2))

  # Analytic scores; also returns a square root of each component covariance
  # for the quadrature below.
  analytic <- pi_analytic_scores(mixture, weights, pairs, n_regions)
  score <- analytic$score
  rest_score <- rowSums(replace(score, is.na(score), 0) * rest_weights)

  if (n_regions == 2) {
    evidence <- pi_two_region_evidence(mixture, weights)

    tables <- pi_result_tables(
      joint$region_labels, counts,
      score = score, evidence = evidence$pairwise,
      boundary = evidence$boundary_mass,
      rest_score = rest_score, rest_evidence = evidence$overall,
      rest_boundary = evidence$boundary_mass, alpha = alpha,
      report_p_value = report_p_value
    )
    diagnostics <- list(
      converged = TRUE, nodes = 0,
      integration = "analytic two-region Gaussian-mixture calculation",
      refinement = data.frame(nodes = 0,
                              evidence_refinement_difference = 0,
                              analytic_score_check = 0),
      configurations = length(weights),
      approximation = joint$approximation,
      estimand = PI_ESTIMAND, boundary = 0.5, alpha = alpha,
      strict_inequality = TRUE, deterministic = TRUE,
      posterior_sampling = FALSE,
      seconds = proc.time()[["elapsed"]] - started
    )
    return(c(tables, list(diagnostics = diagnostics)))
  }

  evidence <- pi_quadrature_evidence(mixture, weights, analytic$roots,
                                     score, pairs, rest_weights,
                                     nodes, max_nodes, tolerance)
  if (!evidence$converged) {
    warning("Pi quadrature reached max_nodes without two stable refinements.")
  }

  tables <- pi_result_tables(
    joint$region_labels, counts,
    score = score, evidence = evidence$pairwise,
    boundary = evidence$pairwise_tie,
    rest_score = rest_score, rest_evidence = evidence$overall,
    rest_boundary = evidence$overall_tie, alpha = alpha,
    report_p_value = report_p_value
  )
  diagnostics <- list(
    converged = evidence$converged, nodes = evidence$nodes,
    tolerance = tolerance,
    refinement = evidence$refinement,
    configurations = length(weights),
    approximation = joint$approximation,
    estimand = PI_ESTIMAND, boundary = 0.5, alpha = alpha,
    strict_inequality = TRUE,
    boundary_roundoff_tolerance = PI_TIE_TOLERANCE,
    deterministic = TRUE, posterior_sampling = FALSE,
    numerical_difference_is_not_error_bound = TRUE,
    seconds = proc.time()[["elapsed"]] - started
  )
  c(tables, list(diagnostics = diagnostics))
}


#' Validate the mixture and quadrature settings passed to pathway_pi_integrate()
#' @noRd
pi_validate_mixture <- function(joint, nodes, max_nodes, tolerance) {
  n_regions <- length(joint$region_labels)
  weights <- joint$weights

  if (n_regions < 2 ||
      length(weights) != length(joint$configurations) ||
      any(!is.finite(weights) | weights < 0) ||
      abs(sum(weights) - 1) > 1e-8 ||
      length(joint$region_spots) != n_regions ||
      any(!is.finite(joint$region_spots) | joint$region_spots <= 0)) {
    stop("Invalid low-dimensional mixture.")
  }

  # Checked left to right, so log2() is only reached for nodes >= 16.
  if (nodes < 16 || max_nodes < 2 * nodes ||
      any(log2(c(nodes, max_nodes)) != round(log2(c(nodes, max_nodes)))) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("Use powers of two, nodes>=16, max_nodes>=2*nodes, and positive tolerance.")
  }
  invisible(TRUE)
}


#' Row r: spot-count weights of all regions other than r (rows sum to one)
#' @noRd
pi_rest_weights <- function(counts) {
  n_regions <- length(counts)
  rest_weights <- matrix(rep(counts, each = n_regions), n_regions, n_regions)
  diag(rest_weights) <- 0
  rest_weights / rowSums(rest_weights)
}


#' Check one mixture component and return a square root of its covariance
#'
#' The returned matrix L satisfies t(L) %*% L = covariance, so that
#' `normal %*% L` turns standard normal rows into draws with that covariance.
#' @noRd
pi_component_root <- function(component, n_regions) {
  covariance <- component$covariance
  if (length(component$mu) != n_regions ||
      length(component$tau2) != n_regions ||
      !identical(dim(covariance), c(n_regions, n_regions)) ||
      !identical(dim(component$kappa2), c(n_regions, n_regions)) ||
      any(!is.finite(c(component$mu, component$tau2, covariance,
                       component$kappa2))) ||
      any(component$tau2 <= 0) || any(component$kappa2 < 0) ||
      max(abs(covariance - t(covariance))) > 1e-8) {
    stop("Invalid conditional moments.")
  }

  eig <- eigen(covariance, symmetric = TRUE)
  if (min(eig$values) < -1e-8 * max(1, max(eig$values))) {
    stop("Conditional covariance is not positive semidefinite.")
  }
  diag(sqrt(pmax(eig$values, 0)), n_regions) %*% t(eig$vectors)
}


#' Posterior mean of pi_rs for all ordered pairs, computed analytically
#'
#' Given a component's regional means mu ~ N(m, C), E[pi_rs] equals
#' Phi((m_r - m_s) / sqrt(tau2_r + tau2_s + kappa2_rs + C_rr + C_ss - 2 C_rs)),
#' a closed form per component.
#'
#' @return List with `score` (R x R, NA diagonal) and `roots` (covariance
#'   square roots, one per component).
#' @noRd
pi_analytic_scores <- function(mixture, weights, pairs, n_regions) {
  score <- matrix(0, n_regions, n_regions)
  roots <- vector("list", length(mixture))

  for (k in seq_along(mixture)) {
    component <- mixture[[k]]
    roots[[k]] <- pi_component_root(component, n_regions)
    mu <- component$mu
    tau2 <- component$tau2
    covariance <- component$covariance

    for (j in seq_len(nrow(pairs))) {
      r <- pairs[j, 1]
      s <- pairs[j, 2]
      p <- pi_gaussian_directional(mu[r], mu[s],
                                   tau2[r] + covariance[r, r],
                                   tau2[s] + covariance[s, s],
                                   covariance[r, s],
                                   component$kappa2[r, s])
      score[r, s] <- score[r, s] + weights[k] * p
      score[s, r] <- score[s, r] + weights[k] * (1 - p)
    }
  }

  diag(score) <- NA
  list(score = score, roots = roots)
}


#' Evidence for exactly two regions, computed analytically
#'
#' With two regions, pi_12 > 0.5 iff mu_1 > mu_2. Integrate this event
#' analytically over each CCD conditional Gaussian and mix the weights.
#'
#' @return List with `pairwise` (2 x 2 evidence matrix, NA diagonal),
#'   `overall` (evidence per region) and `boundary_mass` (P(mu_1 == mu_2)).
#' @noRd
pi_two_region_evidence <- function(mixture, weights) {
  per_component <- lapply(mixture, function(component) {
    covariance <- component$covariance
    difference_variance <- covariance[1, 1] + covariance[2, 2] -
      2 * covariance[1, 2]
    if (!is.finite(difference_variance) || difference_variance < 0) {
      stop("Invalid conditional mean-difference variance.")
    }

    difference <- component$mu[1] - component$mu[2]
    if (difference_variance == 0) {
      return(c(first = as.numeric(difference > 0),
               second = as.numeric(difference < 0),
               tie = as.numeric(difference == 0)))
    }
    p <- stats::pnorm(difference / sqrt(difference_variance))
    c(first = p, second = 1 - p, tie = 0)
  })

  mass <- colSums(do.call(rbind, per_component) * weights)
  evidence_12 <- mass[["first"]]
  evidence_21 <- mass[["second"]]

  list(
    pairwise = matrix(c(NA, evidence_21,
                        evidence_12, NA), 2, 2),
    overall = c(evidence_12, evidence_21),
    boundary_mass = mass[["tie"]]
  )
}


#' Antithetic standard-normal quadrature nodes from an unrandomised Sobol set
#'
#' Returns an n x d matrix whose rows come in pairs (x, -x). The first n / 2
#' rows are therefore a valid (coarser) rule on their own, which is used to
#' measure the refinement error. Never changes the RNG state.
#' @noRd
pi_antithetic_sobol_normals <- function(n, d) {
  half_n <- n / 2
  u <- qrng::sobol(half_n, d = d, randomize = "none", skip = 1)
  normal <- matrix(stats::qnorm(u), nrow = half_n, ncol = d)
  normal <- normal[rep(seq_len(half_n), each = 2), , drop = FALSE]
  normal[seq(2, n, 2), ] <- -normal[seq(2, n, 2), , drop = FALSE]
  normal
}


#' Evidence for three or more regions by adaptive quasi-Monte Carlo
#'
#' Doubles the number of nodes until the difference between the full rule and
#' its first half (for both the evidence and a check of the analytic scores) is
#' within `tolerance` for two consecutive refinements, or `max_nodes` is
#' reached.
#'
#' @return List with pairwise and overall evidence and tie mass, convergence
#'   status, the final node count and the refinement history.
#' @noRd
pi_quadrature_evidence <- function(mixture, weights, roots, score, pairs,
                                   rest_weights, nodes, max_nodes, tolerance) {
  n_regions <- nrow(score)
  history <- list()
  stable_refinements <- 0
  n <- nodes

  repeat {
    normal <- pi_antithetic_sobol_normals(n, n_regions)
    first_half <- seq_len(n / 2)

    # "coarse" uses only the first half of the nodes.
    evidence <- coarse <- tie <- matrix(0, n_regions, n_regions)
    rest_evidence <- rest_coarse <- rest_tie <- numeric(n_regions)
    quadrature_score <- matrix(0, n_regions, n_regions)

    for (k in seq_along(mixture)) {
      component <- mixture[[k]]
      w <- weights[k]
      tau2 <- component$tau2

      # Quadrature draws of the absolute regional means (n x R).
      mu <- sweep(normal %*% roots[[k]], 2, component$mu, "+")
      rest_pi <- matrix(0, n, n_regions)

      for (j in seq_len(nrow(pairs))) {
        r <- pairs[j, 1]
        s <- pairs[j, 2]
        p <- pathway_pi(mu[, r], mu[, s], tau2[r], tau2[s],
                        component$kappa2[r, s])

        evidence[r, s] <- evidence[r, s] + w * mean(p > 0.5)
        evidence[s, r] <- evidence[s, r] + w * mean(p < 0.5)
        coarse[r, s] <- coarse[r, s] + w * mean(p[first_half] > 0.5)
        coarse[s, r] <- coarse[s, r] + w * mean(p[first_half] < 0.5)
        tie[r, s] <- tie[r, s] + w * mean(p == 0.5)
        tie[s, r] <- tie[r, s]
        quadrature_score[r, s] <- quadrature_score[r, s] + w * mean(p)
        quadrature_score[s, r] <- quadrature_score[s, r] + w * mean(1 - p)

        # Weighted average of pi against every other region.
        rest_pi[, r] <- rest_pi[, r] + rest_weights[r, s] * p
        rest_pi[, s] <- rest_pi[, s] + rest_weights[s, r] * (1 - p)
      }

      rest_pi[abs(rest_pi - 0.5) <= PI_TIE_TOLERANCE] <- 0.5
      rest_evidence <- rest_evidence + w * colMeans(rest_pi > 0.5)
      rest_coarse <- rest_coarse +
        w * colMeans(rest_pi[first_half, , drop = FALSE] > 0.5)
      rest_tie <- rest_tie + w * colMeans(rest_pi == 0.5)
    }

    refinement_difference <- max(abs(evidence - coarse),
                                 abs(rest_evidence - rest_coarse))
    score_check <- max(abs(quadrature_score - score), na.rm = TRUE)
    history[[length(history) + 1]] <- data.frame(
      nodes = n,
      evidence_refinement_difference = refinement_difference,
      analytic_score_check = score_check
    )

    is_stable <- refinement_difference <= tolerance && score_check <= tolerance
    stable_refinements <- if (is_stable) stable_refinements + 1 else 0
    if (stable_refinements >= 2 || n >= max_nodes) break
    n <- n * 2
  }

  list(
    pairwise = evidence, pairwise_tie = tie,
    overall = rest_evidence, overall_tie = rest_tie,
    converged = stable_refinements >= 2,
    nodes = n,
    refinement = do.call(rbind, history)
  )
}


#' Assemble the pairwise and overall result tables
#'
#' @param evidence R x R matrix of P(pi_rs > 0.5); P(pi_rs < 0.5) is its
#'   transpose.
#' @param boundary Scalar or R x R matrix of tie mass for each ordered pair.
#' @param rest_evidence P(pi > 0.5) of each region against the rest.
#' @param alpha Level of the call.
#' @param report_p_value Add the `p_value` column.
#' @noRd
pi_result_tables <- function(region_labels, counts, score, evidence, boundary,
                             rest_score, rest_evidence, rest_boundary, alpha,
                             report_p_value = FALSE) {
  ordered <- which(row(score) != col(score), arr.ind = TRUE)
  pair_boundary <- if (is.matrix(boundary)) boundary[ordered] else boundary
  pair_higher <- evidence[ordered]
  pair_lower <- t(evidence)[ordered]
  rest_lower <- 1 - rest_evidence - rest_boundary

  pairwise <- data.frame(
    region = region_labels[ordered[, 1]],
    reference_region = region_labels[ordered[, 2]],
    score = score[ordered],
    evidence = pair_higher,
    evidence_higher = pair_higher,
    evidence_lower = pair_lower,
    call = pi_call(pair_higher, pair_lower, alpha),
    boundary_mass = pair_boundary,
    row.names = NULL
  )
  if (report_p_value) pairwise$p_value <- pi_p_value(pair_higher, pair_lower)
  overall <- data.frame(
    region = region_labels,
    spots = counts,
    score = rest_score,
    evidence = rest_evidence,
    evidence_higher = rest_evidence,
    evidence_lower = rest_lower,
    call = pi_call(rest_evidence, rest_lower, alpha),
    boundary_mass = rest_boundary
  )
  if (report_p_value) overall$p_value <- pi_p_value(rest_evidence, rest_lower)
  list(pairwise = pairwise, overall = overall)
}


#' Doubled posterior probability of the less likely direction (not a
#' frequentist p-value)
#' @noRd
pi_p_value <- function(higher, lower) {
  pmin(1, 2 * pmin(higher, lower))
}


#' "higher" when P(pi > 0.5) >= 1 - alpha / 2, "lower" when P(pi < 0.5) >=
#' 1 - alpha / 2, otherwise "none"
#' @noRd
pi_call <- function(higher, lower, alpha) {
  ifelse(higher >= 1 - alpha / 2, "higher",
         ifelse(lower >= 1 - alpha / 2, "lower", "none"))
}
