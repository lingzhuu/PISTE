# Directional probability engine for the current conditional-iid model.
#
# Notation: for regions r and s, Z_r ~ N(mu_r, tau2_r) is the pathway state
# of a gene in region r. The directional score pi_rs = P(Z_r > Z_s).


# -----------------------------------------------------------------------------
# Conditional directional probability
# -----------------------------------------------------------------------------

#' P(Z_r > Z_s) for (possibly correlated) Gaussian states
#'
#' Vectorised over all arguments; each must be a scalar or share a common
#' length. Values within 32 machine epsilons of 0.5 are snapped to exactly 0.5
#' so that ties are detected reliably downstream.
#'
#' @param mr,ms Means of Z_r and Z_s.
#' @param vr,vs Variances of Z_r and Z_s (must be > 0).
#' @param covariance Cov(Z_r, Z_s).
#' @return Probabilities in [0, 1].
#' @noRd
pi_gaussian_directional <- function(mr, ms, vr, vs, covariance = 0) {
  sizes <- lengths(list(mr, ms, vr, vs, covariance))
  n <- max(sizes)
  if (n == 0L || any(!sizes %in% c(1L, n)))
    stop("Probability inputs must be scalar or have a common positive length.")

  mr <- rep_len(mr, n); ms <- rep_len(ms, n)
  vr <- rep_len(vr, n); vs <- rep_len(vs, n)
  covariance <- rep_len(covariance, n)

  if (any(!is.finite(c(mr, ms, vr, vs, covariance))) || any(vr <= 0 | vs <= 0))
    stop("Means must be finite and variances positive.")
  if (any(abs(covariance) > sqrt(vr * vs) + 1e-8))
    stop("Invalid covariance.")

  vdiff <- vr + vs - 2 * covariance
  if (any(vdiff <= 0)) stop("Degenerate difference variance.")

  # Zr-Zs is Gaussian; a positive difference is directional enrichment.
  ans <- stats::pnorm((mr - ms) / sqrt(vdiff))
  ans[abs(ans - 0.5) <= 32 * .Machine$double.eps] <- 0.5
  pmax(0, pmin(1, ans))
}


#' Conditional probability that the first pathway state exceeds the second
#'
#' Computes pi = P(Z_r > Z_s) for independent Gaussian pathway states
#' Z_r ~ N(mu_r, tau2_r) and Z_s ~ N(mu_s, tau2_s).
#'
#' @param mu_r,mu_s Regional means.
#' @param tau2_r,tau2_s Regional (between-gene) variances.
#' @return Numeric vector of probabilities in [0, 1].
#' @export
pathway_pi <- function(mu_r, mu_s, tau2_r, tau2_s) {
  pi_gaussian_directional(mu_r, mu_s, tau2_r, tau2_s)
}


# -----------------------------------------------------------------------------
# Extraction of the low-dimensional joint posterior from INLA configurations
# -----------------------------------------------------------------------------

#' Covariance of selected latent entries from a sparse precision matrix
#'
#' Select a marginal covariance by sparse solves, not by inverting
#' Q[selected, selected]. Accepts INLA's full, upper- or lower-triangular
#' storage.
#'
#' @param Q Precision matrix (dense or sparse).
#' @param indices Latent indices to extract.
#' @return Dense symmetric covariance matrix of the selected entries.
#' @noRd
pi_selected_covariance <- function(Q, indices) {
  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)

  if (isSymmetric(Q)) {
    Q <- Matrix::forceSymmetric(Q)
  } else if (Matrix::nnzero(Matrix::tril(Q, -1)) == 0) {
    Q <- Matrix::forceSymmetric(Q, uplo = "U")
  } else if (Matrix::nnzero(Matrix::triu(Q, 1)) == 0) {
    Q <- Matrix::forceSymmetric(Q, uplo = "L")
  } else stop("Unsupported nonsymmetric INLA precision storage.")

  # Solve Q X = B where B holds the unit vectors of the selected indices.
  B <- Matrix::sparseMatrix(i = indices, j = seq_along(indices), x = 1,
                            dims = c(nrow(Q), length(indices)))
  solution <- Matrix::solve(Matrix::Cholesky(Q, LDL = FALSE), B, system = "A")
  V <- as.matrix(solution[indices, , drop = FALSE])
  (V + t(V)) / 2
}


#' Extract the CCD mixture of regional means from an INLA fit
#'
#' For every saved CCD hyperparameter configuration, extracts the conditional
#' Gaussian of (baseline intercept, relative regional offsets), maps it to
#' absolute regional means, and records the regional variances tau2_r.
#' Requires the fit to have been run with `control.compute$config = TRUE`.
#'
#' @param fit INLA fit object.
#' @param region_labels Region labels, baseline first.
#' @param region_spots Number of spots per region.
#' @param intercept_name Name of the baseline intercept latent effect.
#' @param regional_offset_names Names of the relative offsets mu_r2..mu_rR.
#' @return A list describing the mixture (configurations, weights, diagnostics).
#' @noRd
pi_extract_joint <- function(fit, region_labels, region_spots,
                             intercept_name, regional_offset_names) {
  started <- proc.time()[["elapsed"]]

  cs <- fit$misc$configs
  if (is.null(cs) || length(cs$config) < 2L)
    stop("Pi inference requires multiple saved CCD configurations.")
  if (!is.null(cs$constr) && length(cs$constr$A) > 0)
    stop("Pi extraction currently requires the package's unconstrained latent model.")

  R <- length(region_labels)

  # Locate the intercept and regional offsets in the latent vector.
  ct <- cs$contents
  starts <- ct$start
  if (isTRUE(cs$.preopt)) {
    # With .preopt the saved Q and mean exclude the predictor blocks, so
    # shift the start offsets accordingly.
    predictors <- ct$tag %in% c("APredictor", "Predictor")
    if (any(predictors) && !identical(which(predictors), seq_len(sum(predictors))))
      stop("Unrecognized INLA predictor layout.")
    starts <- starts - sum(ct$length[predictors])
  }
  selected <- match(c(intercept_name, regional_offset_names), ct$tag)
  if (anyNA(selected) || any(ct$length[selected] != 1L))
    stop("Cannot identify intercept and regional offsets in INLA configuration.")
  indices <- starts[selected]

  # Map (intercept, offset_2, ..., offset_R) to absolute means:
  # mu_1 = intercept, mu_r = intercept + offset_r.
  transform <- matrix(0, R, R)
  transform[, 1L] <- 1
  for (r in seq.int(2L, R)) transform[r, r] <- 1

  # Normalised CCD configuration weights.
  logweights <- vapply(cs$config, function(z) z$log.posterior, numeric(1))
  if (any(!is.finite(logweights))) stop("Nonfinite CCD configuration weights.")
  weights <- exp(logweights - max(logweights))
  weights <- weights / sum(weights)

  configs <- vector("list", length(weights))
  skew_max <- 0
  for (k in seq_along(weights)) {
    z <- cs$config[[k]]
    if (length(z$improved.mean) != nrow(z$Q) ||
        any(indices < 1L | indices > nrow(z$Q)))
      stop("INLA compact precision and mean layout do not agree.")

    # Regional variances tau2_r = 1 / precision of gene_r.
    tags <- paste0("Log precision for gene_r", seq_len(R))
    hidx <- match(tags, names(z$theta))
    if (anyNA(hidx)) stop("Cannot identify regional log precisions by name.")
    tau2 <- exp(-z$theta[hidx])
    if (any(!is.finite(tau2) | tau2 <= 0)) stop("Invalid regional variance.")

    V <- pi_selected_covariance(z$Q, indices)
    if (any(!is.finite(V)) || any(diag(V) <= 0)) stop("Invalid selected covariance.")

    # Transform (baseline intercept, relative offsets) to absolute regional
    # states only for the directional probability integration.
    absolute_mu <- as.numeric(transform %*% z$improved.mean[indices])
    absolute_covariance <- transform %*% V %*% t(transform)
    configs[[k]] <- list(mu = absolute_mu,
                         covariance = absolute_covariance,
                         tau2 = as.numeric(tau2), theta = z$theta)

    sk <- z$skewness[indices]
    if (any(is.finite(sk))) skew_max <- max(skew_max, abs(sk[is.finite(sk)]))
  }

  # Compare mixture moments with INLA's own marginals as a diagnostic.
  mixmean <- Reduce("+", Map(function(z, w) z$mu * w, configs, weights))
  mixsecond <- Reduce("+", Map(function(z, w)
    (z$covariance + tcrossprod(z$mu)) * w, configs, weights))
  mixsd <- sqrt(pmax(0, diag(mixsecond) - mixmean^2))
  reference <- fit$summary.lincomb.derived[paste0("absolute_mu_r", seq_len(R)), ,
                                           drop = FALSE]
  diagnostic <- data.frame(region = region_labels,
                           mixture_mean = mixmean,
                           INLA_marginal_mean = reference$mean,
                           mixture_sd = mixsd,
                           INLA_marginal_sd = reference$sd)

  list(configurations = configs,
       weights = weights,
       region_labels = region_labels,
       region_spots = region_spots,
       estimand = "P(Z_r > Z_s)",
       baseline_region = region_labels[1L],
       approximation = "CCD mixture of improved-mean conditional multivariate Gaussians; no skew transform",
       weight_source = "normalized exp(config$log.posterior), as used by INLA's configuration mixture",
       max_reported_mean_skewness = skew_max,
       marginal_diagnostic = diagnostic,
       extraction_seconds = proc.time()[["elapsed"]] - started)
}


# -----------------------------------------------------------------------------
# Integration of scores and evidence over the mixture
# -----------------------------------------------------------------------------

#' Integrate pi scores and evidence from a low-dimensional CCD Gaussian mixture
#'
#' For every ordered region pair (r, s) computes
#' * **score**: the posterior mean of pi_rs = P(Z_r > Z_s), obtained
#'   analytically per configuration;
#' * **evidence**: the posterior probability that pi_rs > 0.5.
#'
#' "Overall" results compare each region with a spot-count-weighted average of
#' all other regions. With two regions evidence is computed analytically;
#' otherwise it uses antithetic, nested, unrandomised Sobol quadrature that is
#' doubled until two consecutive refinements are stable within `tolerance`.
#' The quadrature is deterministic and never touches the RNG state.
#'
#' @param joint Output of `pi_extract_joint()` (the `joint` element returned by
#'   [screen_joint_fit()]).
#' @param nodes Initial number of quadrature nodes (power of two, >= 16).
#' @param max_nodes Maximum number of nodes (power of two, >= 2 * nodes).
#' @param tolerance Convergence tolerance for refinement differences.
#' @return A list with data frames `pairwise` and `overall`, and a list of
#'   `diagnostics`.
#' @export
pathway_pi_integrate <- function(joint, nodes = 4096L, max_nodes = 65536L,
                                 tolerance = 0.002) {
  started <- proc.time()[["elapsed"]]
  R <- length(joint$region_labels)
  w <- joint$weights
  cfgs <- joint$configurations

  # --- Input validation ------------------------------------------------------
  if (R < 2L || length(w) != length(cfgs) || any(!is.finite(w) | w < 0) ||
      abs(sum(w) - 1) > 1e-8 || length(joint$region_spots) != R ||
      any(!is.finite(joint$region_spots) | joint$region_spots <= 0))
    stop("Invalid low-dimensional mixture.")
  if (nodes < 16 || max_nodes < 2 * nodes ||
      any(log2(c(nodes, max_nodes)) != round(log2(c(nodes, max_nodes)))) ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("Use powers of two, nodes>=16, max_nodes>=2*nodes, and positive tolerance.")

  # Row r of rest_weights: spot-count weights of all regions other than r.
  counts <- joint$region_spots
  rest_weights <- matrix(rep(counts, each = R), R, R)
  diag(rest_weights) <- 0
  rest_weights <- rest_weights / rowSums(rest_weights)

  pairs <- t(combn(seq_len(R), 2))

  # --- Analytic scores (posterior mean of pi) --------------------------------
  # Also precompute a square root of each conditional covariance for
  # the quadrature below.
  score <- matrix(0, R, R)
  roots <- vector("list", length(cfgs))
  for (k in seq_along(cfgs)) {
    z <- cfgs[[k]]
    if (length(z$mu) != R || length(z$tau2) != R ||
        !identical(dim(z$covariance), c(R, R)) ||
        any(!is.finite(c(z$mu, z$tau2, z$covariance))) || any(z$tau2 <= 0) ||
        max(abs(z$covariance - t(z$covariance))) > 1e-8)
      stop("Invalid conditional moments.")

    eig <- eigen(z$covariance, symmetric = TRUE)
    if (min(eig$values) < -1e-8 * max(1, max(eig$values)))
      stop("Conditional covariance is not positive semidefinite.")
    roots[[k]] <- diag(sqrt(pmax(eig$values, 0)), R) %*% t(eig$vectors)

    for (j in seq_len(nrow(pairs))) {
      r <- pairs[j, 1]; s <- pairs[j, 2]
      p <- pi_gaussian_directional(z$mu[r], z$mu[s],
                                   z$tau2[r] + z$covariance[r, r],
                                   z$tau2[s] + z$covariance[s, s],
                                   z$covariance[r, s])
      score[r, s] <- score[r, s] + w[k] * p
      score[s, r] <- score[s, r] + w[k] * (1 - p)
    }
  }
  diag(score) <- NA_real_
  rest_score <- rowSums(replace(score, is.na(score), 0) * rest_weights)

  # --- Two regions: analytic evidence ----------------------------------------
  if (R == 2L) {
    # With two regions, pi_12 > 0.5 iff mu_1 > mu_2. Integrate this
    # event analytically over each CCD conditional Gaussian and mix weights.
    conditional <- lapply(seq_along(cfgs), function(k) {
      z <- cfgs[[k]]
      vdiff <- z$covariance[1, 1] + z$covariance[2, 2] -
        2 * z$covariance[1, 2]
      if (!is.finite(vdiff) || vdiff < 0)
        stop("Invalid conditional mean-difference variance.")
      difference <- z$mu[1] - z$mu[2]
      if (vdiff == 0)
        return(c(first = as.numeric(difference > 0),
                 second = as.numeric(difference < 0),
                 tie = as.numeric(difference == 0)))
      p <- stats::pnorm(difference / sqrt(vdiff))
      c(first = p, second = 1 - p, tie = 0)
    })
    mass <- colSums(do.call(rbind, conditional) * w)
    evidence_12 <- mass[["first"]]
    evidence_21 <- mass[["second"]]
    boundary_mass <- mass[["tie"]]
    evidence <- matrix(c(NA_real_, evidence_21,
                         evidence_12, NA_real_), 2, 2)

    ordered <- which(row(score) != col(score), arr.ind = TRUE)
    pairwise <- data.frame(region = joint$region_labels[ordered[, 1]],
                           reference_region = joint$region_labels[ordered[, 2]],
                           score = score[ordered],
                           evidence = evidence[ordered],
                           boundary_mass = boundary_mass,
                           row.names = NULL)
    overall <- data.frame(region = joint$region_labels,
                          spots = counts,
                          score = rest_score,
                          evidence = c(evidence_12, evidence_21),
                          boundary_mass = boundary_mass)

    return(list(pairwise = pairwise, overall = overall,
                diagnostics = list(
                  converged = TRUE, nodes = 0L,
                  integration = "analytic two-region Gaussian-mixture calculation",
                  refinement = data.frame(nodes = 0L,
                                          evidence_refinement_difference = 0,
                                          analytic_score_check = 0),
                  configurations = length(w),
                  approximation = joint$approximation,
                  estimand = "P(Z_r > Z_s)", boundary = 0.5,
                  strict_inequality = TRUE, deterministic = TRUE,
                  posterior_sampling = FALSE,
                  seconds = proc.time()[["elapsed"]] - started)))
  }

  # --- Three or more regions: adaptive quasi-Monte Carlo evidence ------------
  history <- list()
  stable <- 0L
  n <- as.integer(nodes)
  repeat {
    # Antithetic, nested, unrandomized Sobol quadrature; never changes the RNG state.
    u <- qrng::sobol(n / 2, d = R, randomize = "none", skip = 1)
    normal <- matrix(stats::qnorm(u), nrow = n / 2, ncol = R)
    normal <- normal[rep(seq_len(n / 2), each = 2), , drop = FALSE]
    normal[seq(2, n, 2), ] <- -normal[seq(2, n, 2), , drop = FALSE]

    # "coarse" uses only the first half of the nodes, to measure refinement.
    evidence <- coarse <- tie <- matrix(0, R, R)
    rest_evidence <- rest_coarse <- rest_tie <- numeric(R)
    qscore <- matrix(0, R, R)
    half <- seq_len(n / 2)

    for (k in seq_along(cfgs)) {
      z <- cfgs[[k]]
      # Quadrature draws of the absolute regional means (n x R).
      mu <- sweep(normal %*% roots[[k]], 2, z$mu, "+")
      rest_values <- matrix(0, n, R)

      for (j in seq_len(nrow(pairs))) {
        r <- pairs[j, 1]; s <- pairs[j, 2]
        p <- pathway_pi(mu[, r], mu[, s], z$tau2[r], z$tau2[s])
        evidence[r, s] <- evidence[r, s] + w[k] * mean(p > 0.5)
        evidence[s, r] <- evidence[s, r] + w[k] * mean(p < 0.5)
        coarse[r, s] <- coarse[r, s] + w[k] * mean(p[half] > 0.5)
        coarse[s, r] <- coarse[s, r] + w[k] * mean(p[half] < 0.5)
        tie[r, s] <- tie[r, s] + w[k] * mean(p == 0.5)
        tie[s, r] <- tie[r, s]
        qscore[r, s] <- qscore[r, s] + w[k] * mean(p)
        qscore[s, r] <- qscore[s, r] + w[k] * mean(1 - p)
        rest_values[, r] <- rest_values[, r] + rest_weights[r, s] * p
        rest_values[, s] <- rest_values[, s] + rest_weights[s, r] * (1 - p)
      }

      rest_values[abs(rest_values - 0.5) <= 32 * .Machine$double.eps] <- 0.5
      rest_evidence <- rest_evidence + w[k] * colMeans(rest_values > 0.5)
      rest_coarse <- rest_coarse + w[k] * colMeans(rest_values[half, , drop = FALSE] > 0.5)
      rest_tie <- rest_tie + w[k] * colMeans(rest_values == 0.5)
    }

    # Stop after two consecutive stable refinements, or at max_nodes.
    error <- max(abs(evidence - coarse), abs(rest_evidence - rest_coarse))
    score_check <- max(abs(qscore - score), na.rm = TRUE)
    history[[length(history) + 1L]] <- data.frame(
      nodes = n,
      evidence_refinement_difference = error,
      analytic_score_check = score_check)
    stable <- if (error <= tolerance && score_check <= tolerance) stable + 1L else 0L
    if (stable >= 2L || n >= max_nodes) break
    n <- n * 2L
  }

  converged <- stable >= 2L
  if (!converged) warning("Pi quadrature reached max_nodes without two stable refinements.")

  # --- Assemble output -------------------------------------------------------
  ordered <- which(row(score) != col(score), arr.ind = TRUE)
  pairwise <- data.frame(region = joint$region_labels[ordered[, 1]],
                         reference_region = joint$region_labels[ordered[, 2]],
                         score = score[ordered],
                         evidence = evidence[ordered],
                         boundary_mass = tie[ordered],
                         row.names = NULL)
  overall <- data.frame(region = joint$region_labels,
                        spots = counts,
                        score = rest_score,
                        evidence = rest_evidence,
                        boundary_mass = rest_tie)

  list(pairwise = pairwise, overall = overall,
       diagnostics = list(
         converged = converged, nodes = n, tolerance = tolerance,
         refinement = do.call(rbind, history),
         configurations = length(w),
         approximation = joint$approximation,
         estimand = "P(Z_r > Z_s)", boundary = 0.5, strict_inequality = TRUE,
         boundary_roundoff_tolerance = 32 * .Machine$double.eps,
         deterministic = TRUE, posterior_sampling = FALSE,
         numerical_difference_is_not_error_bound = TRUE,
         seconds = proc.time()[["elapsed"]] - started))
}
