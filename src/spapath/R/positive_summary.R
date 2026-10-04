# =============================================================================
# positive_summary.R
#
# Checked summaries of strictly positive INLA marginals (e.g. precisions).
#
# The ORIGINAL INLA-smoothed density is integrated in log coordinates
# z = log(x). The Jacobian changes the integration variable, not the marginal
# interpolant:
#
#   integral f(x) g(x) dx  =  integral f(exp(z)) exp(z) g(exp(z)) dz
#
# Every result is cross-checked against an independent dense trapezoidal
# integration of the same density.
# =============================================================================


# Number of sub-intervals used when the whole-range integral fails or fails
# the check.
POSITIVE_SUMMARY_SEGMENTS <- 64L

# Where failing marginals are saved for offline debugging.
POSITIVE_SUMMARY_FAILURE_DIR <- "logs/positive_marginal_failures"


#' Summarise a strictly positive posterior marginal
#'
#' Computes the posterior mean of `transform(x)` and the transformed 2.5%,
#' 50% and 97.5% quantiles of x by adaptive integration on the log scale. If
#' the independent trapezoidal check disagrees by a relative error of 5e-4 or
#' more, the calculation is repeated once with the range split into 64
#' segments; if it still disagrees an error is raised. A whole-range integral
#' that fails outright is split the same way and still checked.
#'
#' @param m Two-column INLA marginal (x, density) with all x > 0.
#' @param transform Function applied to x (e.g. `function(x) 1 / x` to turn a
#'   precision into a variance).
#' @param decreasing TRUE when `transform` is decreasing, so the lower and
#'   upper quantiles are swapped.
#' @param split_integral Internal: integrate over 64 sub-intervals.
#' @return Named numeric vector `c(mean, lower, median, upper)` with attributes
#'   `independent_relative_error` and `integration_segments`.
#' @noRd
positive_marginal_summary <- function(m, transform = identity,
                                      decreasing = FALSE,
                                      split_integral = FALSE) {
  stopifnot(all(is.finite(m)), all(m[, 1] > 0), all(m[, 2] >= 0))

  # Smoothed interpolant of the marginal: f$fun(x) is the log density and
  # f$range its support.
  f <- INLA:::inla.sfmarginal(INLA::inla.smarginal(m))
  log_range <- log(f$range)

  # Density on the log scale, including the Jacobian exp(z).
  log_scale_density <- function(z) exp(f$fun(exp(z)) + z)

  # Same integral on POSITIVE_SUMMARY_SEGMENTS sub-intervals.
  split_mass <- function(integrand, a, b) {
    cuts <- seq(a, b, length.out = POSITIVE_SUMMARY_SEGMENTS + 1L)
    pieces <- vapply(seq_len(POSITIVE_SUMMARY_SEGMENTS), function(i) {
      stats::integrate(integrand, cuts[i], cuts[i + 1L],
                       rel.tol = 1e-9, abs.tol = 1e-12,
                       subdivisions = 2000L, stop.on.error = TRUE)$value
    }, 0.)
    sum(pieces)
  }

  # TRUE once a whole-range integral has failed and been split instead.
  split_fallback <- FALSE

  # Integral over z in [a, b] of log_scale_density(z) * weight(exp(z)).
  mass <- function(a, b, weight = function(x) rep(1, length(x))) {
    if (a == b) return(0)
    integrand <- function(z) log_scale_density(z) * weight(exp(z))

    # Integrate the same density on smaller intervals when a whole-range
    # adaptive integral disagrees with the independent dense-grid check.
    if (split_integral) return(split_mass(integrand, a, b))

    # A whole-range integral that fails outright (e.g. "extremely bad
    # integrand behaviour") is also split; the result is still checked below.
    tryCatch(
      stats::integrate(integrand, a, b,
                       rel.tol = 1e-9, subdivisions = 2000L,
                       stop.on.error = TRUE)$value,
      error = function(e) {
        split_fallback <<- TRUE
        tryCatch(split_mass(integrand, a, b), error = function(e2) {
          positive_summary_save_failure(m, a, b, e2)
          stop(e2)
        })
      }
    )
  }

  # --- Adaptive integration --------------------------------------------------
  # norm: total mass of the smoothed density (the normalising constant).
  norm <- mass(log_range[1], log_range[2])
  stopifnot(is.finite(norm), norm > 0)

  cdf <- function(z) mass(log_range[1], z) / norm
  log_quantiles <- vapply(c(.025, .5, .975), function(p) {
    stats::uniroot(function(z) cdf(z) - p, log_range, tol = 1e-9)$root
  }, 0.)
  quantiles <- exp(log_quantiles)
  if (decreasing) quantiles <- rev(quantiles)

  values <- c(
    mean = mass(log_range[1], log_range[2], transform) / norm,
    lower = transform(quantiles[1]),
    median = transform(quantiles[2]),
    upper = transform(quantiles[3])
  )
  stopifnot(all(is.finite(values)),
            values["lower"] <= values["median"],
            values["median"] <= values["upper"])

  # --- Independent check -----------------------------------------------------
  check <- positive_summary_trapezoid(log_scale_density, log_range, transform,
                                      decreasing)
  err <- max(abs(check - values) / pmax(abs(values), 1e-12))

  # err: largest relative disagreement between the two methods (limit 5e-4).
  if (is.finite(err) && err >= 5e-4 &&
      !split_integral) {
    return(positive_marginal_summary(m, transform, decreasing,
                                     split_integral = TRUE))
  }
  stopifnot(is.finite(err), err < 5e-4)

  attr(values, "independent_relative_error") <- err
  attr(values, "integration_segments") <-
    if (split_integral || split_fallback) POSITIVE_SUMMARY_SEGMENTS else 1L
  values
}


#' Independent trapezoidal version of the same summary
#'
#' Integrates on a dense grid of 32769 points with an explicit CDF of 0 at the
#' lower end of the support; quantiles are found by linear interpolation of
#' the CDF, with no extrapolation.
#'
#' @return Unnamed vector (mean, lower, median, upper) on the transformed
#'   scale.
#' @noRd
positive_summary_trapezoid <- function(log_scale_density, log_range, transform,
                                       decreasing) {
  trapezoid <- function(y, dx) dx * (utils::head(y, -1) + utils::tail(y, -1)) / 2

  z <- seq(log_range[1], log_range[2], length.out = 32769L)
  density <- log_scale_density(z)
  dx <- diff(z)

  cdf <- c(0, cumsum(trapezoid(density, dx)))
  area <- utils::tail(cdf, 1)

  quantiles <- exp(approx(cdf / area, z, xout = c(.025, .5, .975),
                          ties = "ordered")$y)
  if (decreasing) quantiles <- rev(quantiles)

  transformed_density <- density * transform(exp(z))
  c(sum(trapezoid(transformed_density, dx)) / area, transform(quantiles))
}


#' Save a marginal whose integration failed, for offline debugging
#' @noRd
positive_summary_save_failure <- function(m, a, b, error) {
  dir.create(POSITIVE_SUMMARY_FAILURE_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(m = m, a = a, b = b, message = conditionMessage(error)),
          file.path(POSITIVE_SUMMARY_FAILURE_DIR,
                    paste0("marginal_", Sys.getpid(), ".rds")))
}
