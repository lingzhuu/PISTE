# Checked summaries of strictly positive INLA marginals.
#
# The ORIGINAL INLA-smoothed density is integrated in log coordinates.
# The Jacobian changes the integration variable, not the marginal interpolant.


#' Summarise a strictly positive posterior marginal
#'
#' Computes the posterior mean of `transform(x)` and the 2.5%, 50% and 97.5%
#' quantiles of a positive marginal by adaptive integration on the log scale.
#' The result is cross-checked against an independent dense trapezoidal
#' integration; if the two disagree by a relative error of 5e-4 or more, the
#' calculation is repeated with the range split into 64 segments.
#'
#' @param m Two-column INLA marginal (x, density) with all x > 0.
#' @param transform Function applied to x (e.g. `function(x) 1 / x` to turn a
#'   precision into a variance).
#' @param decreasing TRUE when `transform` is decreasing, so the lower and
#'   upper quantiles must be swapped.
#' @param split_integral Internal: integrate over 64 sub-intervals.
#' @return Named numeric vector `c(mean, lower, median, upper)` with attributes
#'   `independent_relative_error` and `integration_segments`.
#' @noRd
positive_marginal_summary <- function(m, transform = identity,
                                      decreasing = FALSE,
                                      split_integral = FALSE) {
  stopifnot(all(is.finite(m)), all(m[, 1] > 0), all(m[, 2] >= 0))

  # Smoothed log-density interpolant of the marginal and its support.
  f <- INLA:::inla.sfmarginal(INLA::inla.smarginal(m))
  zr <- log(f$range)

  # Integral over z = log(x) in [a, b] of density(x) * x * trans(x).
  mass <- function(a, b, trans = function(x) rep(1, length(x))) {
    if (a == b) return(0)

    if (split_integral) {
      # Integrate the same density on smaller intervals when a whole-range
      # adaptive integral disagrees with the independent dense-grid check.
      cuts <- seq(a, b, length.out = 65L)
      return(sum(vapply(seq_len(64L), function(i)
        stats::integrate(
          function(z) exp(f$fun(exp(z)) + z) * trans(exp(z)),
          cuts[i], cuts[i + 1L],
          rel.tol = 1e-9, abs.tol = 1e-12,
          subdivisions = 2000L, stop.on.error = TRUE
        )$value, 0.)))
    }

    tryCatch(
      stats::integrate(
        function(z) exp(f$fun(exp(z)) + z) * trans(exp(z)),
        a, b,
        rel.tol = 1e-9, subdivisions = 2000L, stop.on.error = TRUE
      )$value,
      error = function(e) {
        # Save the failing marginal for offline debugging, then re-raise.
        dir.create("logs/positive_marginal_failures",
                   recursive = TRUE, showWarnings = FALSE)
        saveRDS(
          list(m = m, a = a, b = b, message = conditionMessage(e)),
          file.path("logs/positive_marginal_failures",
                    paste0("marginal_", Sys.getpid(), ".rds"))
        )
        stop(e)
      }
    )
  }

  # --- Adaptive integration --------------------------------------------------
  norm <- mass(zr[1], zr[2])
  stopifnot(is.finite(norm), norm > 0)

  qlog <- vapply(c(.025, .5, .975), function(p)
    stats::uniroot(function(z) mass(zr[1], z) / norm - p, zr, tol = 1e-9)$root,
    0.)
  q <- exp(qlog)
  if (decreasing) q <- rev(q)

  values <- c(
    mean   = mass(zr[1], zr[2], transform) / norm,
    lower  = transform(q[1]),
    median = transform(q[2]),
    upper  = transform(q[3])
  )
  stopifnot(all(is.finite(values)),
            values["lower"] <= values["median"],
            values["median"] <= values["upper"])

  # --- Independent check -----------------------------------------------------
  # Independent trapezoidal integration with explicit CDF(1)=0; no extrapolation.
  z <- seq(zr[1], zr[2], length.out = 32769L)
  den <- exp(f$fun(exp(z)) + z)
  dx <- diff(z)
  cdf <- c(0, cumsum(dx * (utils::head(den, -1) + utils::tail(den, -1)) / 2))
  area <- utils::tail(cdf, 1)

  q2 <- exp(approx(cdf / area, z, xout = c(.025, .5, .975),
                   ties = "ordered")$y)
  if (decreasing) q2 <- rev(q2)

  fy <- den * transform(exp(z))
  other <- c(sum(dx * (utils::head(fy, -1) + utils::tail(fy, -1)) / 2) / area,
             transform(q2))
  err <- max(abs(other - values) / pmax(abs(values), 1e-12))

  # Retry once on split intervals if the two methods disagree.
  if (is.finite(err) && err >= 5e-4 && !split_integral)
    return(positive_marginal_summary(m, transform, decreasing,
                                     split_integral = TRUE))
  stopifnot(is.finite(err), err < 5e-4)

  attr(values, "independent_relative_error") <- err
  attr(values, "integration_segments") <- if (split_integral) 64L else 1L
  values
}
