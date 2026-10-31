# Spatial map of a pathway: one shared spatial field for all genes of a gene set.
#
#   y[i, g] = a_g + z(s_i) + w_g(s_i) + e[i, g]
#
#   a_g     ~ N(0, 100)                      gene intercept
#   z       Matern field (alpha = 2) on the spot mesh, sum-to-zero; the pathway map
#   w_g     Matern field of gene g on a coarser mesh, one range and sd for all
#           genes; only with gene_field = TRUE
#   e[i, g] ~ N(0, sigma_e^2)
#
# Both fields have PC priors P(range < 0.1) = 0.5 and P(sd > 1) = 0.05, with
# coordinates scaled so that the longer side of the section is 1; the noise
# precision has the PC prior P(sigma_e > 2) = 0.05. Every gene enters with
# weight +1, so genes of a gene set that move in opposite directions cancel.
# The model is fitted with INLA and z is sampled jointly at the spots.
#
# Regions of the map are joint excursion sets: the largest set of spots on
# which z > u everywhere with posterior probability 1 - alpha, and likewise
# z < -u, with the minimum effect u = d * sigma_e.

# Meshes on the scaled coordinates (max.edge inside / outside, cutoff, offset
# inside / outside) and the number of joint posterior draws of z.
FIELD_MESH <- list(max_edge = c(.02, .1), cutoff = .01, offset = c(.05, .2))
GENE_FIELD_MESH <- list(max_edge = c(.1, .3), cutoff = .3 / 9, offset = c(.15, .45))
FIELD_DRAWS <- 1000


field_spde <- function(mesh) {
  INLA::inla.spde2.pcmatern(mesh, alpha = 2, prior.range = c(.1, .5),
                            prior.sigma = c(1, .05), constr = TRUE)
}


field_mesh <- function(xy, m) {
  fmesher::fm_mesh_2d(loc = xy, max.edge = m$max_edge, cutoff = m$cutoff, offset = m$offset)
}


#' Fit the spatial field of one gene set
#'
#' Fits one spatial field z shared by all genes of the gene set, and with
#' `gene_field = TRUE` a spatial field of each gene, and samples z jointly at
#' the spots. Pass the result to [field_regions()].
#'
#' @param Y Spots x genes matrix of normalised expression, with spot rownames.
#' @param coords Spots x 2 matrix of spot coordinates, in the spot order of `Y`.
#' @param gene_field Add a spatial field of each gene, which absorbs
#'   gene-specific spatial patterns. `FALSE` is about 50 times faster but
#'   less accurate: spatial patterns of single genes that are not
#'   shared by the gene set can then appear as regions of the pathway.
#' @param threads INLA `num.threads`.
#' @return A list with `draws` (spots x draws matrix of z), `sigma_e`,
#'   `hyperparameters` (posterior summaries of the noise precision and the
#'   range and sd of each field; ranges on the scaled coordinates) and
#'   `metadata`.
#' @export
field_fit <- function(Y, coords, gene_field = TRUE, threads = "8:1") {
  started <- proc.time()[["elapsed"]]
  Y <- as.matrix(Y)
  xy <- as.matrix(coords)
  xy <- sweep(xy, 2, apply(xy, 2, min))
  xy <- xy / max(xy)
  N <- nrow(Y); K <- ncol(Y)

  mesh <- field_mesh(xy, FIELD_MESH)
  spde <- field_spde(mesh)
  A <- fmesher::fm_basis(mesh, loc = xy)
  gene <- rep(seq_len(K), each = N)
  A_list <- list(do.call(rbind, rep(list(A), K)), 1)
  effects <- list(field = seq_len(spde$n.spde), gene = gene)
  formula <- y ~ -1 + f(gene, model = "iid", n = K, hyper = list(prec = list(initial = log(.01), fixed = TRUE))) +
    f(field, model = spde)
  if (gene_field) {
    mesh_w <- field_mesh(xy, GENE_FIELD_MESH)
    spde_w <- field_spde(mesh_w)
    A_list <- list(A_list[[1]], Matrix::bdiag(rep(list(fmesher::fm_basis(mesh_w, loc = xy)), K)), 1)
    effects <- list(field = seq_len(spde$n.spde),
                    INLA::inla.spde.make.index("gene_field", n.spde = spde_w$n.spde, n.repl = K),
                    gene = gene)
    formula <- stats::update(formula, . ~ . + f(gene_field, model = spde_w, replicate = gene_field.repl))
  }
  stack <- INLA::inla.stack(data = list(y = as.vector(Y)), A = A_list, effects = effects)
  fit <- INLA::inla(formula, family = "gaussian", data = INLA::inla.stack.data(stack),
                    control.family = list(hyper = list(prec = list(prior = "pc.prec", param = c(2, .05)))),
                    control.predictor = list(A = INLA::inla.stack.A(stack), compute = FALSE),
                    control.compute = list(config = TRUE),
                    control.inla = list(strategy = "gaussian", int.strategy = "ccd"),
                    inla.mode = "experimental", num.threads = threads)
  samples <- INLA::inla.posterior.sample(FIELD_DRAWS, fit, selection = list(field = 0),
                                         num.threads = threads, seed = 42, skew.corr = FALSE)
  draws <- as.matrix(A %*% do.call(cbind, lapply(samples, function(s) s$latent[, 1])))
  rownames(draws) <- rownames(Y)

  hyper <- fit$summary.hyperpar[, c("mean", "0.025quant", "0.5quant", "0.975quant")]
  names(hyper) <- c("mean", "lower", "median", "upper")
  rownames(hyper) <- c("Precision for the Gaussian observations" = "noise_precision",
                       "Range for field" = "range_z", "Stdev for field" = "sd_z",
                       "Range for gene_field" = "range_w", "Stdev for gene_field" = "sd_w")[rownames(hyper)]
  list(draws = draws, sigma_e = 1 / sqrt(hyper["noise_precision", "mean"]), hyperparameters = hyper,
       metadata = list(N = N, K = K, gene_labels = colnames(Y), gene_field = gene_field,
                       mesh_nodes = mesh$n, threads = threads,
                       model_version = as.character(utils::packageVersion("PathwayLGM")),
                       fit_seconds = proc.time()[["elapsed"]] - started))
}


# Largest set of spots, in decreasing order of P(B), on which B holds in all
# of a fraction 1 - alpha of the draws.
excursion_set <- function(B, alpha) {
  prob <- rowMeans(B)
  ord <- order(prob, decreasing = TRUE)
  run <- apply(B[ord, , drop = FALSE], 2, function(v) {
    f <- which(!v)
    if (length(f)) f[1] - 1 else length(v)
  })
  m <- sum(vapply(seq_len(nrow(B)), function(k) mean(run >= k) >= 1 - alpha, logical(1)))
  in_set <- logical(nrow(B))
  in_set[ord[seq_len(m)]] <- TRUE
  list(in_set = in_set, prob = prob)
}


#' Regions of higher and lower pathway activity
#'
#' Calls the spots where the pathway field z exceeds the minimum effect
#' u = `d` * sigma_e, jointly: a spot is called `"higher"` when it belongs to
#' the largest set of spots on which z > u everywhere with posterior
#' probability at least 1 - `alpha`, and `"lower"` likewise for z < -u.
#'
#' @param fit Output of [field_fit()].
#' @param d Minimum effect, in units of the noise sd sigma_e.
#' @param alpha One minus the joint posterior probability of a region.
#' @return Data frame with one row per spot: `spot`, `z_mean`, `z_sd`,
#'   `p_higher` = P(z > u), `p_lower` = P(z < -u) and `call`.
#' @export
field_regions <- function(fit, d = .1, alpha = .05) {
  u <- d * fit$sigma_e
  higher <- excursion_set(fit$draws > u, alpha)
  lower <- excursion_set(fit$draws < -u, alpha)
  data.frame(spot = rownames(fit$draws), z_mean = rowMeans(fit$draws),
             z_sd = apply(fit$draws, 1, stats::sd), p_higher = higher$prob, p_lower = lower$prob,
             call = ifelse(higher$in_set, "higher", ifelse(lower$in_set, "lower", "none")),
             stringsAsFactors = FALSE)
}
