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
# The model is fitted with INLA.
#
# Regions of the map are joint excursion sets (Bolin & Lindgren 2015): the
# largest set of spots on which z > 0 everywhere with posterior probability
# 1 - alpha, and likewise z < 0, computed by the excursions package with its
# QC method.

# Meshes on the scaled coordinates (max.edge inside / outside, cutoff, offset
# inside / outside).
FIELD_MESH <- list(max_edge = c(.02, .1), cutoff = .01, offset = c(.05, .2))
GENE_FIELD_MESH <- list(max_edge = c(.1, .3), cutoff = .3 / 9, offset = c(.15, .45))


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
#' `gene_field = TRUE` a spatial field of each gene. Pass the result to
#' [field_regions()].
#'
#' @param Y Spots x genes matrix of normalised expression, with spot rownames.
#' @param coords Spots x 2 matrix of spot coordinates, in the spot order of `Y`.
#' @param gene_field Add a spatial field of each gene, which absorbs
#'   gene-specific spatial patterns. `FALSE` is about 50 times faster but
#'   less accurate: spatial patterns of single genes that are not
#'   shared by the gene set can then appear as regions of the pathway.
#' @param threads INLA `num.threads`.
#' @return A list with `inla` (the INLA fit), `node` (the mesh node of each
#'   spot), `spots`, `sigma_e`, `hyperparameters` (posterior summaries of the
#'   noise precision and the range and sd of each field; ranges on the scaled
#'   coordinates) and `metadata`.
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

  hyper <- fit$summary.hyperpar[, c("mean", "0.025quant", "0.5quant", "0.975quant")]
  names(hyper) <- c("mean", "lower", "median", "upper")
  rownames(hyper) <- c("Precision for the Gaussian observations" = "noise_precision",
                       "Range for field" = "range_z", "Stdev for field" = "sd_z",
                       "Range for gene_field" = "range_w", "Stdev for gene_field" = "sd_w")[rownames(hyper)]
  # The mesh is built on the spots, so each spot is a mesh node; spots closer
  # than the cutoff share the node with the largest basis weight.
  list(inla = fit, node = apply(A, 1, which.max), spots = rownames(Y),
       sigma_e = 1 / sqrt(hyper["noise_precision", "mean"]), hyperparameters = hyper,
       metadata = list(N = N, K = K, gene_labels = colnames(Y), gene_field = gene_field,
                       mesh_nodes = mesh$n, threads = threads,
                       model_version = as.character(utils::packageVersion("PISTE")),
                       fit_seconds = proc.time()[["elapsed"]] - started))
}


#' Regions of higher and lower pathway activity
#'
#' Calls the spots where the pathway field z is above or below the baseline of
#' the pathway, jointly: a spot is called `"higher"` when it belongs to the
#' largest set of spots on which z > 0 everywhere with posterior probability
#' at least 1 - `alpha`, and `"lower"` likewise for z < 0. The
#' sets are computed with the QC method of [excursions::excursions.inla()].
#'
#' @param fit Output of [field_fit()].
#' @param alpha One minus the joint posterior probability of a region.
#' @return Data frame with one row per spot: `spot`, `z_mean`, `z_sd`,
#'   `p_higher` = P(z > 0), `p_lower` = P(z < 0) and `call`.
#' @export
field_regions <- function(fit, alpha = .05) {
  threads <- as.integer(sub(":.*", "", fit$metadata$threads))
  excursion <- function(type)
    excursions::excursions.inla(fit$inla, name = "field", ind = unique(fit$node), method = "QC", alpha = alpha,
                                u = 0, type = type, max.threads = threads, seed = 42)
  higher <- excursion(">")
  lower <- excursion("<")
  z <- fit$inla$summary.random$field[fit$node, ]
  data.frame(spot = fit$spots, z_mean = z$mean, z_sd = z$sd,
             p_higher = higher$rho[fit$node], p_lower = lower$rho[fit$node],
             call = ifelse(higher$E[fit$node] == 1, "higher", ifelse(lower$E[fit$node] == 1, "lower", "none")),
             stringsAsFactors = FALSE)
}
