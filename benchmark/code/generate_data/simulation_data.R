# Benchmark raw-count simulation: 3,639 spots x 10,000 genes, for evaluating
# STpathway / GSDensity / PathwayLGM.
#
#   Y_gi ~ NB(mu_gi, theta_g),  mu_gi = c_i * m_g * exp(s_g * b_g * A_pi)
#   Var(Y_gi) = mu_gi + mu_gi^2 / theta_g
#
# (m_g, theta_g, c_i) come from real DLPFC 151673 via estimate_reference.R.
# Genes 1-5,000: 50 disjoint pathways of 100 (5 response sizes x 5 activity
# patterns x 2 effect profiles). Genes 5,001-10,000: background, 1,000 with
# independent radial effects. sigma_eta = 0.

# ---- config ---------------------------------------------------------------------
ref_dir  <- file.path("..", "data", "reference_151673")
out_root <- file.path("..", "data", "simulated")
seeds    <- c(20260926L, 20260929L)   # benchmark seeds

n_genes      <- 10000L
n_pathways   <- 50L
pathway_size <- 100L
n_bg_genes   <- n_genes - n_pathways * pathway_size   # 5,000 background genes
n_bg_svg     <- 1000L    # background genes with a radial spatial effect
matern_ell   <- 0.25     # range of the Matern-3/2 field used for gradient maps

# ---- helpers --------------------------------------------------------------------

# Scale a vector to [0, 1].
unit <- function(v) {
  span <- max(v) - min(v)
  if (!is.finite(span) || span <= 0) stop("Activity map has zero range")
  (v - min(v)) / span
}

# Solve for b so the depth-corrected target/rest expression ratio is 2^L.
coefficient_for_fc <- function(A, target, sign, L) {
  objective <- function(b) {
    sign * (log(mean(exp(sign * b * A[target]))) -
            log(mean(exp(sign * b * A[!target])))) / log(2) - L
  }
  upper <- 2
  while (objective(upper) < 0 && upper < 128) upper <- upper * 2
  if (objective(upper) < 0) stop("Could not calibrate target/rest effect")
  uniroot(objective, c(0, upper), tol = 1e-8)$root
}

# Number of responsive genes, drawn uniformly within the pathway's group.
draw_K <- function(group) {
  bounds <- switch(group,
    K0_4    = 0:4,   K5_10   = 5:10,  K11_30 = 11:30,
    K31_60  = 31:60, K61_100 = 61:100,
    stop("Unknown response group"))
  sample(bounds, 1L)
}

# Number of down-regulated genes among the K responsive ones: each
# responsive gene is down-regulated with probability p_down(K), flat at
# 0.4 up to K = 10, then decaying exponentially.
draw_down_count <- function(K) {
  if (K == 0L) return(0L)
  p_down <- 0.4 * exp(-pmax(0, K - 10) / 18)
  rbinom(1L, K, p_down)
}

# ---- main -----------------------------------------------------------------------

generate_dataset <- function(seed, out) {
  stopifnot(!dir.exists(out))
  set.seed(seed)
  dir.create(out, recursive = TRUE)

  # -- reference from real data: spots, depth factors, per-gene (m, theta) --------
  ref <- read.csv(file.path(ref_dir, "reference_parameters.csv"))
  ref <- ref[is.finite(ref$m) & ref$m > 0 & is.finite(ref$theta) &
               ref$theta > 0 & ref$detected >= 100L, , drop = FALSE]
  spots <- read.csv(file.path(ref_dir, "spots.csv"))
  depth <- read.csv(file.path(ref_dir, "reference_depth.csv"))
  stopifnot(nrow(spots) == 3639L, identical(spots$spot, depth$spot))

  n <- nrow(spots)
  spot <- spots$spot

  # Shift pixel coordinates to the origin and scale to [0, 1]
  # (one common scale for x and y keeps the tissue shape).
  xy <- as.matrix(spots[, c("pixel_col", "pixel_row")])
  xy <- sweep(xy, 2L, apply(xy, 2L, min), "-")
  xy <- xy / max(xy)
  dimnames(xy) <- list(spot, c("x", "y"))

  # Shuffle depth factors across spots; normalize to geometric mean 1.
  c_i <- depth$c_i[sample.int(n)]
  c_i <- c_i / exp(mean(log(c_i)))

  # Each simulated gene inherits one intact empirical (m, theta) pair.
  param <- ref[sample.int(nrow(ref), n_genes, replace = FALSE),
               c("m", "theta"), drop = FALSE]
  gene_ids <- paste0("g", seq_len(n_genes))
  rownames(param) <- gene_ids

  # The first 5,000 genes form 50 pathways of 100 genes each.
  pathways <- split(gene_ids[seq_len(n_pathways * pathway_size)],
                    rep(sprintf("P%02d", seq_len(n_pathways)),
                        each = pathway_size))

  # -- pathway design: 5 response sizes x 5 activity patterns x 2 effect profiles --
  design <- expand.grid(
    response_group = c("K0_4", "K5_10", "K11_30", "K31_60", "K61_100"),
    spatial_class  = c("domain", "gradient", "radial", "hotspot", "interface"),
    effect_profile = c("broad_weak", "narrow_strong"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  design$pathway <- names(pathways)
  design <- design[, c("pathway", "response_group", "spatial_class",
                       "effect_profile")]

  # -- spatial activity maps A_pi --------------------------------------------------
  # All gradient maps share one Matern-3/2 Cholesky factor; compute it once
  # (no random sampling involved in this step).
  d <- as.matrix(dist(xy))
  u <- sqrt(3) * d / matern_ell
  cov_mat <- (1 + u) * exp(-u)
  diag(cov_mat) <- diag(cov_mat) + 1e-7
  gradient_factor <- chol(cov_mat)
  rm(d, u, cov_mat)

  # Draw one activity map of the requested class, scaled to [0, 1].
  make_map <- function(cl) {
    raw <- switch(cl,
      domain = {          # three horizontal bands with jittered boundaries
        b1 <- .33 + runif(1L, -.08, .08) + .04 * sin(2*pi*xy[, 2L])
        b2 <- .65 + runif(1L, -.08, .08) + .04 * sin(2*pi*xy[, 2L])
        c(.1, .8, .3)[1L + (xy[, 1L] > b1) + (xy[, 1L] > b2)]
      },
      gradient = as.vector(crossprod(gradient_factor, rnorm(n))),
      radial = {          # exponential decay around a random spot
        center <- xy[sample.int(n, 1L), ]
        distance <- sqrt((xy[, 1L]-center[1L])^2 + (xy[, 2L]-center[2L])^2)
        exp(-distance / runif(1L, .12, .24))
      },
      hotspot = {         # two Gaussian hotspots of different widths
        center <- xy[sample.int(n, 2L), , drop = FALSE]
        pmax(exp(-((xy[,1L]-center[1L,1L])^2 +
                     (xy[,2L]-center[1L,2L])^2)/(2*.06^2)),
             exp(-((xy[,1L]-center[2L,1L])^2 +
                     (xy[,2L]-center[2L,2L])^2)/(2*.09^2)))
      },
      interface = {       # decaying band along a wavy boundary
        distance <- abs(xy[,1L]-(.5 + .06*sin(2*pi*xy[,2L])))
        as.numeric(distance <= .14) * exp(-distance/.06)
      },
      stop("Unknown spatial class: ", cl))
    unit(raw)
  }

  # Oracle target region of a map, used only to calibrate the log2FC.
  target_mask <- function(A, cl) {
    top20 <- order(A, decreasing = TRUE)[seq_len(ceiling(.2 * length(A)))]
    switch(cl,
      domain    = A == max(A),              # the most active band
      gradient  = seq_along(A) %in% top20,  # the top 20% spots
      radial    = seq_along(A) %in% top20,
      hotspot   = A >= .5,
      interface = A > 0,
      stop("Unknown spatial class: ", cl))
  }

  activity <- matrix(NA_real_, n_pathways, n,
                     dimnames = list(design$pathway, spot))
  target <- matrix(FALSE, n, n_pathways,
                   dimnames = list(spot, design$pathway))
  for (p in seq_len(n_pathways)) {
    activity[p, ] <- make_map(design$spatial_class[p])
    target[, p] <- target_mask(activity[p, ], design$spatial_class[p])
  }

  # -- pathway genes: assign response, direction and effect, then sample counts ---
  responsive <- logical(n_genes)    # is the gene responsive
  direction  <- integer(n_genes)    # +1 up / -1 down / 0 not responsive
  effect_L   <- numeric(n_genes)    # target/rest |log2FC|
  effect_b   <- numeric(n_genes)    # calibrated activity coefficient
  Y <- matrix(0L, n_genes, n, dimnames = list(gene_ids, spot))

  design$K <- integer(n_pathways)
  design$n_up <- integer(n_pathways)
  design$n_down <- integer(n_pathways)
  design$target_spots <- as.integer(colSums(target))
  expected_target <- numeric(n_pathways)
  expected_rest   <- numeric(n_pathways)

  for (p in seq_len(n_pathways)) {
    K <- draw_K(design$response_group[p])
    n_down <- draw_down_count(K)
    ids <- ((p-1L)*pathway_size + 1L):(p*pathway_size)   # the 100 member genes
    chosen <- if (K) sample(ids, K) else integer()
    responsive[chosen] <- TRUE
    direction[chosen] <- 1L                              # start all up
    if (n_down) direction[chosen[sample.int(K, n_down)]] <- -1L
    if (K) {          # |log2FC| from a log-normal, weaker or stronger profile
      if (design$effect_profile[p] == "broad_weak")
        effect_L[chosen] <- rlnorm(K, meanlog = log(.4), sdlog = .3)
      else
        effect_L[chosen] <- rlnorm(K, meanlog = log(.9), sdlog = .4)
    }
    A <- activity[p, ]
    mask <- target[, p]
    for (g in ids) {
      if (responsive[g])
        effect_b[g] <- coefficient_for_fc(A, mask, direction[g], effect_L[g])
      corrected_mean <- param$m[g] * exp(direction[g] * effect_b[g] * A)
      expected_target[p] <- expected_target[p] + mean(corrected_mean[mask])
      expected_rest[p]   <- expected_rest[p]   + mean(corrected_mean[!mask])
      Y[g, ] <- rnbinom(n, size = param$theta[g], mu = c_i * corrected_mean)
    }
    design$K[p] <- K
    design$n_up[p] <- K - n_down
    design$n_down[p] <- n_down
  }

  # -- background genes: 1,000 with independent radial effects, 4,000 flat --------
  background_svg <- logical(n_genes)
  bg_center_x <- rep(NA_real_, n_genes)
  bg_center_y <- rep(NA_real_, n_genes)
  bg_scale    <- rep(NA_real_, n_genes)

  bg_first <- n_pathways * pathway_size + 1L
  bg_ids <- sample.int(n_bg_genes, n_bg_svg) + n_pathways * pathway_size
  background_svg[bg_ids] <- TRUE
  for (g in bg_first:n_genes) {
    if (background_svg[g]) {
      center <- xy[sample.int(n, 1L), ]
      scale <- runif(1L, .08, .25)
      A <- unit(exp(-sqrt((xy[,1L]-center[1L])^2 +
                          (xy[,2L]-center[2L])^2) / scale))
      mask <- seq_len(n) %in% order(A, decreasing = TRUE)[seq_len(ceiling(.2*n))]
      bg_center_x[g] <- center[1L]
      bg_center_y[g] <- center[2L]
      bg_scale[g] <- scale
      direction[g] <- sample(c(-1L, 1L), 1L)
      effect_L[g] <- rlnorm(1L, meanlog = log(.4), sdlog = .3)
      effect_b[g] <- coefficient_for_fc(A, mask, direction[g], effect_L[g])
      corrected_mean <- param$m[g] * exp(direction[g]*effect_b[g]*A)
    } else {
      corrected_mean <- rep(param$m[g], n)
    }
    Y[g, ] <- rnbinom(n, size = param$theta[g], mu = c_i * corrected_mean)
  }

  # -- checks and output -----------------------------------------------------------
  design$all100_expected_log2FC <- log2(expected_target / expected_rest)
  design$all100_expected_fold_change <- expected_target / expected_rest
  design$down_fraction_among_responsive <-
    ifelse(design$K > 0L, design$n_down / design$K, NA_real_)

  # Pathway genes must stay usable downstream: total count >= 10 and
  # detected in at least 20 spots.
  pw <- seq_len(n_pathways * pathway_size)
  if (any(rowSums(Y[pw, , drop = FALSE]) < 10L) ||
      any(rowSums(Y[pw, , drop = FALSE] > 0) < 20L))
    stop("Generated dataset failed pathway-gene eligibility")

  gene_truth <- data.frame(
    gene = gene_ids,
    pathway = c(rep(design$pathway, each = pathway_size),
                rep(NA_character_, n_bg_genes)),
    responsive = responsive, direction = direction,
    target_rest_abs_log2FC = effect_L, activity_coefficient = effect_b,
    background_svg = background_svg,
    background_center_x = bg_center_x,
    background_center_y = bg_center_y,
    background_scale = bg_scale)
  dataset <- list(
    counts = Y, xy = xy, c_i = c_i, pathways = pathways, activity = activity,
    oracle_target = target, gene_truth = gene_truth, pathway_truth = design,
    reference_parameters = param, seed = seed, nuisance_strength = 0,
    generator = "simulate_response_grid.R")
  saveRDS(dataset, file.path(out, "dataset.rds"), compress = "gzip")
  write.csv(design, file.path(out, "pathway_truth.csv"), row.names = FALSE)
  write.csv(gene_truth, file.path(out, "gene_truth.csv"), row.names = FALSE)
  write.csv(data.frame(spot = spot, x = xy[,1L], y = xy[,2L], c_i = c_i),
            file.path(out, "spots.csv"), row.names = FALSE)
  write.csv(data.frame(spot = spot, target*1L, check.names = FALSE),
            file.path(out, "oracle_target_masks.csv"), row.names = FALSE)
}

# Generate one dataset per benchmark seed.
for (s in seeds)
  generate_dataset(s, file.path(out_root, paste0("seed", s)))
