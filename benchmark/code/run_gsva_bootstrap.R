# Target-vs-rest tests of the GSVA spot scores from run_gsva_score.R.
# The statistic is the mean score of the target spots minus that of the rest.
#   block_bootstrap_test  two-sided spatial block bootstrap: coordinates are
#                         scaled to the unit square and cut into square blocks
#                         of side `width`; blocks are resampled with replacement
#   permutation_test      two-sided permutation test: region labels are
#                         shuffled across spots, ignoring spatial dependence
# P-values are BH-adjusted across pathways (q).

# ---- config ---------------------------------------------------------------------
work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")
out_dir <- file.path(work_dir, "GSVA")
width <- .1
B <- 1999
seed <- 42

# ---- functions ------------------------------------------------------------------
# scores: spots x pathways; target: logical, one value per spot, or a
# spots x pathways matrix when each pathway has its own target region.
target_minus_rest <- function(scores, target) {
  colSums(scores * target) / colSums(target) -
    colSums(scores * !target) / colSums(!target)
}

as_target <- function(target, scores) {
  if (is.null(dim(target))) target <- matrix(target, nrow(scores), ncol(scores))
  stopifnot(dim(target) == dim(scores), is.logical(target))
  target
}

block_bootstrap_test <- function(scores, xy, target, width = .1, B = 1999, seed = 42) {
  target <- as_target(target, scores)
  xy <- sweep(xy, 2, apply(xy, 2, min))
  xy <- xy / max(xy)
  key <- paste(floor(xy[, 1] / width), floor(xy[, 2] / width))
  block <- match(key, unique(key))
  n_block <- max(block)
  # Block sums, so one draw is a weighted sum over blocks.
  rest <- !target
  nt <- rowsum(target * 1, block)
  nr <- rowsum(rest * 1, block)
  st <- rowsum(scores * target, block)
  sr <- rowsum(scores * rest, block)
  observed <- target_minus_rest(scores, target)

  set.seed(seed)
  boot <- matrix(NA, B, ncol(scores), dimnames = list(NULL, colnames(scores)))
  b <- 0
  rejected <- 0
  while (b < B) {
    w <- tabulate(sample(n_block, n_block, replace = TRUE), n_block)
    bt <- colSums(nt * w)
    br <- colSums(nr * w)
    # Redraw when a resample has no target or no rest spot for some pathway.
    if (any(bt == 0 | br == 0)) {
      rejected <- rejected + 1
      stopifnot(rejected < 1e5)
      next
    }
    b <- b + 1
    boot[b, ] <- colSums(st * w) / bt - colSums(sr * w) / br
  }

  centered <- sweep(boot, 2, observed)
  p <- (1 + colSums(abs(centered) >= rep(abs(observed), each = B))) / (B + 1)
  ci <- apply(boot, 2, quantile, c(.975, .025))
  data.frame(pathway = colnames(scores), target_minus_rest = observed,
             ci_lower = 2 * observed - ci[1, ], ci_upper = 2 * observed - ci[2, ],
             p = p, q = p.adjust(p, "BH"), n_blocks = n_block,
             target_blocks = colSums(nt > 0), rest_blocks = colSums(nr > 0),
             replicates = B, rejected_draws = rejected, row.names = NULL)
}

permutation_test <- function(scores, target, B = 1999, seed = 42) {
  target <- as_target(target, scores)
  observed <- target_minus_rest(scores, target)
  set.seed(seed)
  perm <- t(replicate(B, target_minus_rest(scores, target[sample(nrow(target)), , drop = FALSE])))
  p <- (1 + colSums(abs(perm) >= rep(abs(observed), each = B))) / (B + 1)
  data.frame(pathway = colnames(scores), target_minus_rest = observed,
             p = p, q = p.adjust(p, "BH"), replicates = B, row.names = NULL)
}

# ---- load -----------------------------------------------------------------------
# Rows: spots, columns: pathways, as written by run_gsva_score.R.
scores <- as.matrix(read.csv(file.path(out_dir, "gsva_scores.csv"),
                             row.names = 1, check.names = FALSE))
# One row per spot: spot, x, y, region ("target" or "rest").
spots <- read.csv(file.path(data_dir, "spots.csv"))
spots <- spots[match(rownames(scores), spots$spot), ]
stopifnot(!anyNA(spots$spot), spots$region %in% c("target", "rest"))
xy <- as.matrix(spots[, c("x", "y")])
target <- spots$region == "target"

# ---- run ------------------------------------------------------------------------
bootstrap <- block_bootstrap_test(scores, xy, target, width, B, seed)
permutation <- permutation_test(scores, target, B, seed)

# ---- save -----------------------------------------------------------------------
write.csv(bootstrap, file.path(out_dir, "gsva_block_bootstrap.csv"), row.names = FALSE)
write.csv(permutation, file.path(out_dir, "gsva_permutation.csv"), row.names = FALSE)
