# SpaPath evidence on the simulated data sets under two rules, from the saved fits (no refits).
#   old:   evidence = P(pi > 0.5 | Y); call up if > 0.95, down if < 0.05.
#   delta: call up if P(pi > 0.5 + delta | Y) > 0.95, down if P(pi < 0.5 - delta | Y) > 0.95.
# pi is the conditional target-vs-rest pi of SpaPath 0.6.3,
#   pi = Phi(D / sqrt(tau2_t + tau2_r + kappa2 + V)),  D = mu_target - mu_rest,
# and D is Gaussian within each CCD configuration, so both probabilities are
# closed-form per configuration and mixed with the CCD weights.
# Usage: Rscript spapath_evidence_delta.R SIM_ROOT [delta]
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
delta <- if (length(args) > 1) as.numeric(args[[2]]) else 0.1

tail_probs <- function(sm, delta) {
  stopifnot(identical(sm$metadata$region_labels, c("rest", "target")))
  comps <- sm$joint$configurations; w <- sm$joint$weights
  per <- vapply(comps, function(k) {
    m <- k$mu[2] - k$mu[1]
    s <- sqrt(max(k$covariance[1, 1] + k$covariance[2, 2] - 2 * k$covariance[1, 2], 0))
    den <- sqrt(sum(k$tau2) + k$kappa2[2, 1] + k$field_variance[2, 1])
    c(gt_half = pnorm(m / s),
      up = 1 - pnorm((qnorm(0.5 + delta) * den - m) / s),
      down = pnorm((qnorm(0.5 - delta) * den - m) / s))
  }, numeric(3))
  colSums(t(per) * w)
}

truth <- read.csv(file.path(root, "pi_true.csv"))
rows <- list()
for (s in 2020:2026) {
  base <- file.path(root, paste0("seed", s))
  pt <- read.csv(file.path(base, "pathway_truth.csv"))
  fit_dir <- file.path(base, "SpaPath_dev_0.6.3", "Gaussian")
  st <- read.csv(file.path(fit_dir, "status.csv"))
  done <- unique(st$pathway[st$status == "ok"])
  for (p in intersect(pt$pathway, done)) {
    sm <- readRDS(file.path(fit_dir, p, "screen_summary.rds"))
    ov <- sm$pi$overall
    tp <- tail_probs(sm, delta)
    i <- match(p, pt$pathway)
    rows[[length(rows) + 1]] <- data.frame(
      seed = s, pathway = p, response_group = pt$response_group[i], K = pt$K[i],
      effect_profile = pt$effect_profile[i], spatial_class = pt$spatial_class[i],
      pi_true = truth$pi_true[truth$seed == s & truth$pathway == p],
      pi = ov$score[ov$region == "target"],
      evidence_old = ov$evidence[ov$region == "target"],
      evidence_old_check = tp[["gt_half"]],
      p_up_delta = tp[["up"]], p_down_delta = tp[["down"]])
  }
}
out <- do.call(rbind, rows)
call <- function(up, down) ifelse(up, "up", ifelse(down, "down", "none"))
out$call_old <- call(out$evidence_old > 0.95, out$evidence_old < 0.05)
out$call_delta <- call(out$p_up_delta > 0.95, out$p_down_delta > 0.95)
out$true_direction <- ifelse(out$K == 0, "none", ifelse(out$pi_true > 0.5, "up", "down"))
cat("Max |evidence_old - closed form| =", max(abs(out$evidence_old - out$evidence_old_check)), "\n")
write.csv(out, file.path(root, sprintf("spapath_evidence_delta%.2f.csv", delta)), row.names = FALSE)

grp <- c("K0", "K1_10", "K11_30", "K31_60", "K61_100")
summ <- do.call(rbind, lapply(grp, function(g) {
  d <- out[out$response_group == g, ]
  f <- function(cl) c(called = sum(cl != "none"),
                      correct = sum(cl != "none" & cl == d$true_direction),
                      wrong_direction = sum(cl != "none" & d$true_direction != "none" & cl != d$true_direction))
  o <- f(d$call_old); n <- f(d$call_delta)
  data.frame(K_group = g, n = nrow(d),
             old_called = o[["called"]], old_correct = o[["correct"]], old_wrong_dir = o[["wrong_direction"]],
             delta_called = n[["called"]], delta_correct = n[["correct"]], delta_wrong_dir = n[["wrong_direction"]])
}))
write.csv(summ, file.path(root, sprintf("spapath_evidence_delta%.2f_summary.csv", delta)), row.names = FALSE)
print(summ, row.names = FALSE)
k0 <- out[out$K == 0, ]
cat("\nK = 0 sets with a call (false positives): old", sum(k0$call_old != "none"), "/", nrow(k0),
    " delta", sum(k0$call_delta != "none"), "/", nrow(k0), "\n")
print(k0[k0$call_old != "none" | k0$call_delta != "none",
         c("seed", "pathway", "pi", "evidence_old", "p_up_delta", "p_down_delta", "call_old", "call_delta")],
      row.names = FALSE, digits = 3)
