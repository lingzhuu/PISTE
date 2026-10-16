# Summarise the FindMarkers + fgsea runs on simulated_v4: within each seed, BH over the 50
# pathway p-values gives q; detected if q < 0.05; direction from the sign of NES, expected up
# when the responsive genes are mostly up. Writes per-pathway and per-group csv files.
# Usage: Rscript summarize_findmarkers_fgsea.R <simulated_v4 root> <out_dir of run_findmarkers_fgsea.R>
args <- commandArgs(TRUE)
root <- args[1]; out <- args[2]
rows <- list()
for (s in 2020:2026) {
  pt <- read.csv(file.path(root, paste0("seed", s), "pathway_truth.csv"))
  st <- read.csv(file.path(out, paste0("seed", s), "status.csv"))
  fg <- do.call(rbind, lapply(pt$pathway, function(p) {
    f <- readRDS(file.path(out, paste0("seed", s), p, "fgsea.rds"))
    data.frame(ES = f$ES, NES = f$NES, pval = f$pval, size = f$size, n_leading_edge = length(f$leadingEdge[[1]]))
  }))
  d <- data.frame(seed = s, pathway = pt$pathway, spatial_class = pt$spatial_class, effect_profile = pt$effect_profile,
                  response_group = pt$response_group, K = pt$K, n_up = pt$n_up, n_down = pt$n_down, fg)
  d$q <- p.adjust(d$pval, "BH")
  d$detected <- d$q < .05
  d$direction <- ifelse(d$NES > 0, "up", "down")
  d$direction_correct <- ifelse(d$K > 0, d$direction == ifelse(d$n_up >= d$n_down, "up", "down"), NA)
  d$de_seconds <- st[[grep("^(findmarkers|auc)_seconds$", names(st), value = TRUE)]][match(d$pathway, st$pathway)]
  rows[[length(rows) + 1]] <- d
}
all <- do.call(rbind, rows)
stopifnot(nrow(all) == 350)
write.csv(all, file.path(out, "fgsea_per_pathway.csv"), row.names = FALSE)
grp <- c("K0", "K1_10", "K11_30", "K31_60", "K61_100")
summ <- do.call(rbind, lapply(grp, function(g) do.call(rbind, lapply(c("broad_weak", "narrow_strong"), function(e) {
  x <- all[all$response_group == g & all$effect_profile == e, ]
  q <- quantile(x$NES, c(.25, .5, .75))
  data.frame(K_group = g, effect = e, n = nrow(x), detected = sum(x$detected),
             wrong_direction_detected = sum(x$detected & x$direction_correct %in% FALSE),
             wrong_direction_any = sum(x$direction_correct %in% FALSE),
             NES_median = q[2], NES_q1 = q[1], NES_q3 = q[3])
}))))
write.csv(summ, file.path(out, "fgsea_summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
cat("DE seconds: mean", round(mean(all$de_seconds), 2), "total", round(sum(all$de_seconds) / 60, 1), "min\n")
