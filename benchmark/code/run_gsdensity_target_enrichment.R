# Chi-square enrichment test of high-activity spots in the target region.
# Per pathway, a 2x2 table (target/rest x positive/negative label) is tested
# with chisq.test(correct = TRUE), the Pearson chi-square with Yates
# continuity correction (R default for 2x2 tables), then BH-adjusted across
# pathways. Constant labels make the test unavailable and give NA, not a
# fabricated p. Saves gsdensity_target_enrichment_scores.csv, one row per
# pathway.

# ---- config ---------------------------------------------------------------------
work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")

# ---- load -----------------------------------------------------------------------
x <- readRDS(file.path(data_dir, "raw_data.rds"))
labels <- read.csv(file.path(data_dir, "gsdensity_PAL_binary_labels.csv"),
                   row.names = 1, check.names = FALSE)

# ---- run ------------------------------------------------------------------------
rows <- list()
for (p in names(labels)) {
  target <- as.logical(x$oracle_target[, p])
  tab <- table(factor(ifelse(target, "target", "rest"), levels = c("target", "rest")),
               factor(labels[[p]], levels = c("positive", "negative")))
  row <- data.frame(pathway = p, status = "ok",
    target_positive = tab["target", "positive"], target_negative = tab["target", "negative"],
    rest_positive = tab["rest", "positive"], rest_negative = tab["rest", "negative"],
    target_positive_fraction = tab["target", "positive"] / sum(tab["target", ]),
    rest_positive_fraction = tab["rest", "positive"] / sum(tab["rest", ]),
    odds_ratio = NA_real_, chi_square = NA_real_, df = NA_integer_,
    p_value = NA_real_, min_expected_count = NA_real_)
  row$target_minus_rest_fraction <- row$target_positive_fraction - row$rest_positive_fraction
  if (all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
    test <- chisq.test(tab, correct = TRUE)
    row$odds_ratio <- (tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1])
    row$chi_square <- unname(test$statistic)
    row$df <- unname(test$parameter)
    row$p_value <- test$p.value
    row$min_expected_count <- min(test$expected)
  } else {
    row$status <- "constant_labels"
  }
  rows[[p]] <- row
  cat(p, row$status, "p =", row$p_value, "\n")
}

# ---- save -----------------------------------------------------------------------
scores <- do.call(rbind, rows)
scores$q_BH <- p.adjust(scores$p_value, method = "BH")
write.csv(scores, file.path(data_dir, "gsdensity_target_enrichment_scores.csv"),
          row.names = FALSE)
