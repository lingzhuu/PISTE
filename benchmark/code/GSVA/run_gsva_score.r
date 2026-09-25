library(GSVA)

work_dir <- "../../analysis/DLPFC"
data_dir <- file.path(work_dir, "data")
out_dir <- file.path(work_dir, "GSVA")

x <- readRDS(file.path(data_dir, "processed_data.rds"))

gmt <- strsplit(readLines(file.path(data_dir, "h.all.v2026.1.Hs.symbols.gmt")),"\t")

# Keep gene IDs after the first two fields; name each gene list using the first field.
pathways <- setNames(
  lapply(gmt, function(z) z[-c(1, 2)]),
  sapply(gmt, function(z) z[1])
)

scores <- gsva(gsvaParam(x, pathways, kcdf = "auto"))

# Save result (Row: spots, Column: pathways scores)
write.csv(t(scores), file.path(out_dir, "gsva_scores.csv"))
