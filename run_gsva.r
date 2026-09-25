library(GSVA)

x <- readRDS("../data/DLPFC/processed_data.rds")
gmt <- strsplit(readLines("../data/pathway_databases/h.all.v2026.1.Hs.symbols.gmt"), "\t")
# Keep gene IDs after the first two fields; name each gene list using the first field.
pathways <- setNames(lapply(gmt, function(z) z[-c(1, 2)]),
                     sapply(gmt, function(z) z[1]))

scores <- gsva(gsvaParam(x, pathways, kcdf = "auto"))

# Row: spots, Column: pathways scores
write.csv(t(scores), "../GSVA/DLPFC_gsva_scores.csv")