# Adjacent sections on DLPFC

<!-- Generated from doc/examples/joint_sections.Rmd by doc/examples/render.R; edit the .Rmd. -->

This example tests pathways in white matter (WM) against the cortical layers
in two adjacent DLPFC Visium sections, 151673 and 151674, from the same donor
(spatialLIBD; Maynard et al., 2021). The two sections are fitted jointly, so
that the score reflects the regional difference reproduced in both.

## Requirements

This example also uses:

| Package | Source | Used for |
|---------|--------|----------|
| spatialLIBD | Bioconductor | Provides the DLPFC Visium data and its layer annotation |
| sctransform | CRAN | Normalizes the counts |

## 1. Load the data

``` r
library(spatialLIBD)
library(SpatialExperiment)
library(PISTE)

spe <- fetch_data(type = "spe")
spe <- spe[, spe$sample_id %in% c("151673", "151674")]
spe <- spe[, !is.na(spe$layer_guess_reordered)]

counts <- t(counts(spe))                       # spots x genes, raw counts
colnames(counts) <- make.unique(rowData(spe)$gene_name)
rownames(counts) <- paste(spe$sample_id, colnames(spe), sep = "_")
region  <- as.character(spe$layer_guess_reordered)
section <- spe$sample_id
table(section, region)
```

## 2. Normalize

Each section is normalized on its own:

``` r
normalize <- function(x) {
  vst_out <- sctransform::vst(t(x), vst.flavor = "v2", verbosity = 0,
                              res_clip_range = c(-1, 1) * sqrt(nrow(x) / 30))
  t(vst_out$y)
}
per_section <- lapply(split(seq_len(nrow(counts)), section),
                      function(rows) normalize(counts[rows, ]))
genes <- Reduce(intersect, lapply(per_section, colnames))
expression <- do.call(rbind, lapply(per_section, function(x) x[, genes]))
expression <- expression[rownames(counts), ]   # spot order of region and section
```

## 3. Gene sets

``` r
gmt <- "benchmark/analysis/DLPFC/data/DLPFC_brain_and_reference.gmt"
```

## 4. Run PISTE

``` r
res <- compare(expression, region, gmt, baseline = "WM",
               section = section, joint_sections = TRUE,
               out = "results_joint")
```

## 5. Read the results

``` r
wm <- res$overall[res$overall$region == "WM", ]
wm <- wm[order(wm$evidence, decreasing = TRUE), ]
rownames(wm) <- NULL
wm[, c("gene_set", "genes", "score", "evidence")]
```

`score` and `evidence` are read as in the
[observed-mode example](observed.md#5-read-the-results), here for the
difference shared by the two sections.
