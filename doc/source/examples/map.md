# Spatial map on breast cancer

<!-- Generated from doc/examples/map.Rmd by doc/examples/render.R; edit the .Rmd. -->

This example maps where T-cell exhaustion and the interferon-γ response are
more or less active in one triple-negative breast cancer Visium section
(GSM6433597; Bassiouni et al., 2023), processed by Zhu et al. (2025). No region labels are used.

## Requirements

This example also uses:

| Package | Source | Used for |
|---------|--------|----------|
| ggplot2 | CRAN | Plots the maps |

## 1. Load the data

Download `GSM_giotto_processed.tar.gz` from
[Zenodo 15247064](https://zenodo.org/records/15247064) and extract it. The
expression in `Giotto_norm_expr.txt` was already normalized by the original
authors with Giotto, so it is used as is.

``` r
library(PISTE)

sample_dir <- "GSM_giotto_processed/dir_GSM6433597"
expression <- t(as.matrix(read.delim(file.path(sample_dir, "Giotto_norm_expr.txt"),
                                     row.names = 1, check.names = FALSE)))
locs   <- read.delim(file.path(sample_dir, "Giotto_spatial_locs.txt"))
coords <- as.matrix(locs[match(rownames(expression), locs$cell_ID), c("sdimx", "sdimy")])
dim(expression)
```

## 2. Gene sets

``` r
gene_sets <- list(
  T_cell_exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX",
                        "ENTPD1", "CXCL13", "LAYN"),
  IFNG_response     = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA")
)
```

## 3. Run PISTE

``` r
res <- map_activity(expression, coords, gene_sets, out = "results_map")
res$gene_sets[, c("gene_set", "genes_used", "sd_z", "spots_higher", "spots_lower")]
```

## 4. Plot the maps

`z_mean` is the posterior mean of the pathway field at each spot; `call` marks
the regions of higher and lower activity.

``` r
library(ggplot2)

spots <- do.call(rbind, lapply(names(res$results), function(set) {
  data.frame(gene_set = set, res$results[[set]]$spots,
             x = coords[, 1], y = coords[, 2])
}))

ggplot(spots, aes(x, y, colour = z_mean)) +
  geom_point(size = 0.6) +
  scale_colour_gradient(low = rgb(247, 224, 221, maxColorValue = 255),
                        high = rgb(210, 59, 43, maxColorValue = 255)) +
  facet_wrap(~ gene_set) + coord_equal() + theme_void()

ggplot(spots, aes(x, y, colour = call)) +
  geom_point(size = 0.6) +
  scale_colour_manual(values = c(higher = rgb(210, 59, 43, maxColorValue = 255),
                                 lower = rgb(83, 104, 166, maxColorValue = 255),
                                 none = "grey85")) +
  facet_wrap(~ gene_set) + coord_equal() + theme_void()
```
