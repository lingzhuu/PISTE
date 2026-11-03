# PISTE

**P**athway **I**nference from **S**patial **T**ranscriptomic **E**xpression

PISTE is an R package that tests whether a gene pathway is more active in one
tissue region than in another in spatial transcriptomics data. It models all
genes of a gene set jointly, accounts for the spatial dependence between
neighbouring spots with gene-specific spatial random fields, and fits the
model with INLA. For every pair of regions, and for each region against the
rest of the tissue, PISTE reports a directional score π between 0 and 1 and
its posterior evidence.

PISTE can be run on normalized expression (observed mode) or on raw counts
adjusted for cell-type composition (composition-adjusted mode).

![Overview of PISTE](figure/method_overview.png)

## Installation

PISTE requires R 4.1.0 or later.

```r
install.packages("remotes")
remotes::install_github("lingzhuu/pathwayLGM", subdir = "src/piste",
                        dependencies = TRUE)
```

If INLA fails to install, install it first and run the command again:

```r
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dependencies = TRUE)
```

A Docker image with R, INLA, PISTE and the `piste` command can be built from
the repository:

```bash
git clone https://github.com/lingzhuu/pathwayLGM.git
cd pathwayLGM
docker build -t piste .
```

## Quick start

### R

```r
library(PISTE)

# counts: spots x genes raw counts; xy: spots x 2 spot coordinates;
# region: region label of each spot; gene_sets: named list or GMT file
res <- run_piste(counts, xy, region, "pathways.gmt", baseline = "WM",
                 out = "results")
head(res$overall)
```

### Command line

```bash
piste run --input sample.h5ad --region-key layer --target WM \
          --gmt pathways.gmt --out results
```

## Citation

A manuscript describing PISTE is in preparation.

