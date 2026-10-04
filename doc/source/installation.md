# Installation

SpaPath requires R 4.1.0 or later. Install it from GitHub with:

```r
install.packages("remotes")
remotes::install_github("lingzhuu/pathwayLGM", subdir = "src/spapath",
                        dependencies = TRUE)
```

Check the installation:

```r
library(PathwayLGM)
```

## Dependencies

The command above installs SpaPath together with all the packages it depends
on:

| Package | Source | Used for |
|---------|--------|----------|
| INLA | [INLA repository](https://www.r-inla.org) | Model fitting |
| fmesher | CRAN | Building the mesh for the spatial fields |
| Matrix | CRAN | Sparse matrices |
| qrng | CRAN | Quadrature for the posterior evidence |
| DESeq2 | Bioconductor | Dispersion estimation in the composition-adjusted mode |
| S4Vectors, SummarizedExperiment | Bioconductor | Required by the DESeq2 step |
| sctransform | CRAN | SCTransform normalization in `run_spapath()` and the [command line](cli.md) |
| hdf5r | CRAN | Reading `.h5ad` and 10x `.h5` files in the command line |
| glmGamPoi | Bioconductor | Faster SCTransform (optional) |

Packages that are only needed to run the [examples](examples/observed.md), such
as the data package and the deconvolution method, are listed at the start of
each example.

## If INLA fails to install

If the command above stops while installing INLA, install INLA on its own
first and then run the command again:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dependencies = TRUE
)
```

See the [INLA website](https://www.r-inla.org/download-install) for
platform-specific notes.

## Docker

SpaPath can also be installed using the Docker image, which contains R,
SpaPath, all its dependencies, the `spapath` [command](cli.md) and the
packages used in the examples.

```bash
git clone https://github.com/lingzhuu/pathwayLGM.git
cd pathwayLGM
docker build -t spapath .
docker run --rm -it -v "$PWD":/work spapath
```

## Install from source

```bash
git clone https://github.com/lingzhuu/pathwayLGM.git
R CMD INSTALL pathwayLGM/src/spapath
```

`R CMD INSTALL` does not install dependencies; install the packages in the
table above first.
