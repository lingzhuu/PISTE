# Installation

SpaPath requires R 4.1.0 or later. Install it from GitHub with:

```r
install.packages("remotes")
remotes::install_github("lingzhuu/pathwayLGM", ref = "dev",
                        subdir = "src/spapath", dependencies = TRUE)
```

This installs SpaPath together with its dependencies, including INLA from the
INLA repository and fmesher for building meshes.

Check the installation:

```r
library(PathwayLGM)
```

## If INLA fails to install

INLA is a large package that is not on CRAN. If the command above stops while
installing INLA, install INLA on its own first and then run the command again:

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

## Install from source

```bash
git clone -b dev https://github.com/lingzhuu/pathwayLGM.git
R CMD INSTALL pathwayLGM/src/spapath
```

When installing from source, install the dependencies (INLA, Matrix, qrng and
fmesher) first.
