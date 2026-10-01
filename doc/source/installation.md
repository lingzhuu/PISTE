# Installation

SpaPath requires R 4.1.0 or later. Install it from GitHub with:

```r
install.packages("remotes")
remotes::install_github("lingzhuu/pathwayLGM", subdir = "src/spapath",
                        dependencies = TRUE)
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

## Docker

The Docker image contains R, INLA and SpaPath, so nothing needs to be
installed by hand.

```bash
git clone https://github.com/lingzhuu/pathwayLGM.git
cd pathwayLGM
docker build -t spapath .
docker run --rm -it -v "$PWD":/work spapath
```

The last command opens R with the current folder available at `/work`. On
Apple Silicon or other ARM machines, add `--platform linux/amd64` to the build
and run commands.

## Install from source

```bash
git clone https://github.com/lingzhuu/pathwayLGM.git
R CMD INSTALL pathwayLGM/src/spapath
```

When installing from source, install the dependencies (INLA, Matrix, qrng and
fmesher) first.
