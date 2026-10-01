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

## If INLA fails to install

If the command above stops while
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

SpaPath can also be installed using docker image.

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

