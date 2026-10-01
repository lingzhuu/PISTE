# Installation

SpaPath is a pure R package. It needs R 4.1.0 or later.

## 1. Install INLA

INLA is not on CRAN and is installed from its own repository:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dependencies = TRUE
)
```

Check that it loads:

```r
library(INLA)
packageVersion("INLA")
```

See the [INLA download page](https://www.r-inla.org/download-install) if the
installation fails on your platform.

## 2. Install the other dependencies

```r
install.packages(c("Matrix", "qrng", "fmesher", "remotes"))
```

`Matrix` and `qrng` are required. `fmesher` is suggested and is the usual way
to build the mesh for the spatial fields. `remotes` is only needed for the
GitHub installation below.

## 3. Install SpaPath

The package source lives in `src/spapath/` on the `dev` branch of
[lingzhuu/pathwayLGM](https://github.com/lingzhuu/pathwayLGM).

### From GitHub

```r
remotes::install_github("lingzhuu/pathwayLGM", ref = "dev",
                        subdir = "src/spapath")
```

If the repository is private, create a GitHub personal access token with read
access and pass it as `auth_token = "<token>"`, or set the `GITHUB_PAT`
environment variable.

### From a local clone

```bash
git clone -b dev https://github.com/lingzhuu/pathwayLGM.git
R CMD INSTALL pathwayLGM/src/spapath
```

If installation stops at the lazy-loading step on some Windows systems, try
`R CMD INSTALL --no-byte-compile pathwayLGM/src/spapath`.

## 4. Check the installation

```r
library(PathwayLGM)
packageVersion("PathwayLGM")
```

The exported functions are `screen_joint_fit()`, `screen_summarize()`,
`screen_write()`, `pathway_pi()` and `pathway_pi_integrate()`.
