# Installation

Requires R ≥ 4.1.0.

```r
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dependencies = TRUE)
install.packages("remotes")
remotes::install_github("lingzhuu/PISTE", subdir = "src/piste",
                        dependencies = TRUE)
```

## Docker

```bash
git clone https://github.com/lingzhuu/PISTE.git
cd PISTE
docker build -t piste .
docker run --rm -it -v "$PWD":/work piste
```

## Install from source

```bash
git clone https://github.com/lingzhuu/PISTE.git
R CMD INSTALL PISTE/src/piste
```
