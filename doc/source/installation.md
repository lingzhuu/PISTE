# Installation

PISTE needs R 4.1 or later and INLA. Install INLA first, then PISTE:

```r
inla_repo <- "https://inla.r-inla-download.org/R/stable"
install.packages("INLA", repos = c(getOption("repos"), INLA = inla_repo),
                 dependencies = TRUE)
install.packages("remotes")
remotes::install_github("lingzhuu/PISTE", subdir = "src/piste",
                        dependencies = TRUE)
```

Or build the Docker image, which has R, INLA, PISTE and the `piste` command:

```bash
git clone https://github.com/lingzhuu/PISTE.git
cd PISTE
docker build -t piste .
docker run --rm -it -v "$PWD":/work piste
```

Or install from a clone of the repository:

```bash
git clone https://github.com/lingzhuu/PISTE.git
R CMD INSTALL PISTE/src/piste
```
