# PathwayLGM reproducible environment
# Build:  docker build -t pathwaylgm -f Dockerfile .
# Run:    docker run --rm -it -v ${PWD}/output:/work/output pathwaylgm
# Note: on Apple Silicon / ARM machines add --platform linux/amd64
#       (INLA ships x86_64 binaries; emulation will be slow).

FROM rocker/r-ver:4.5.3

# ---- system libraries (rhdf5, XML, magick, curl/ssl, git) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev libhdf5-dev libpng-dev libjpeg62-turbo-dev \
        libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
        libmagick++-dev libssl-dev libcurl4-openssl-dev libgfortran5 git \
    && rm -rf /var/lib/apt/lists/*

# ---- CRAN deps ----
RUN install2.r --error --skipinstalled \
        Matrix pbivnorm qrng digest jsonlite energy spdep fmesher \
        BiocManager remotes

# ---- INLA (own repository, not on CRAN) ----
RUN R -e "install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dependencies=TRUE)" \
    && R -e "library(INLA); cat('INLA', as.character(packageVersion('INLA')), 'OK\n')"

# ---- Bioconductor (DESeq2 for dispersion pre-estimation, GSVA for benchmark) ----
RUN R -e "BiocManager::install(c('DESeq2','GSVA','S4Vectors','SummarizedExperiment'), update=FALSE, ask=FALSE)" \
    && R -e "suppressMessages(library(DESeq2)); suppressMessages(library(GSVA)); cat('Bioc OK\n')"

# ---- PathwayLGM itself ----
# Assumes the build context contains the package source in ./PathwayLGM
# (with DESCRIPTION, R/, NAMESPACE, ...). Build context = parent dir of the package.
COPY PathwayLGM /opt/PathwayLGM
RUN R CMD INSTALL /opt/PathwayLGM \
    && R -e "library(PathwayLGM); cat('PathwayLGM', as.character(packageVersion('PathwayLGM')), 'installed OK\n')"

# If R CMD INSTALL fails at "lazy loading" (as seen on some Windows builds),
# replace the line above with:
#   RUN R CMD INSTALL --no-byte-compile /opt/PathwayLGM
# or source the R files directly at runtime instead of installing.

WORKDIR /work
CMD ["R"]
