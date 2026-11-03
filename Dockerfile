# PISTE with INLA preinstalled.
#
# Build from the repository root:
#   docker build -t piste .
# Run the command-line interface on files in the current folder:
#   docker run --rm -v "$PWD":/work piste piste run --help
# Or open R with the current folder mounted at /work:
#   docker run --rm -it -v "$PWD":/work piste
#
# INLA ships x86_64 Linux binaries; on Apple Silicon / ARM machines add
# --platform linux/amd64 to both commands (runs under emulation).

FROM rocker/r-ver:4.5.3

# System libraries needed by INLA and by packages built from source.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgfortran5 libglpk40 libgsl27 libxml2-dev libssl-dev \
        libcurl4-openssl-dev libudunits2-dev libgdal-dev libgeos-dev \
        libproj-dev libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN dependencies of PISTE and INLA.
RUN install2.r --error --skipinstalled \
        Matrix qrng fmesher sp foreach BiocManager remotes sctransform hdf5r

# Bioconductor dependencies of PISTE (composition-adjusted mode).
RUN R -q -e "BiocManager::install(c('DESeq2', 'S4Vectors', 'SummarizedExperiment', 'glmGamPoi'), update = FALSE, ask = FALSE)"

# INLA, from its own repository (not on CRAN).
RUN R -q -e "install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'))" \
    && R -q -e "library(INLA); cat('INLA', as.character(packageVersion('INLA')), 'OK\n')"

# PISTE itself.
COPY src/piste /opt/piste
RUN R CMD INSTALL /opt/piste \
    && R -q -e "library(PISTE); cat('PISTE', as.character(packageVersion('PISTE')), 'OK\n')" \
    && ln -s "$(Rscript -e "cat(system.file('bin/piste', package = 'PISTE'))")" /usr/local/bin/piste \
    && piste --version

# Packages used in the documentation examples: DLPFC data (spatialLIBD) and
# RCTD deconvolution (spacexr).
RUN R -q -e "BiocManager::install(c('spatialLIBD', 'HDF5Array'), update = FALSE, ask = FALSE)" \
    && R -q -e "remotes::install_github('dmcable/spacexr', build_vignettes = FALSE, upgrade = 'never')"

WORKDIR /work
CMD ["R"]
