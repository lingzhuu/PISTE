# SpaPath (R package PathwayLGM) with INLA preinstalled.
#
# Build from the repository root:
#   docker build -t spapath .
# Run the command-line interface on files in the current folder:
#   docker run --rm -v "$PWD":/work spapath spapath run --help
# Or open R with the current folder mounted at /work:
#   docker run --rm -it -v "$PWD":/work spapath
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

# CRAN dependencies of SpaPath and INLA.
RUN install2.r --error --skipinstalled \
        Matrix qrng fmesher sp foreach BiocManager remotes sctransform hdf5r

# Bioconductor dependencies of SpaPath (composition-adjusted mode).
RUN R -q -e "BiocManager::install(c('DESeq2', 'S4Vectors', 'SummarizedExperiment', 'glmGamPoi'), update = FALSE, ask = FALSE)"

# INLA, from its own repository (not on CRAN).
RUN R -q -e "install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'))" \
    && R -q -e "library(INLA); cat('INLA', as.character(packageVersion('INLA')), 'OK\n')"

# SpaPath itself.
COPY src/spapath /opt/spapath
RUN R CMD INSTALL /opt/spapath \
    && R -q -e "library(PathwayLGM); cat('PathwayLGM', as.character(packageVersion('PathwayLGM')), 'OK\n')" \
    && ln -s "$(Rscript -e "cat(system.file('bin/spapath', package = 'PathwayLGM'))")" /usr/local/bin/spapath \
    && spapath --version

# Packages used in the documentation examples: DLPFC data (spatialLIBD),
# SCTransform (Seurat) and RCTD deconvolution (spacexr).
RUN install2.r --error --skipinstalled Seurat \
    && R -q -e "BiocManager::install(c('spatialLIBD', 'HDF5Array'), update = FALSE, ask = FALSE)" \
    && R -q -e "remotes::install_github('dmcable/spacexr', build_vignettes = FALSE, upgrade = 'never')"

WORKDIR /work
CMD ["R"]
