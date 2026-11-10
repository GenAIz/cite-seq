# syntax=docker/dockerfile:1
FROM rocker/tidyverse:4.5.1

ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/2025-11-01
ENV BIOC_VERSION=3.21
ENV BIOC_SOFT=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/bioc
ENV BIOC_ANN=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/data/annotation
ENV BIOC_EXP=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/data/experiment
ENV BIOC_WORK=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/workflows

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    curl wget \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    libhdf5-dev

RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    options(Ncpus=8); \
    pkgs <- c('optparse','jsonlite','Matrix'); \
    install.packages(pkgs, repos=Sys.getenv('CRAN_REPO')); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('CRAN install failed: ', paste(missing, collapse=', '))" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    options(Ncpus=8); \
    install.packages('BiocManager', repos=Sys.getenv('CRAN_REPO')); \
    options(repos = c( \
        CRAN = Sys.getenv('CRAN_REPO'), \
        BioCsoft = Sys.getenv('BIOC_SOFT'), \
        BioCann  = Sys.getenv('BIOC_ANN'), \
        BioCexp  = Sys.getenv('BIOC_EXP'), \
        BioCworkflows = Sys.getenv('BIOC_WORK') \
    )); \
    BiocManager::install(version=Sys.getenv('BIOC_VERSION'), ask=FALSE, update=TRUE); \
    pkgs <- c('UCell','BiocParallel'); \
    BiocManager::install(pkgs, ask=FALSE, update=FALSE); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('Bioconductor install failed: ', paste(missing, collapse=', '))" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    options(Ncpus=8); \
    install.packages('Seurat', repos=Sys.getenv('CRAN_REPO')); \
    stopifnot(requireNamespace('Seurat', quietly=TRUE))" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

ARG BPCELLS_REF=v0.3.1
ENV BPCELLS_REF=${BPCELLS_REF}
RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos=Sys.getenv('CRAN_REPO')); \
    remotes::install_github(paste0('bnprks/BPCells/r@', Sys.getenv('BPCELLS_REF')), upgrade='never'); \
    stopifnot(requireNamespace('BPCells', quietly=TRUE)); \
    cat('BPCells', as.character(utils::packageVersion('BPCells')), '\n')" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

WORKDIR /app

COPY main.R .
COPY error_codes.R .
COPY portability.R .
COPY disk_backed.R .
COPY classify_functional_state.R .
COPY cell_state_markers.json .
COPY cell_state_markers_adt.json .
COPY entrypoint.sh .
COPY configs/ configs/

RUN chmod +x entrypoint.sh

ENTRYPOINT ["/bin/bash", "/app/entrypoint.sh"]
