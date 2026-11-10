# syntax=docker/dockerfile:1
FROM rocker/tidyverse:4.5.1

ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/2025-11-01

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
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
    pkgs <- c('optparse','jsonlite','readxl','dplyr','tibble'); \
    install.packages(pkgs, repos=Sys.getenv('CRAN_REPO')); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('CRAN install failed: ', paste(missing, collapse=', '))" \
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
COPY add_metadata.R .
COPY entrypoint.sh .
COPY configs/ configs/

RUN chmod +x entrypoint.sh

ENTRYPOINT ["/bin/bash", "/app/entrypoint.sh"]
