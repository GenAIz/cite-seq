# syntax=docker/dockerfile:1
FROM rocker/tidyverse:4.5.1

ARG BPCELLS_REF=v0.3.1

ENV CRAN_REPO=https://packagemanager.posit.co/cran/__linux__/jammy/2025-11-01
ENV BIOC_VERSION=3.21
ENV BIOC_SOFT=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/bioc
ENV BIOC_ANN=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/data/annotation
ENV BIOC_EXP=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/data/experiment
ENV BIOC_WORK=https://packagemanager.posit.co/bioconductor/__linux__/jammy/2025-11-01/workflows

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    curl wget git \
    patch \
    libhdf5-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev

RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MINIFORGE=Miniforge3-Linux-x86_64.sh ; \
    else \
        MINIFORGE=Miniforge3-Linux-aarch64.sh ; \
    fi && \
    curl -L https://github.com/conda-forge/miniforge/releases/latest/download/$MINIFORGE \
         -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"

COPY environment.yml /tmp/environment.yml
RUN --mount=type=cache,target=/opt/conda/pkgs set -eu; \
    conda env create -f /tmp/environment.yml; \
    rm /tmp/environment.yml; \
    conda clean --tarballs --index-cache -y || true; \
    find /opt/conda/envs -name "*.pyc" -delete 2>/dev/null || true; \
    find /opt/conda/envs -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true; \
    test -x /opt/conda/envs/scanpy-env/bin/python \
      || { echo "FATAL: conda env 'scanpy-env' was not created"; exit 1; }

RUN /opt/conda/envs/scanpy-env/bin/python - <<'PYCHK'
import importlib, sys
mods = ["numpy", "h5py", "anndata", "mudata", "scanpy", "muon"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:
        bad.append(f"{m}: {type(e).__name__}: {e}")
if bad:
    sys.exit("FATAL: python imports failed in scanpy-env:\n  " + "\n  ".join(bad))
import numpy
if int(numpy.__version__.split(".")[0]) >= 2:
    sys.exit(f"FATAL: numpy {numpy.__version__} - environment.yml pins 1.26.4; the solver drifted")
print("python env OK:", ", ".join(mods))
PYCHK

RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    options(Ncpus=8); \
    pkgs <- c('remotes','optparse','jsonlite','reticulate','Matrix'); \
    install.packages(pkgs, repos=Sys.getenv('CRAN_REPO')); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('CRAN install failed: ', paste(missing, collapse=', '))" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

RUN R -e " \
    Sys.setenv(RETICULATE_PYTHON='/opt/conda/envs/scanpy-env/bin/python'); \
    reticulate::py_run_string('import h5py, anndata'); \
    cat('reticulate in-process h5py/anndata import OK\n')"

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
    for (i in 1:3) { \
      try(BiocManager::install(version=Sys.getenv('BIOC_VERSION'), ask=FALSE, update=TRUE)); \
      pin_ok <- requireNamespace('BiocVersion', quietly=TRUE) && \
        startsWith(as.character(utils::packageVersion('BiocVersion')), Sys.getenv('BIOC_VERSION')); \
      if (pin_ok) break; \
      Sys.sleep(15); \
    }; \
    if (!pin_ok) stop('Bioconductor version pin FAILED after 3 attempts: BiocVersion ', \
      Sys.getenv('BIOC_VERSION'), ' could not be installed (transient network?). ', \
      'Without this pin BiocManager silently drifts to the CURRENT Bioconductor ', \
      'release for the rest of the layer, pulling unpinned package versions from ', \
      'live bioconductor.org instead of the dated Posit snapshot.'); \
    pkgs <- c('HDF5Array','AnnotationDbi'); \
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

RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    options(Ncpus=8); \
    options(repos = c( \
        CRAN = Sys.getenv('CRAN_REPO'), \
        BioCsoft = Sys.getenv('BIOC_SOFT'), \
        BioCann  = Sys.getenv('BIOC_ANN'), \
        BioCexp  = Sys.getenv('BIOC_EXP'), \
        BioCworkflows = Sys.getenv('BIOC_WORK') \
    )); \
    pkgs <- c('org.Hs.eg.db','org.Mm.eg.db'); \
    BiocManager::install(pkgs, ask=FALSE, update=FALSE); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('Annotation package install failed: ', paste(missing, collapse=', '))" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

ENV BPCELLS_REF=${BPCELLS_REF}
RUN --mount=type=cache,target=/root/.cache/R \
    R -e " \
    install.packages('hdf5r', repos=Sys.getenv('CRAN_REPO')); \
    remotes::install_github(paste0('bnprks/BPCells/r@', Sys.getenv('BPCELLS_REF')), upgrade='never'); \
    pkgs <- c('hdf5r','BPCells'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('install failed: ', paste(missing, collapse=', ')); \
    cat('BPCells', as.character(utils::packageVersion('BPCells')), '\n')" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

ENV RETICULATE_PYTHON=/opt/conda/envs/scanpy-env/bin/python
ENV PIPELINE_PY_ENV=scanpy-env
ENV RETICULATE_USE_MANAGED_VENV=no

RUN R -q -e " \
    py <- reticulate::py_config()\$python; \
    cat('reticulate python:', py, '\n'); \
    if (!startsWith(py, '/opt/conda/envs/scanpy-env/')) \
      stop('reticulate bound ', py, ' instead of the pinned scanpy-env interpreter'); \
    v <- reticulate::py_eval('__import__(\"importlib.metadata\", fromlist=[\"x\"]).version(\"anndata\")'); \
    cat('anndata:', v, '\n'); \
    if (utils::compareVersion(v, '0.12') >= 0) \
      stop('anndata ', v, ' writes obs/_index as an HDF5 group; BPCells cannot read it'); \
    cat('OK: reticulate pinned, anndata ', v, ' is BPCells-compatible\n', sep='')"

WORKDIR /app

COPY main.R .
COPY error_codes.R .
COPY portability.R .
COPY disk_backed.R .
COPY load_data.R .
COPY entrypoint.sh .
COPY configs/ configs/

RUN chmod +x entrypoint.sh

ENTRYPOINT ["/bin/bash", "/app/entrypoint.sh"]
