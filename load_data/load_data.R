library(Seurat)
library(Matrix)

#' Convert Ensembl gene IDs to symbols on a matrix's rownames
#'
#' Uses AnnotationDbi with org.Hs.eg.db (human) or org.Mm.eg.db (mouse).
#' Gene names that are already symbols are passed through unchanged, and the
#' matrix is returned as-is when conversion is not possible.
#'
#' @param mat Matrix whose rownames are gene identifiers.
#' @param species `"human"` or `"mouse"`.
#' @param verbose Print what was converted or skipped.
#' @return The matrix with converted rownames (unchanged on passthrough).
.convert_gene_ids <- function(mat, species = "human", verbose = FALSE) {
  genes <- rownames(mat)

  if (is.null(genes) || length(genes) == 0) {
    if (verbose) message("No gene names found - skipping conversion.")
    return(mat)
  }

  idx <- seq_len(min(100, length(genes)))
  is_ensembl <- all(grepl("^ENS[A-Z]*[0-9]+", genes[idx]))

  if (!is_ensembl) {
    if (verbose) message("Gene IDs appear to be symbols already - skipping conversion.")
    return(mat)
  }

  if (verbose) message("Detected Ensembl gene IDs - converting to gene symbols using AnnotationDbi...")

  db_pkg <- if (tolower(species) %in% c("human", "homo sapiens", "hs")) {
    "org.Hs.eg.db"
  } else if (tolower(species) %in% c("mouse", "mus musculus", "mm")) {
    "org.Mm.eg.db"
  } else {
    stop("Unsupported species for gene conversion: ", species)
  }

  if (!requireNamespace(db_pkg, quietly = TRUE))
    stop(sprintf("Required package %s is not installed.", db_pkg))
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Required package AnnotationDbi is not installed.")

  db <- get(db_pkg, envir = asNamespace(db_pkg))
  syms <- AnnotationDbi::mapIds(
    db,
    keys = genes,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  syms[is.na(syms)] <- genes[is.na(syms)]

  if (anyDuplicated(syms)) {
    if (verbose) message("Duplicate gene symbols detected - making unique...")
    syms <- make.unique(syms)
  }

  rownames(mat) <- syms

  if (verbose) {
    n_converted <- sum(syms != genes)
    message(sprintf("Gene ID conversion complete: %d / %d genes mapped.",
                    n_converted, length(genes)))
  }

  return(mat)
}

#' Bind reticulate to the intended Python before anything else does
#'
#' Resolution order: (1) `RETICULATE_PYTHON` - set by the Dockerfile to the
#' bundled conda env, always the answer inside a container; (2) the
#' `PIPELINE_PY_ENV` conda env name, for runs outside the image; (3) neither
#' set - leave reticulate to its own discovery, whose caller-side import
#' probe reports problems with a clearer message.
#'
#' @param verbose Print the interpreter that was bound.
#' @return Invisibly `NULL`.
.ensure_python <- function(verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(NULL))

  rp <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(rp) && file.exists(rp)) {
    if (verbose) message(" Python interpreter: ", rp, " (RETICULATE_PYTHON)")
    return(invisible(rp))
  }

  env_name <- Sys.getenv("PIPELINE_PY_ENV", "")
  if (!nzchar(env_name)) {
    if (verbose)
      message(" RETICULATE_PYTHON / PIPELINE_PY_ENV unset - using reticulate's default interpreter.")
    return(invisible(NULL))
  }

  conda_path <- Sys.which("conda")
  if (!nzchar(conda_path))
    conda_path <- Find(file.exists, path.expand(c(
      "/opt/conda/bin/conda",
      "~/miniforge3/bin/conda", "~/mambaforge/bin/conda",
      "~/miniconda3/bin/conda", "~/anaconda3/bin/conda")))

  py <- NULL
  if (!is.null(conda_path) && nzchar(conda_path)) {
    py <- tryCatch(reticulate::conda_python(env_name, conda = conda_path),
                   error = function(e) NULL)
    if (is.null(py) || !file.exists(py)) {
      cand <- file.path(dirname(dirname(conda_path)), "envs", env_name, "bin", "python")
      py <- if (file.exists(cand)) cand else NULL
    }
  }

  if (!is.null(py) && file.exists(py)) {
    Sys.setenv(RETICULATE_PYTHON = py)
    if (verbose) message(" Python interpreter: ", py, " (conda env '", env_name, "')")
    return(invisible(py))
  }

  if (verbose)
    message(" Could not resolve conda env '", env_name,
            "' - using reticulate's default interpreter.")
  invisible(NULL)
}

#' Keys of a Python mapping, read safely through reticulate
#'
#' With `convert = FALSE`, `names()` on a Python mapping yields the object's
#' METHOD names (`keys`, `items`, `copy`, ...) - never its keys - so a naive
#' loop fabricates entries out of Python attributes with no error raised.
#' Every obsm lookup in this file routes through this helper, which goes
#' through `builtins.list()` on the mapping's own `keys()` view.
#'
#' @param obj A non-converted reticulate reference to a Python mapping.
#' @param verbose Print the failure reason when keys cannot be listed.
#' @return Character vector of keys; empty on any failure.
.py_mapping_keys <- function(obj, verbose = FALSE) {
  tryCatch({
    bi <- reticulate::import_builtins(convert = FALSE)
    as.character(reticulate::py_to_r(bi$list(obj$keys())))
  }, error = function(e) {
    if (verbose) message(sprintf(" could not list mapping keys: %s",
                                 conditionMessage(e)))
    character(0)
  })
}

#' Extract obs metadata and obsm embeddings from an .h5ad
#'
#' BPCells reads only the matrix, so cell metadata and embeddings are read separately through the Python anndata backend.
#'
#' @param path Path to the `.h5ad` file.
#' @param verbose Print progress.
#' @return List with `obs` (data frame or `NULL`) and `obsm` (named list of
#' embedding matrices or `NULL`).
.read_h5ad_metadata <- function(path, verbose = FALSE) {
  result <- list(obs = NULL, obsm = NULL)

  .ensure_python(verbose)

  tryCatch({
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      if (verbose) message("reticulate not available - skipping h5ad metadata extraction.")
      return(result)
    }

    ad <- reticulate::import("anndata", convert = FALSE)
    if (verbose) message("Reading h5ad metadata via reticulate + anndata...")
    adata <- ad$read_h5ad(path)

    obs_py <- adata$obs
    obs_df <- tryCatch(reticulate::py_to_r(obs_py), error = function(e) NULL)
    if (!is.null(obs_df) && is.data.frame(obs_df) && nrow(obs_df) > 0) {
      result$obs <- obs_df
      if (verbose) message(sprintf(" obs: %d cells x %d columns", nrow(obs_df), ncol(obs_df)))
    }

    obsm_keys <- .py_mapping_keys(adata$obsm, verbose)

    if (length(obsm_keys) > 0) {
      obsm_list <- list()
      for (key in obsm_keys) {
        emb <- tryCatch({
          m <- reticulate::py_to_r(adata$obsm[[key]])
          if (!is.matrix(m)) m <- as.matrix(m)
          m
        }, error = function(e) NULL)
        if (!is.null(emb)) obsm_list[[key]] <- emb
      }
      if (length(obsm_list) > 0) {
        result$obsm <- obsm_list
        if (verbose)
          message(sprintf(" obsm: %d embeddings (%s)",
                          length(obsm_list), paste(names(obsm_list), collapse = ", ")))
      }
    }

    rm(adata); gc()

  }, error = function(e) {
    if (verbose)
      message(sprintf("Could not extract h5ad metadata: %s", e$message))
  })

  return(result)
}

#' Load and preprocess CITE-seq RNA and optionally ADT data into a Seurat object
#'
#' This function loads paired scRNA-seq and (optionally) ADT UMI count matrices
#' from multiple formats:
#' \itemize{
#' \item AnnData (.h5ad)
#' \item MuData (.h5mu) - multi-modal data format
#' \item RDS (a pre-processed Seurat object; bare matrices are rejected)
#' }
#'
#' Matrices are mapped onto Seurat layers as the file declares them: an
#' integer raw-counts layer (`layers/counts`, `layers/raw`, `raw/X`) becomes
#' `counts` and a non-integer `X` becomes `data`; no normalization or other
#' computation is performed.
#'
#' For .h5ad inputs, matrices are read using \pkg{BPCells}
#' and Ensembl gene IDs are converted to symbols using \pkg{AnnotationDbi}
#' with org.Hs.eg.db (human) or org.Mm.eg.db (mouse). Gene names that are
#' already symbols are passed through unchanged.
#'
#' For **.h5ad inputs**, the function automatically:
#' \itemize{
#' \item **Prefers raw counts** over normalized data by probing for
#' \code{layers/raw} or \code{raw/X} before falling back to \code{X}.
#' Requires \pkg{hdf5r} for the probe; falls back to \code{X} if absent.
#' \item Extracts cell metadata (\code{adata.obs}) and attaches to
#' \code{@meta.data} (via \pkg{anndata} R package or \pkg{reticulate}).
#' \item Extracts embeddings (\code{adata.obsm}) and adds as Seurat
#' dimensional reductions (UMAP, PCA, etc.).
#' \item Detects embedded ADT/protein data in three storage patterns:
#' (a) features with \code{AB_} prefix interleaved in \code{X},
#' (b) sparse matrix in \code{uns/antibody_X} (via \pkg{hdf5r}),
#' (c) dense matrix in \code{obsm/protein} or \code{obsm/ADT}.
#' }
#' Raw inputs (.csv, .csv.gz, .tar.gz, 10X .h5) are REJECTED: this module loads
#' an object that has already been through QC, normalisation and clustering.
#'
#' For **.h5mu inputs** (MuData format):
#' \itemize{
#' \item Extracts both RNA and ADT modalities automatically
#' \item Preserves all embeddings stored in the RNA modality's `obsm` slot
#' \item Uses temporary directory for intermediate files with automatic cleanup
#' \item Performs garbage collection to optimize memory usage
#' }
#'
#' For **.rds inputs**:
#' \itemize{
#' \item If a single .rds file is provided, the object is returned **as is**
#' (Seurat object, matrix, or list).
#' \item If a directory containing multiple .rds files is provided, the files
#' are loaded as Seurat objects and **merged using JoinLayers** to create
#' a single integrated object with one layer per sample.
#' \item Supports both RNA-only and RNA+ADT multi-modal data merging.
#' \item Automatically handles assay matching and layer naming.
#' }
#'
#' @param rna_path Character. Path to ONE pre-processed object
#' (.rds/.h5ad/.h5mu). One file, not a directory. An .rds must contain a
#' Seurat object; a bare matrix is rejected.
#' @param adt_path RETIRED. Supplying it is an error. Every accepted format
#' carries RNA and ADT in one object and splits them internally, so a second
#' matrix path has nothing to point at.
#' @param top_mouse_genes Integer. Number of top expressed mouse genes to retain
#' when collapsing species (default = 100).
#' @param project_name Character. Project name to assign in Seurat object
#' (default = "SeuratProject").
#' @param species_name Character. Species name for Ensembl -> symbol conversion
#' via AnnotationDbi (default = "human"; also accepts "mouse").
#' @param rds_rna_key Character. When loading from RDS directory, the name/key
#' of the RNA matrix object (default = "rna"). Ignored for single .rds files.
#' @param rds_adt_key Character. When loading from RDS directory, the name/key
#' of the ADT matrix object (default = "adt"). Ignored for single .rds files.
#' merges them using JoinLayers (default = TRUE).
#' @param output_dir Character. Directory to save output RDS file. If NULL,
#' results are returned without saving (default = NULL). Directory will be
#' created if it doesn't exist.
#' @param output_name Character. Name for output RDS file when output_dir is
#' specified (default = "loaded_data.rds"). Ignored if output_dir is NULL.
#' @param save_output Logical. If TRUE and output_dir is specified, saves the
#' result to output_dir/output_name (default = TRUE when output_dir is provided).
#' @param use_disk_backed Logical. If TRUE, rewrite the RNA counts to an on-disk
#' BPCells matrix so QC / normalization stay within RAM on very large datasets
#' (default = FALSE). ADT is always kept in memory.
#' @param bpcells_dir Character or NULL. Directory for the on-disk BPCells matrix
#' when use_disk_backed = TRUE (default = NULL -> a subdirectory of tempdir()).
#' @param slim_to_counts Logical. If TRUE, drop non-counts layers (data /
#' scale.data) after loading to minimise the saved object size (default = FALSE).
#' @param verbose Logical. If TRUE, prints progress messages (default = FALSE).
#'
#' @return
#' If \code{rna_path} is .h5mu: a Seurat object with RNA assay, ADT assay,
#' and all embeddings from the RNA modality.
#'
#' If \code{rna_path}/\code{adt_path} are CSV/GZ/TAR/HDF5/AnnData:
#' a Seurat object with RNA assay (always) and ADT assay (if provided).
#'
#' If \code{rna_path} is a .rds file:
#' the Seurat object is returned **as is**. A bare matrix is rejected.
#'
#' If \code{output_dir} is specified, the result is also saved to disk and the
#' file path is returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Load a pre-processed Seurat object
#' cbmc_rna <- load_data(
#' rna_path = "./data/preprocessed_object.rds",
#' project_name = "PBMC",
#' species_name = "human",
#' verbose = TRUE
#' )
#'
#' # Load RNA (.h5) + ADT (.h5ad)
#' cbmc_full <- load_data(
#' rna_path = "./raw_files/COMBAT_rna.h5",
#' adt_path = "./raw_files/COMBAT_adt.h5ad",
#' project_name = "COMBAT_Project",
#' species_name = "mouse",
#' verbose = TRUE
#' )
#'
#' # Load from .h5mu file (automatically extracts RNA and ADT)
#' cbmc_mudata <- load_data(
#' rna_path = "./raw_files/multimodal_data.h5mu",
#' project_name = "MultiModal_Project",
#' verbose = TRUE
#' )
#'
#' # Load and save the result
#' cbmc <- load_data(
#' rna_path = "./data/preprocessed_object.rds",
#' project_name = "Merged_Project",
#' output_dir = "./results",
#' output_name = "loaded.rds",
#' verbose = TRUE
#' )
#'
#' # Load from config file
#' config <- yaml::read_yaml("load_data_config.yml")
#' cbmc_config <- do.call(load_data, config)
#' }
#' @export
load_data <- function(
    rna_path,
    adt_path = NULL,
    top_mouse_genes = 100,
    project_name = "SeuratProject",
    species_name = "human",
    rds_rna_key = "rna",
    rds_adt_key = "adt",
    output_dir = NULL,
    output_name = "loaded_data.rds",
    save_output = !is.null(output_dir),
    use_disk_backed = FALSE,
    bpcells_dir = NULL,
    slim_to_counts = FALSE,
    verbose = FALSE
) {

  #' Move RNA counts to an on-disk BPCells matrix (opt-in)
  #'
  #' Only used when `use_disk_backed = TRUE`; RNA only, fully guarded with a fallback to in-memory counts.
  #'
  #' @param obj Seurat object.
  #' @param assays Assays whose counts to move on disk.
  #' @return The object, with counts disk-backed when possible.
  .disk_back_counts <- function(obj, assays = "RNA") {
    if (!inherits(obj, "Seurat")) return(obj)
    do_disk <- isTRUE(use_disk_backed) && requireNamespace("BPCells", quietly = TRUE)
    if (isTRUE(use_disk_backed) && !requireNamespace("BPCells", quietly = TRUE) && verbose)
      message(" BPCells not installed - keeping counts in memory.")
    cache_dir <- if (!is.null(bpcells_dir) && nzchar(bpcells_dir)) bpcells_dir
                 else file.path(if (is.null(output_dir)) tempdir() else output_dir,
                                "bpcells_cache")
    if (do_disk) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    for (assay_name in intersect(assays, names(obj@assays))) {
      obj <- tryCatch({

        counts <- tryCatch(
          SeuratObject::LayerData(obj, assay = assay_name, layer = "counts"),
          error = function(e) tryCatch(
            SeuratObject::GetAssayData(obj, assay = assay_name, slot = "counts"),
            error = function(e2) NULL))
        if (is.null(counts)) {
          if (verbose) message(sprintf(" %s: no counts layer found - skipping.", assay_name))
          obj
        } else if (inherits(counts, "IterableMatrix")) {
          if (verbose) message(sprintf(" %s counts already on-disk - skipping.", assay_name))
          obj
        } else {
          mat_dir <- file.path(cache_dir, sprintf("%s_counts", assay_name))
          if (dir.exists(mat_dir)) unlink(mat_dir, recursive = TRUE)

          BPCells::write_matrix_dir(mat = counts, dir = mat_dir, overwrite = TRUE)
          disk <- BPCells::open_matrix_dir(dir = mat_dir)
          suppressWarnings(obj[[assay_name]] <- SeuratObject::CreateAssay5Object(counts = disk))
          if (verbose) message(sprintf(" %s counts moved on-disk (BPCells): %s",
                                       assay_name, mat_dir))
          gc()
          obj
        }
      }, error = function(e) {
        if (verbose) message(sprintf(
          " Disk-backing %s failed (%s) - kept in memory.",
          assay_name, conditionMessage(e)))
        obj
      })
    }
    }

    if (isTRUE(slim_to_counts)) {
      obj <- tryCatch({
        cl <- unique(unlist(lapply(
          names(obj@assays),
          function(a) tryCatch(
            SeuratObject::Layers(obj, assay = a, search = "counts"),
            error = function(e) character(0)))))
        if (length(cl) == 0) cl <- "counts"
        slim <- Seurat::DietSeurat(obj, layers = cl)
        if (verbose) message(sprintf(
          " Slimmed to counts layers (%s); data layer will be rebuilt by normalization.",
          paste(cl, collapse = ", ")))
        slim
      }, error = function(e) {
        if (verbose) message(sprintf(" Slim-to-counts skipped (%s)", conditionMessage(e)))
        obj
      })
    }
    gc()
    obj
  }

  #' Lower-case file extension of a path.
  get_ext <- function(path) {

    if (dir.exists(path)) return("")
    tolower(tools::file_ext(path))
  }

  #' Convert an .h5mu file to a Seurat object
  #'
  #' Reads the MuData through the Python muon backend, exports each modality
  #' (RNA and PROTEIN/ADT) to a temporary `.h5ad` carrying `X` plus the
  #' modality's integer counts layer when it has one, opens the matrices with
  #' BPCells, and assembles the Seurat object: each assay's `counts` comes
  #' from the integer layer (or `X` when none exists) and its `data` from a
  #' non-integer `X`. Global obs metadata and RNA obsm embeddings are carried
  #' over. The temporary files stay alive for the rest of the process so the
  #' returned object's disk-backed matrices remain readable until saved.
  #'
  #' @param h5mu_path Path to the `.h5mu` file.
  #' @param project_name Project name for the Seurat object.
  #' @param verbose Print progress.
  #' @return A Seurat object with RNA and ADT assays.
  load_h5mu_to_seurat <- function(h5mu_path, project_name, verbose) {
    if (verbose) message(sprintf("[%s] Loading MuData from: %s", Sys.time(), h5mu_path))

    .ensure_python(verbose)
    py_copy <- reticulate::import("copy", convert = FALSE)
    muon <- reticulate::import("muon", convert = FALSE)

    h5mu_py <- reticulate::py_run_string("
import pandas as _pd
import anndata as _ad

def _sf_plain_index(idx):
    return _pd.Index([str(x) for x in idx], name=idx.name)

def sf_sanitize_strings(df):
    df = df.copy()
    df.index = _sf_plain_index(df.index)
    for c in df.columns:
        if isinstance(df[c].dtype, _pd.StringDtype):
            df[c] = df[c].astype(object).where(df[c].notna(), None)
    return df

def sf_minimal_h5ad(mod, path):
    import numpy as _np
    a = _ad.AnnData(X=mod.X,
                    obs=_pd.DataFrame(index=_sf_plain_index(mod.obs.index)),
                    var=_pd.DataFrame(index=_sf_plain_index(mod.var.index)))
    for k in ('counts', 'raw'):
        if k in mod.layers:
            L = mod.layers[k]
            v = L.data[:2000] if hasattr(L, 'data') else _np.asarray(L)[:20].ravel()[:2000]
            v = _np.asarray(v)
            if v.size and bool(_np.all(v == _np.round(v))):
                a.layers['counts'] = L
                break
    a.write_h5ad(path)
", convert = FALSE)

    temp_dir <- tempfile("seurat_temp_")
    dir.create(temp_dir)

    on.exit(gc(full = TRUE))

    tryCatch({
      mdata <- muon$read_h5mu(h5mu_path)

      rna_mod <- mdata$mod$get("RNA")
      rna_mod$obs <- py_copy$deepcopy(mdata$obs)
      rna_path <- file.path(temp_dir, "rna_only.h5ad")
      h5mu_py$sf_minimal_h5ad(rna_mod, rna_path)

      adt_mod <- tryCatch({
        mdata$mod$get("PROTEIN")
      }, error = function(e) {
        tryCatch({
          mdata$mod$get("ADT")
        }, error = function(e2) {
          stop("Could not find ADT/PROTEIN modality in MuData file")
        })
      })

      adt_mod$obs <- py_copy$deepcopy(mdata$obs)
      adt_path <- file.path(temp_dir, "adt_only.h5ad")
      h5mu_py$sf_minimal_h5ad(adt_mod, adt_path)

      obs_path <- file.path(temp_dir, "rna_obs.csv")
      rna_obs_df <- reticulate::py_to_r(h5mu_py$sf_sanitize_strings(rna_mod$obs))
      write.csv(rna_obs_df, file = obs_path, row.names = TRUE)

      .guard_square_h5ad(rna_path)
      rna_x <- BPCells::open_matrix_anndata_hdf5(path = rna_path)
      rna_cnt <- tryCatch(
        BPCells::open_matrix_anndata_hdf5(path = rna_path, group = "layers/counts"),
        error = function(e) NULL)
      rna_dat <- NULL
      if (!is.null(rna_cnt)) {
        rna_dat <- rna_x
        if (verbose) message(" RNA layer mapping: counts <- layers/counts, data <- X")
      } else if (isFALSE(.mat_values_integer(rna_x))) {
        rna_cnt <- rna_x
        rna_dat <- rna_x
        if (verbose) message(" RNA layer mapping: counts <- X, data <- X (no counts layer; X is non-integer)")
      } else {
        rna_cnt <- rna_x
        if (verbose) message(" RNA layer mapping: counts <- X (integer values)")
      }
      .guard_square_h5ad(adt_path)
      adt_x <- BPCells::open_matrix_anndata_hdf5(path = adt_path)
      adt_cnt <- tryCatch(
        BPCells::open_matrix_anndata_hdf5(path = adt_path, group = "layers/counts"),
        error = function(e) NULL)

      seurat_obj <- CreateSeuratObject(counts = rna_cnt, assay = "RNA", project = project_name)
      if (!is.null(rna_dat)) {
        dimnames(rna_dat) <- list(rownames(seurat_obj), colnames(seurat_obj))
        seurat_obj[["RNA"]]$data <- rna_dat
      }
      if (!is.null(adt_cnt)) {
        seurat_obj[["ADT"]] <- CreateAssayObject(counts = methods::as(adt_cnt, "dgCMatrix"))
        adt_dm <- methods::as(adt_x, "dgCMatrix")
        dimnames(adt_dm) <- list(rownames(seurat_obj[["ADT"]]), colnames(seurat_obj))
        seurat_obj[["ADT"]] <- SeuratObject::SetAssayData(seurat_obj[["ADT"]],
                                                          layer = "data", new.data = adt_dm)
        if (verbose) message(" ADT layer mapping: counts <- layers/counts, data <- X")
      } else {
        seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_x)
      }

      metadata <- read.csv(obs_path, row.names = 1)
      seurat_obj@meta.data <- metadata

      obsm_keys <- .py_mapping_keys(rna_mod$obsm, verbose)

      for (key in obsm_keys) {
        embedding <- tryCatch({
          reticulate::py_to_r(rna_mod$obsm[[key]])
        }, error = function(e) {
          reticulate::py_to_r(rna_mod$obsm$get(key))
        })

        if (!is.matrix(embedding)) {
          embedding <- as.matrix(embedding)
        }

        if (nrow(embedding) != ncol(seurat_obj)) {
          if (verbose) message(sprintf("Skipping embedding %s: dimension mismatch", key))
          next
        }

        rownames(embedding) <- colnames(seurat_obj)
        if (is.null(colnames(embedding))) {
          colnames(embedding) <- paste0(key, "_", seq_len(ncol(embedding)))
        }

        if (ncol(embedding) >= 2) {
          reduc_key <- paste0(toupper(key), "_")
          seurat_obj[[key]] <- CreateDimReducObject(
            embeddings = embedding,
            key = reduc_key,
            assay = "RNA"
          )
          if (verbose) message(sprintf("Added embedding: %s with %d dimensions", key, ncol(embedding)))
        }
      }

      gc(full = TRUE)

      if (verbose) message("Successfully loaded MuData and converted to Seurat object.")
      return(seurat_obj)

    }, error = function(e) {
      unlink(temp_dir, recursive = TRUE)
      gc(full = TRUE)
      stop("Error creating Seurat object from MuData: ", e$message)
    })
  }

  #' Dimensions of an .h5ad matrix group, read without BPCells
  #'
  #' @param path Path to the `.h5ad` file.
  #' @param group HDF5 group of the matrix (default `"X"`).
  #' @return Integer vector `c(rows, cols)`, or `NULL` when hdf5r is missing
  #' or the group cannot be read.
  .h5ad_x_dims <- function(path, group = "X") {
    if (!requireNamespace("hdf5r", quietly = TRUE)) return(NULL)
    tryCatch({
      h5 <- hdf5r::H5File$new(path, mode = "r")
      on.exit(h5$close_all(), add = TRUE)
      if (!(group %in% names(h5))) return(NULL)
      node <- h5[[group]]
      d <- if (inherits(node, "H5Group")) {
        hdf5r::h5attributes(node)[["shape"]]
      } else {
        node$dims
      }
      if (length(d) == 2L) as.integer(d) else NULL
    }, error = function(e) NULL)
  }

  #' Refuse square .h5ad matrices before BPCells can segfault
  #'
  #' BPCells 0.3.1 does not raise a catchable R error on a square anndata
  #' matrix: `dims_matrix_anndata_hdf5_cpp` dereferences a null pointer and
  #' the whole process dies with SIGSEGV, printing nothing a caller can act
  #' on. Measured in this image: 5x5 through 100x100 square matrices all
  #' crash; every non-square shape reads correctly. Real CITE-seq inputs are
  #' never square, but a subset or toy file can be, so fail loudly instead of
  #' dying silently.
  #'
  #' @param path Path to the `.h5ad` file.
  #' @param group HDF5 group of the matrix to check.
  #' @return Invisibly `NULL`; stops with an explanatory error on a square
  #' matrix.
  .guard_square_h5ad <- function(path, group = "X") {
    d <- .h5ad_x_dims(path, group)
    if (is.null(d) || anyNA(d)) return(invisible(NULL))
    if (d[1] == d[2]) {
      stop(sprintf(paste0(
        "h5ad '%s' has a square X matrix (%d x %d). BPCells %s segfaults on ",
        "square anndata matrices, which would terminate this process with no ",
        "recoverable error. Subset the input so cells != genes, or convert it ",
        "to 10x h5/mtx first."),
        basename(path), d[1], d[2],
        as.character(utils::packageVersion("BPCells"))), call. = FALSE)
    }
    invisible(NULL)
  }

  #' Are a matrix's stored values all whole numbers?
  #'
  #' Samples stored values column-block by column-block (never materializing the matrix) to classify it as raw counts (integer) versus processed expression (non-integer).
  #'
  #' @param m Matrix of any supported class (BPCells, dgCMatrix, base).
  #' @param need Number of stored values to sample.
  #' @return `TRUE`, `FALSE`, or `NA` when the matrix holds no values.
  .mat_values_integer <- function(m, need = 500L) {
    vals <- numeric(0)
    j <- 1L
    while (length(vals) < need && j <= ncol(m)) {
      k <- min(j + 199L, ncol(m))
      vals <- c(vals, methods::as(m[, j:k, drop = FALSE], "dgCMatrix")@x)
      j <- k + 1L
    }
    if (!length(vals)) return(NA)
    all(vals == round(vals))
  }

  #' Find the raw-counts group of an .h5ad
  #'
  #' Checks `layers/counts`, then `layers/raw` (COMBAT pattern), then `raw/X` (Stephenson / CELLxGENE pattern), falling back to `"X"`.
  #'
  #' @param path Path to the `.h5ad` file.
  #' @param verbose Print which group was chosen.
  #' @return The HDF5 group to read counts from.
  .detect_h5ad_raw_group <- function(path, verbose = FALSE) {

    group <- "X"

    if (!requireNamespace("hdf5r", quietly = TRUE)) {
      if (verbose) message(" hdf5r not installed - reading from X (may be normalized)")
      return(group)
    }

    tryCatch({
      h5 <- hdf5r::H5File$new(path, mode = "r")
      on.exit(h5$close_all(), add = TRUE)

      if (h5$exists("layers") && h5[["layers"]]$exists("counts")) {
        cnt_grp <- h5[["layers"]][["counts"]]
        if (cnt_grp$exists("data") && cnt_grp$exists("indptr")) {
          group <- "layers/counts"
          if (verbose) message(" Found raw counts in layers/counts - using instead of X")
          return(group)
        }
      }

      if (h5$exists("layers") && h5[["layers"]]$exists("raw")) {
        raw_grp <- h5[["layers"]][["raw"]]
        if (raw_grp$exists("data") && raw_grp$exists("indptr")) {
          group <- "layers/raw"
          if (verbose) message(" Found raw counts in layers/raw - using instead of X")
          return(group)
        }
      }

      if (h5$exists("raw") && h5[["raw"]]$exists("X")) {
        raw_x <- h5[["raw"]][["X"]]
        has_sparse <- raw_x$exists("data") && raw_x$exists("indptr")
        has_dense <- inherits(raw_x, "H5D")
        if (has_sparse || has_dense) {
          group <- "raw/X"
          if (verbose) message(" Found raw counts in raw/X - using instead of X")
          return(group)
        }
      }

      if (verbose) message(" No raw count layer found - reading from X")

    }, error = function(e) {
      if (verbose) message(sprintf(" Could not probe for raw layers: %s", e$message))
    })

    return(group)
  }

  #' Extract ADT/protein data stored outside the main X matrix
  #'
  #' Supported storage patterns: `uns['antibody_X']` + `uns['antibody_features']` (Stephenson et al.), `uns['protein_X']` + `uns['protein_features']`, and dense `obsm['protein']` / `obsm['ADT']` / `obsm['antibody']` matrices.
  #'
  #' @param path Path to the `.h5ad` file.
  #' @param verbose Print which pattern matched.
  #' @return Proteins-by-cells matrix, or `NULL` when the file carries no
  #' separate ADT data.
  .extract_h5ad_adt <- function(path, verbose = FALSE) {

    if (requireNamespace("hdf5r", quietly = TRUE)) {
      tryCatch({
        h5 <- hdf5r::H5File$new(path, mode = "r")
        on.exit(h5$close_all(), add = TRUE)

        if (h5$exists("uns")) {
          uns <- h5[["uns"]]
          for (prefix in c("antibody", "protein", "ADT", "adt")) {
            mat_key <- paste0(prefix, "_X")
            feat_key <- paste0(prefix, "_features")

            if (uns$exists(mat_key) && uns$exists(feat_key)) {

              raw_mat_key <- paste0(prefix, "_raw.X")
              if (uns$exists(raw_mat_key)) {
                if (verbose) message(sprintf(" Found raw ADT data in uns/%s (preferred over uns/%s)",
                                             raw_mat_key, mat_key))
                mat_grp <- uns[[raw_mat_key]]
              } else {
                if (verbose) message(sprintf(" Found ADT data in uns/%s", mat_key))
                mat_grp <- uns[[mat_key]]
              }
              data_vec <- mat_grp[["data"]]$read()
              indices <- mat_grp[["indices"]]$read()
              indptr <- mat_grp[["indptr"]]$read()

              grp_attrs <- tryCatch(hdf5r::h5attr_names(mat_grp), error = function(e) character(0))
              shape <- if ("shape" %in% grp_attrs) as.integer(hdf5r::h5attr(mat_grp, "shape")) else NULL
              encoding <- if ("encoding-type" %in% grp_attrs) tolower(hdf5r::h5attr(mat_grp, "encoding-type")) else NULL

              feat_grp <- uns[[feat_key]]
              feat_names <- if (feat_grp$exists("_index")) {
                feat_grp[["_index"]]$read()
              } else if (feat_grp$exists(feat_key)) {
                feat_grp[[feat_key]]$read()
              } else NULL
              n_feat <- if (!is.null(feat_names)) length(feat_names) else NA_integer_

              ip_len <- length(indptr) - 1L
              if (!is.null(shape) && length(shape) == 2L) {
                R <- shape[1]; C <- shape[2]
              } else {
                other <- max(indices) + 1L
                if (is.null(encoding) || grepl("csc", encoding)) { R <- other; C <- ip_len }
                else { R <- ip_len; C <- other }
              }

              if (is.null(encoding))
                encoding <- if (ip_len == C) "csc_matrix" else if (ip_len == R) "csr_matrix" else "csr_matrix"

              if (grepl("csc", encoding)) {
                stored <- Matrix::sparseMatrix(
                  i = indices + 1L, p = indptr, x = as.numeric(data_vec), dims = c(R, C)
                )
              } else {
                stored <- Matrix::sparseMatrix(
                  i = rep(seq_len(R), diff(indptr)), j = indices + 1L,
                  x = as.numeric(data_vec), dims = c(R, C)
                )
              }

              adt_mat <- if (!is.na(n_feat) && nrow(stored) == n_feat) stored
                         else if (!is.na(n_feat) && ncol(stored) == n_feat) Matrix::t(stored)
                         else Matrix::t(stored)
              if (is.null(feat_names) || length(feat_names) != nrow(adt_mat))
                feat_names <- paste0("ADT_", seq_len(nrow(adt_mat)))
              rownames(adt_mat) <- feat_names

              if (verbose)
                message(sprintf(" ADT matrix: %d proteins x %d cells (%s)",
                                nrow(adt_mat), ncol(adt_mat), encoding))
              return(adt_mat)
            }
          }
        }
      }, error = function(e) {
        if (verbose) message(sprintf(" hdf5r scan for ADT failed: %s", e$message))
      })
    }

    tryCatch({
      ad <- reticulate::import("anndata", convert = FALSE)
      adata <- ad$read_h5ad(path)
      obsm_keys <- .py_mapping_keys(adata$obsm, verbose)
      adt_keys <- c("protein", "ADT", "adt", "antibody", "X_protein",
                    "X_ADT", "X_adt", "X_antibody", "cite", "CITE")
      for (key in intersect(adt_keys, obsm_keys)) {
        adt_obsm <- as.matrix(reticulate::py_to_r(adata$obsm[[key]]))
        if (ncol(adt_obsm) > 0) {
          adt_mat <- Matrix::Matrix(t(adt_obsm), sparse = TRUE)
          if (is.null(rownames(adt_mat)))
            rownames(adt_mat) <- paste0("ADT_", seq_len(nrow(adt_mat)))
          if (verbose)
            message(sprintf(" Found ADT in obsm['%s']: %d proteins x %d cells",
                            key, nrow(adt_mat), ncol(adt_mat)))
          return(adt_mat)
        }
      }
    }, error = function(e) {
      if (verbose) message(sprintf(" obsm scan for ADT failed: %s", e$message))
    })

    if (verbose) message(" No ADT/protein data found in h5ad file")
    return(NULL)
  }

  if (!is.null(adt_path)) {
    stop(paste0(
      "adt_path is retired. Every format this module accepts (.rds, .h5ad, ",
      ".h5mu) carries RNA and ADT in one object and splits ",
      "them internally, so a second matrix path has nothing to point at. ",
      "Pass the single pre-processed object as rna_path."), call. = FALSE)
  }

  #' Load one input file into matrices ready for object assembly
  #'
  #' Dispatches on the file extension: `.h5mu` returns a finished Seurat object, `.h5ad` opens counts (and, when the file declares one, the processed matrix for the `data` layer) with BPCells and splits interleaved ADT features, `.rds` returns the pre-processed Seurat object as-is.
  #'
  #' @param path Input file path.
  #' @param type `"rna"` or `"adt"` (historical; ADT now arrives inside the
  #' one input file).
  #' @param rds_key Unused for single-object inputs.
  #' @return List with `mat`, `ext` and, for `.h5ad`, `obs`, `obsm`,
  #' `adt_mat`, `data_mat`, `adt_data`.
  load_matrix <- function(path, type = c("rna", "adt"), rds_key = NULL) {
    type <- match.arg(type)
    ext <- get_ext(path)
    if (verbose) message(sprintf("[%s] Loading %s data from: %s", Sys.time(), type, path))

    if (type == "rna" && ext == "h5mu") {
      return(list(mat = load_h5mu_to_seurat(path, project_name, verbose), ext = "h5mu_seurat"))
    }

    if (ext == "h5ad") {

      h5ad_group <- .detect_h5ad_raw_group(path, verbose = verbose)
      .guard_square_h5ad(path, h5ad_group)
      mat <- BPCells::open_matrix_anndata_hdf5(path = path, group = h5ad_group)
      h5ad_data <- NULL
      h5ad_adt_data <- NULL
      if (!identical(h5ad_group, "X")) {
        h5ad_data <- tryCatch(
          BPCells::open_matrix_anndata_hdf5(path = path, group = "X"),
          error = function(e) NULL)
        if (!is.null(h5ad_data) &&
            (!identical(dim(h5ad_data), dim(mat)) ||
             isTRUE(.mat_values_integer(h5ad_data)))) h5ad_data <- NULL
        if (verbose) message(sprintf(
          " Layer mapping: counts <- %s%s", h5ad_group,
          if (is.null(h5ad_data)) "" else ", data <- X"))
      } else if (isFALSE(.mat_values_integer(mat))) {
        h5ad_data <- mat
        if (verbose) message(" Layer mapping: counts <- X, data <- X (no raw layer; X is non-integer)")
      } else if (verbose) {
        message(" Layer mapping: counts <- X (integer values)")
      }
      mat <- .convert_gene_ids(mat, species = species_name, verbose = verbose)
      if (!is.null(h5ad_data)) {
        h5ad_data <- if (identical(h5ad_group, "X")) mat
                     else .convert_gene_ids(h5ad_data, species = species_name, verbose = FALSE)
      }

      h5ad_meta <- .read_h5ad_metadata(path, verbose = verbose)

      all_features <- rownames(mat)
      adt_mask <- grepl("^AB_", all_features, ignore.case = FALSE)

      if (sum(adt_mask) > 0 && sum(!adt_mask) > 0) {
        if (verbose) message(sprintf(
          " Detected %d interleaved ADT features (AB_ prefix) among %d total features - splitting.",
          sum(adt_mask), length(all_features)))

        h5ad_adt <- mat[adt_mask, ]
        mat <- mat[!adt_mask, ]
        if (!is.null(h5ad_data)) {
          h5ad_adt_data <- h5ad_data[adt_mask, ]
          h5ad_data <- h5ad_data[!adt_mask, ]
        }

        if (verbose) {
          message(sprintf(" RNA matrix: %d genes x %d cells", nrow(mat), ncol(mat)))
          message(sprintf(" ADT matrix: %d proteins x %d cells", nrow(h5ad_adt), ncol(h5ad_adt)))
        }
      } else {

        h5ad_adt <- .extract_h5ad_adt(path, verbose = verbose)
      }

    } else if (ext == "rds") {
      mat <- readRDS(path)

      if (!inherits(mat, "Seurat")) {
        stop(sprintf(paste0(
          "%s contains a '%s', not a Seurat object. This module loads a ",
          "pre-processed Seurat object; a bare matrix has not been through QC, ",
          "normalisation or clustering."), basename(path), class(mat)[1]),
          call. = FALSE)
      }
      if (verbose) message("Detected RDS format. Returning object as is.")

    } else {
      stop(sprintf(paste0(
        "Unsupported input format '.%s' for %s data. This module loads a ",
        "PRE-PROCESSED object; supported formats are .rds, .h5ad and .h5mu - ",
        "one file, not a directory. Raw formats (.csv, .csv.gz, .tar.gz, 10x ",
        "CellRanger .h5) were removed deliberately: loading raw counts here ",
        "yields an object that looks valid but has had no QC, normalisation ",
        "or clustering."), ext, type), call. = FALSE)
    }

    obs <- if (exists("h5ad_meta") && !is.null(h5ad_meta)) h5ad_meta$obs else NULL
    obsm <- if (exists("h5ad_meta") && !is.null(h5ad_meta)) h5ad_meta$obsm else NULL
    adt_mat <- if (exists("h5ad_adt")) h5ad_adt else NULL
    data_mat <- if (exists("h5ad_data")) h5ad_data else NULL
    adt_data <- if (exists("h5ad_adt_data")) h5ad_adt_data else NULL

    return(list(mat = mat, ext = ext, obs = obs, obsm = obsm, adt_mat = adt_mat,
                data_mat = data_mat, adt_data = adt_data))
  }

  rna_load <- load_matrix(rna_path, type = "rna", rds_key = rds_rna_key)
  cbmc.rna <- rna_load$mat
  rna_ext <- rna_load$ext

  if (rna_ext == "h5mu_seurat") {
    if (verbose) message("Returning Seurat object from MuData file.")

    if (save_output && !is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        if (verbose) message(sprintf("Created output directory: %s", output_dir))
      }
      output_path <- file.path(output_dir, output_name)
      if (verbose) message(sprintf("[%s] Saving output to: %s", Sys.time(), output_path))
      saveRDS(cbmc.rna, file = output_path, compress = TRUE)
      if (verbose) {
        file_size_mb <- file.size(output_path) / 1024^2
        message(sprintf("Saved successfully (%.2f MB)", file_size_mb))
      }
      attr(cbmc.rna, "output_path") <- output_path
    }

    return(cbmc.rna)
  }

  if (rna_ext == "rds") {

    result <- cbmc.rna

    result <- .disk_back_counts(result)

    if (save_output && !is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        if (verbose) message(sprintf("Created output directory: %s", output_dir))
      }
      output_path <- file.path(output_dir, output_name)
      if (verbose) message(sprintf("[%s] Saving output to: %s", Sys.time(), output_path))
      saveRDS(result, file = output_path, compress = TRUE)
      if (verbose) {
        file_size_mb <- file.size(output_path) / 1024^2
        message(sprintf("Saved successfully (%.2f MB)", file_size_mb))
      }
      attr(result, "output_path") <- output_path
    }

    return(result)
  }

  cbmc <- CreateSeuratObject(counts = cbmc.rna, project = project_name)

  if (!is.null(rna_load$data_mat)) {
    dm <- rna_load$data_mat
    dimnames(dm) <- list(rownames(cbmc), colnames(cbmc))
    cbmc[["RNA"]]$data <- dm
    if (verbose) message("Attached the file's processed expression as the 'data' layer")
  }

  if (!is.null(rna_load$obs)) {

    obs_df <- rna_load$obs

    if (!is.null(rownames(obs_df)) && length(intersect(rownames(obs_df), colnames(cbmc))) > 0) {
      common <- intersect(colnames(cbmc), rownames(obs_df))
      obs_df <- obs_df[common, , drop = FALSE]
      cbmc <- cbmc[, common]
    } else if (nrow(obs_df) == ncol(cbmc)) {
      rownames(obs_df) <- colnames(cbmc)
    } else {
      if (verbose) warning(sprintf(
        "h5ad obs has %d rows but Seurat has %d cells - skipping obs attachment",
        nrow(obs_df), ncol(cbmc)))
      obs_df <- NULL
    }
    if (!is.null(obs_df)) {

      new_cols <- setdiff(colnames(obs_df), colnames(cbmc@meta.data))
      if (length(new_cols) > 0) {
        cbmc@meta.data[, new_cols] <- obs_df[colnames(cbmc), new_cols, drop = FALSE]
        if (verbose) message(sprintf("Added %d metadata columns from h5ad obs: %s",
                                     length(new_cols), paste(head(new_cols, 10), collapse = ", ")))
      }
    }
  }

  if (!is.null(rna_load$obsm)) {
    for (key in names(rna_load$obsm)) {
      emb <- rna_load$obsm[[key]]
      if (!is.matrix(emb)) emb <- as.matrix(emb)
      if (nrow(emb) != ncol(cbmc)) {
        if (verbose) message(sprintf("Skipping obsm '%s': %d rows vs %d cells", key, nrow(emb), ncol(cbmc)))
        next
      }
      rownames(emb) <- colnames(cbmc)
      if (is.null(colnames(emb)))
        colnames(emb) <- paste0(key, "_", seq_len(ncol(emb)))
      if (ncol(emb) >= 2) {
        reduc_key <- gsub("[^A-Za-z0-9]", "", toupper(key))
        reduc_key <- paste0(reduc_key, "_")
        cbmc[[key]] <- CreateDimReducObject(
          embeddings = emb, key = reduc_key, assay = "RNA"
        )
        if (verbose) message(sprintf("Added embedding '%s' (%d dims) from h5ad obsm", key, ncol(emb)))
      }
    }
  }

  if (is.null(adt_path) && !is.null(rna_load$adt_mat)) {
    adt_from_h5ad <- rna_load$adt_mat
    if (!is.matrix(adt_from_h5ad) && !methods::is(adt_from_h5ad, "Matrix"))
      adt_from_h5ad <- methods::as(adt_from_h5ad, "dgCMatrix")

    if (is.null(colnames(adt_from_h5ad)) && ncol(adt_from_h5ad) == ncol(cbmc)) {
      colnames(adt_from_h5ad) <- colnames(cbmc)
    } else if (!is.null(colnames(adt_from_h5ad))) {
      colnames(adt_from_h5ad) <- gsub("_", "-", colnames(adt_from_h5ad))
      common_cells <- intersect(colnames(cbmc), colnames(adt_from_h5ad))
      if (length(common_cells) > 0) {
        cbmc <- cbmc[, common_cells]
        adt_from_h5ad <- adt_from_h5ad[, common_cells]
      } else if (ncol(adt_from_h5ad) == ncol(cbmc)) {

        colnames(adt_from_h5ad) <- colnames(cbmc)
      }
    }

    if (ncol(adt_from_h5ad) == ncol(cbmc)) {
      cbmc[["ADT"]] <- CreateAssayObject(counts = adt_from_h5ad)
      if (verbose) message(sprintf("Created ADT assay from h5ad: %d proteins x %d cells",
                                   nrow(adt_from_h5ad), ncol(adt_from_h5ad)))
      if (!is.null(rna_load$adt_data)) {
        src <- rna_load$adt_data
        orig <- colnames(src)
        pos <- if (is.null(orig)) seq_len(ncol(cbmc))
                else match(colnames(cbmc), gsub("_", "-", orig))
        if (!anyNA(pos) && length(pos) == ncol(cbmc)) {
          adt_dm <- methods::as(src[, pos, drop = FALSE], "dgCMatrix")
          dimnames(adt_dm) <- list(rownames(cbmc[["ADT"]]), colnames(cbmc))
          cbmc[["ADT"]] <- SeuratObject::SetAssayData(cbmc[["ADT"]], layer = "data",
                                                      new.data = adt_dm)
          if (verbose) message("Attached processed ADT values as the ADT 'data' slot")
        }
      }
    } else if (verbose) {
      warning(sprintf("ADT from h5ad has %d cells but Seurat has %d - skipping ADT assay creation",
                      ncol(adt_from_h5ad), ncol(cbmc)))
    }
  }

  if (!is.null(adt_path)) {
    adt_load <- load_matrix(adt_path, type = "adt", rds_key = rds_adt_key)
    cbmc.adt <- adt_load$mat

    colnames(cbmc.adt) <- gsub("_", "-", colnames(cbmc.adt))

    common_cells <- intersect(colnames(cbmc), colnames(cbmc.adt))
    if (length(common_cells) == 0) stop("No overlapping cells between RNA and ADT.")
    cbmc <- cbmc[, common_cells]
    cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt[, common_cells])
  }

  cbmc <- .disk_back_counts(cbmc)

  if (save_output && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) message(sprintf("Created output directory: %s", output_dir))
    }
    output_path <- file.path(output_dir, output_name)
    if (verbose) message(sprintf("[%s] Saving output to: %s", Sys.time(), output_path))
    saveRDS(cbmc, file = output_path, compress = TRUE)
    if (verbose) {
      file_size_mb <- file.size(output_path) / 1024^2
      message(sprintf("Saved successfully (%.2f MB)", file_size_mb))
    }
    attr(cbmc, "output_path") <- output_path
  }

  if (verbose) message("load_data() completed successfully.")
  return(cbmc)
}
