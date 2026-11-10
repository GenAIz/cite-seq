`%||%` <- function(a, b) if (!is.null(a)) a else b

make_valid_config <- function(tmp_dir = tempdir()) {
  list(
    metadata = list(project_name = "TestProject"),
    paths = list(
      input_data = list(
        rna_path      = file.path(tmp_dir, "fake_rna.rds"),
        adt_path      = NULL,
        metadata_file = NULL
      ),
      output_dir = file.path(tmp_dir, "output"),
      logs_dir   = file.path(tmp_dir, "logs")
    ),
    execution = list(
      mode    = "step-by-step",
      steps   = c("load_data"),
      verbose = FALSE,
      n_cores = 1
    ),
    load_data = list(
      project_name = "TestProject",
      species_name = "human"
    ),
    add_metadata = list(enabled = FALSE),
    processing = list(
      qc_filters = list(
        min_features  = 200,
        max_features  = 6000,
        max_percent_mt = 20
      ),
      clustering = list(rna_resolution = 1.5)
    ),
    aggregation = list(
      participant_col = "SampleID",
      celltype_cols   = "CellType"
    )
  )
}

write_valid_config <- function(tmp_dir = NULL) {
  if (is.null(tmp_dir)) {
    tmp_dir <- tempfile("orch_test_")
    dir.create(tmp_dir, recursive = TRUE)
  }

  cfg <- make_valid_config(tmp_dir)

  writeLines("placeholder", cfg$paths$input_data$rna_path)

  dir.create(cfg$paths$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$paths$logs_dir,   recursive = TRUE, showWarnings = FALSE)

  config_path <- file.path(tmp_dir, "test_config.json")
  jsonlite::write_json(cfg, config_path, auto_unbox = TRUE, pretty = TRUE)

  config_path
}

make_cite_seurat <- function(n_cells = 60, n_genes = 80, seed = 1) {
  suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  cells <- paste0("Cell", seq_len(n_cells))

  rna <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes,
                dimnames = list(genes, cells))

  rna[(n_genes - 9):n_genes, ] <- 0L
  rna[(n_genes - 9):n_genes, 1:2] <- 1L
  obj <- CreateSeuratObject(counts = as(rna, "sparseMatrix"))

  adt_names <- c("Ab.CD3", "Ab.CD4", "Ab.CD8", "Ab.CD19", "Ab.CD56",
                 "Ab.Isotype-MOPC.21", "Set1-Sample1", "Set1-Sample2")
  adt <- matrix(rpois(length(adt_names) * n_cells, lambda = 6),
                nrow = length(adt_names), dimnames = list(adt_names, cells))
  obj[["ADT"]] <- CreateAssayObject(counts = as(adt, "sparseMatrix"))

  obj$orig.ident <- c(rep("BatchA", n_cells - 3), rep("BatchTiny", 3))
  obj$Participant <- rep(c("P1", "P2"), length.out = n_cells)
  obj$Timepoint   <- rep(c("T0", "T1"), length.out = n_cells)
  obj$CellType    <- rep(c("CD4", "CD8", "NK"), length.out = n_cells)
  obj
}

.write_test_h5ad <- function(path, n_cells, n_genes, obs_col = NULL,
                             obs_values = NULL, obsm = NULL) {
  stopifnot(n_cells != n_genes)
  ad <- reticulate::import("anndata",       convert = FALSE)
  np <- reticulate::import("numpy",         convert = FALSE)
  pd <- reticulate::import("pandas",        convert = FALSE)
  sp <- reticulate::import("scipy.sparse",  convert = FALSE)
  cells <- paste0("Cell", seq_len(n_cells))
  genes <- paste0("Gene", seq_len(n_genes))
  dense <- np$array(matrix(as.numeric(rpois(n_cells * n_genes, 5)), nrow = n_cells))
  X   <- sp$csr_matrix(dense$astype("float32"))
  obs <- if (is.null(obs_col)) {
    pd$DataFrame(index = reticulate::r_to_py(cells))
  } else {
    pd$DataFrame(reticulate::dict(setNames(list(obs_values), obs_col)),
                 index = reticulate::r_to_py(cells))
  }
  var <- pd$DataFrame(index = reticulate::r_to_py(genes))
  a <- if (is.null(obsm)) {
    ad$AnnData(X = X, obs = obs, var = var)
  } else {
    ad$AnnData(X = X, obs = obs, var = var,
               obsm = reticulate::dict(
                 lapply(obsm, function(m) np$array(m)$astype("float32"))))
  }
  a$write_h5ad(path)
  invisible(path)
}
