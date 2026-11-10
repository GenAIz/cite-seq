FIXTURES <- testthat::test_path("fixtures")

TINY_RNA_CSV <- file.path(FIXTURES, "tiny_rna.csv")
TINY_ADT_CSV <- file.path(FIXTURES, "tiny_adt.csv")
TINY_META_CSV <- file.path(FIXTURES, "tiny_metadata.csv")

make_tiny_seurat <- function(n_cells = 50, n_genes = 200, n_adt = 15,
                              seed = 42, add_adt = TRUE) {
  set.seed(seed)
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(Matrix))

  gene_names <- paste0("GENE", sprintf("%04d", seq_len(n_genes)))
  cell_names <- paste0("Cell", sprintf("%04d", seq_len(n_cells)))
  adt_names  <- c("CD3","CD4","CD8a","CD19","CD56","CD14","CD16",
                  "CD25","CD27","CD45RA","CD45RO","CD127","PD-1","CTLA-4","CD28")[
                  seq_len(n_adt)]

  rna_counts <- matrix(0L, nrow = n_genes, ncol = n_cells,
                       dimnames = list(gene_names, cell_names))
  for (i in seq_len(n_genes)) {
    base <- sample(0:4, n_cells, replace = TRUE)
    if (i <= 20) base[26:n_cells] <- base[26:n_cells] + sample(5:15, n_cells - 25, replace = TRUE)
    rna_counts[i, ] <- base
  }
  rna_sparse <- as(rna_counts, "sparseMatrix")

  obj <- CreateSeuratObject(counts = rna_sparse, project = "TestProject")

  obj$orig.ident <- ifelse(seq_len(n_cells) <= 25, "Sample1", "Sample2")
  obj$Batch      <- ifelse(seq_len(n_cells) <= 25, "BatchA", "BatchB")
  obj$Condition  <- rep(c("Control", "Treated"), length.out = n_cells)
  obj$Donor      <- paste0("D", sprintf("%02d", ceiling(seq_len(n_cells) / 10)))
  obj$CellType   <- rep(c("CD4 T", "CD8 T", "NK", "B cell", "Monocyte"),
                        length.out = n_cells)

  if (add_adt) {
    adt_counts <- matrix(
      sample(0:200, n_adt * n_cells, replace = TRUE),
      nrow = n_adt, ncol = n_cells,
      dimnames = list(adt_names, cell_names)
    )
    obj[["ADT"]] <- CreateAssayObject(counts = as(adt_counts, "sparseMatrix"))
  }

  obj
}

make_processed_seurat <- function(seed = 42, add_adt = TRUE) {
  obj <- make_tiny_seurat(seed = seed, add_adt = add_adt)
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 50, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, npcs = 10, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:5, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:5, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.5, cluster.name = "RNA_clusters",
                        verbose = FALSE)
    if (add_adt) {
      DefaultAssay(obj) <- "ADT"
      obj <- NormalizeData(obj, normalization.method = "CLR",
                           margin = 2, verbose = FALSE)
      obj <- ScaleData(obj, features = rownames(obj[["ADT"]]), verbose = FALSE)
      obj <- RunPCA(obj, assay = "ADT",
                   features = rownames(obj[["ADT"]]),
                   npcs = min(10, nrow(obj[["ADT"]]) - 1),
                   reduction.name = "apca", reduction.key = "apca_",
                   verbose = FALSE)
      DefaultAssay(obj) <- "RNA"
    }
  })
  obj
}

skip_if_no_py_anndata <- function() {
  testthat::skip_if_not(requireNamespace("reticulate", quietly = TRUE),
                        "reticulate not installed")
  if (exists(".ensure_python", mode = "function")) .ensure_python(verbose = FALSE)
  testthat::skip_if_not(reticulate::py_module_available("anndata"),
                        "python anndata not available")
}

.write_test_h5ad <- function(path, n_cells = 10L, n_genes = 5L,
                             obs = NULL, obsm = NULL, seed = 42L) {
  set.seed(seed)
  ad <- reticulate::import("anndata",      convert = FALSE)
  np <- reticulate::import("numpy",        convert = FALSE)
  pd <- reticulate::import("pandas",       convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  cells <- paste0("Cell", seq_len(n_cells))
  genes <- paste0("Gene", seq_len(n_genes))

  dense <- np$array(matrix(as.numeric(rpois(n_cells * n_genes, 5)),
                           nrow = n_cells, ncol = n_genes))
  X <- sp$csr_matrix(dense$astype("float32"))

  obs_py <- if (is.null(obs) || length(obs) == 0) {
    pd$DataFrame(index = reticulate::r_to_py(cells))
  } else {
    pd$DataFrame(reticulate::dict(obs), index = reticulate::r_to_py(cells))
  }
  var_py <- pd$DataFrame(index = reticulate::r_to_py(genes))

  a <- if (is.null(obsm) || length(obsm) == 0) {
    ad$AnnData(X = X, obs = obs_py, var = var_py)
  } else {
    a_obsm <- reticulate::dict(
      lapply(obsm, function(m) np$array(m)$astype("float32")))
    ad$AnnData(X = X, obs = obs_py, var = var_py, obsm = a_obsm)
  }
  a$write_h5ad(path)
  invisible(path)
}

.h5ad_index_storage <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) return("hdf5r not installed")

  tryCatch({
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    if (!("obs" %in% names(h5))) return("no obs group")
    obs <- h5[["obs"]]

    idx <- tryCatch(hdf5r::h5attributes(obs)[["_index"]], error = function(e) NULL)
    if (is.null(idx) || !nzchar(idx))  return("no _index attribute")
    if (!(idx %in% names(obs)))        return(sprintf("index '%s' absent", idx))
    node <- obs[[idx]]
    if (inherits(node, "H5D"))          "dataset"
    else if (inherits(node, "H5Group")) "group"
    else                                class(node)[1]
  }, error = function(e) paste("error:", conditionMessage(e)))
}

.write_test_h5mu <- function(path, n_cells = 12L, n_genes = 20L, n_adt = 5L,
                             obsm = NULL, seed = 42L) {
  stopifnot(n_cells != n_genes, n_cells != n_adt)
  set.seed(seed)
  ad <- reticulate::import("anndata",      convert = FALSE)
  md <- reticulate::import("mudata",       convert = FALSE)
  np <- reticulate::import("numpy",        convert = FALSE)
  pd <- reticulate::import("pandas",       convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  cells <- paste0("Cell", seq_len(n_cells))
  idx   <- reticulate::r_to_py(cells)

  mk <- function(nc, nf, lambda) sp$csr_matrix(
    np$array(matrix(as.numeric(rpois(nc * nf, lambda)), nrow = nc, ncol = nf))$astype("float32"))

  rna <- ad$AnnData(
    X   = mk(n_cells, n_genes, 5),
    obs = pd$DataFrame(index = idx),
    var = pd$DataFrame(index = reticulate::r_to_py(paste0("Gene", seq_len(n_genes)))))
  if (!is.null(obsm) && length(obsm) > 0) {
    rna$obsm <- reticulate::dict(
      lapply(obsm, function(m) np$array(m)$astype("float32")))
  }

  adt <- ad$AnnData(
    X   = mk(n_cells, n_adt, 20),
    obs = pd$DataFrame(index = idx),
    var = pd$DataFrame(index = reticulate::r_to_py(paste0("ADT", seq_len(n_adt)))))

  mdata <- md$MuData(reticulate::dict(list(RNA = rna, PROTEIN = adt)))
  mdata$obs <- pd$DataFrame(
    reticulate::dict(list(donor = rep(c("D1", "D2"), length.out = n_cells))),
    index = idx)
  mdata$write(path)
  invisible(path)
}

skip_if_no_py_mudata <- function() {
  skip_if_no_py_anndata()
  testthat::expect_true(reticulate::py_module_available("mudata"),
                        label = "python mudata importable")
  testthat::expect_true(reticulate::py_module_available("muon"),
                        label = "python muon importable (scanpy/numpy compatible)")
}
