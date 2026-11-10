skip_ld <- function() {
  skip_if_not_installed("Seurat")
  if (!exists(".convert_gene_ids", mode = "function"))
    skip("load_data_updated.R not sourced (Seurat unavailable at runner load)")
}

test_that("gene symbols pass through unchanged", {
  skip_ld()
  mat <- matrix(1:9, nrow = 3,
                dimnames = list(c("CD3E", "CD4", "CD8A"), paste0("C", 1:3)))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_identical(rownames(result), c("CD3E", "CD4", "CD8A"))
})

test_that("non-Ensembl IDs pass through unchanged", {
  skip_ld()
  genes <- c("ENS_FAKE", "NOTENSEMBL", "GENE001")
  mat   <- matrix(0, nrow = 3, ncol = 2, dimnames = list(genes, c("A", "B")))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_identical(rownames(result), genes)
})

test_that("matrix with NULL rownames passes through unchanged", {
  skip_ld()
  mat    <- matrix(1:4, nrow = 2)
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_null(rownames(result))
})

test_that("unsupported species throws an informative error", {
  skip_ld()
  genes <- paste0("ENSG", sprintf("%011d", 1:5))
  mat   <- matrix(0, nrow = 5, ncol = 2, dimnames = list(genes, c("A", "B")))
  expect_error(.convert_gene_ids(mat, species = "zebrafish"), regexp = "species")
})

test_that("AnnotationDbi conversion is attempted for Ensembl IDs", {
  skip_ld()
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  genes  <- c("ENSG00000139618", "ENSG00000141736")
  mat    <- matrix(0, nrow = 2, ncol = 2, dimnames = list(genes, c("A", "B")))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_false(any(grepl("^ENSG", rownames(result))))
  expect_equal(length(rownames(result)), 2L)
})

test_that(".read_h5ad_metadata returns list(obs=, obsm=) even on failure", {
  skip_ld()
  result <- .read_h5ad_metadata("/nonexistent/file.h5ad", verbose = FALSE)
  expect_type(result, "list")
  expect_true(all(c("obs", "obsm") %in% names(result)))
  expect_null(result$obs)
  expect_null(result$obsm)
})

test_that(".read_h5ad_metadata extracts obs via the anndata backend", {
  skip_ld()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE) &&
              reticulate::py_module_available("anndata"),
              "python anndata backend not available")
  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  tryCatch(
    .write_test_h5ad(tmp, n_cells = 10, n_genes = 5,
                     obs_col = "batch", obs_values = rep(c("A", "B"), each = 5)),
    error = function(e) skip(paste("Cannot create test h5ad:", e$message)))
  result <- .read_h5ad_metadata(tmp, verbose = FALSE)
  expect_s3_class(result$obs, "data.frame")
  expect_equal(nrow(result$obs), 10L)
  expect_true("batch" %in% colnames(result$obs))
})

test_that("load_data from h5ad attaches obs metadata to the Seurat object", {
  skip_ld()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE) &&
              reticulate::py_module_available("anndata"),
              "python anndata backend not available")
  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  tryCatch(
    .write_test_h5ad(tmp, n_cells = 20, n_genes = 30,
                     obs_col = "sample_id", obs_values = rep(c("S1", "S2"), each = 10)),
    error = function(e) skip(paste("Cannot create test h5ad:", e$message)))
  obj <- load_data(rna_path = tmp, project_name = "H5AD_Test", verbose = FALSE)
  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 20L)
  expect_true("sample_id" %in% colnames(obj@meta.data))

  expect_identical(as.character(unique(obj$orig.ident)), "H5AD_Test")
})

test_that(".read_h5ad_metadata extracts obsm embeddings", {
  skip_ld()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE) &&
              reticulate::py_module_available("anndata"),
              "python anndata backend not available")
  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  tryCatch(
    .write_test_h5ad(tmp, n_cells = 10, n_genes = 5,
                     obsm = list(X_umap = matrix(rnorm(20), nrow = 10, ncol = 2))),
    error = function(e) skip(paste("Cannot create test h5ad:", e$message)))

  result <- .read_h5ad_metadata(tmp, verbose = FALSE)
  expect_true("X_umap" %in% names(result$obsm))
  expect_equal(nrow(result$obsm$X_umap), 10L)
  expect_equal(ncol(result$obsm$X_umap), 2L)

  expect_false(any(c("attrname", "axes", "dim", "is_view", "as_dict") %in%
                   names(result$obsm)))
})

test_that("load_data from h5ad attaches obsm as reductions", {
  skip_ld()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE) &&
              reticulate::py_module_available("anndata"),
              "python anndata backend not available")
  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  tryCatch(
    .write_test_h5ad(tmp, n_cells = 20, n_genes = 30,
                     obsm = list(X_umap = matrix(rnorm(40), nrow = 20, ncol = 2),
                                 X_pca  = matrix(rnorm(200), nrow = 20, ncol = 10))),
    error = function(e) skip(paste("Cannot create test h5ad:", e$message)))

  obj <- load_data(rna_path = tmp, project_name = "H5AD_Obsm", verbose = FALSE)
  expect_true("X_umap" %in% names(obj@reductions))
  expect_true("X_pca"  %in% names(obj@reductions))
  expect_equal(ncol(Embeddings(obj, "X_umap")), 2L)
  expect_equal(ncol(Embeddings(obj, "X_pca")), 10L)
})
