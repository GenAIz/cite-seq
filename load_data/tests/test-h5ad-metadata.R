test_that(".read_h5ad_metadata returns list(obs=, obsm=) even when it fails", {
  result <- .read_h5ad_metadata("/nonexistent/file.h5ad", verbose = FALSE)
  expect_type(result, "list")
  expect_true("obs" %in% names(result))
  expect_true("obsm" %in% names(result))
  expect_null(result$obs)
  expect_null(result$obsm)
})

test_that("obs/_index is written as a dataset, not a group", {
  skip_if_no_py_anndata()
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 6, n_genes = 4,
                   obs = list(batch = rep(c("A", "B"), 3)))

  expect_equal(.h5ad_index_storage(tmp), "dataset")
})

test_that(".read_h5ad_metadata extracts obs", {
  skip_if_no_py_anndata()

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 10, n_genes = 5, obs = list(
    batch     = rep(c("A", "B"), each = 5),
    condition = rep(c("ctrl", "stim"), 5)
  ))

  result <- .read_h5ad_metadata(tmp, verbose = FALSE)

  expect_false(is.null(result$obs))
  expect_s3_class(result$obs, "data.frame")
  expect_equal(nrow(result$obs), 10L)
  expect_true(all(c("batch", "condition") %in% colnames(result$obs)))
  expect_equal(sort(unique(as.character(result$obs$batch))), c("A", "B"))
})

test_that(".read_h5ad_metadata extracts obsm", {
  skip_if_no_py_anndata()

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 10, n_genes = 5,
                   obsm = list(X_umap = matrix(rnorm(20), nrow = 10, ncol = 2)))

  result <- .read_h5ad_metadata(tmp, verbose = FALSE)

  expect_false(is.null(result$obsm))
  expect_true("X_umap" %in% names(result$obsm))
  expect_equal(nrow(result$obsm$X_umap), 10L)
  expect_equal(ncol(result$obsm$X_umap), 2L)
})

test_that("load_data from h5ad attaches obs metadata to the Seurat object", {
  skip_if_no_py_anndata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 20, n_genes = 30, obs = list(
    sample_id = rep(c("S1", "S2"), each = 10),
    treatment = rep(c("ctrl", "drug"), 10)
  ))

  obj <- load_data(rna_path = tmp, project_name = "H5AD_Test", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 20L)
  expect_true("sample_id" %in% colnames(obj@meta.data),
              info = "obs column 'sample_id' should reach @meta.data")
  expect_true("treatment" %in% colnames(obj@meta.data),
              info = "obs column 'treatment' should reach @meta.data")

  expect_equal(sort(unique(as.character(obj$sample_id))), c("S1", "S2"))
})

test_that("load_data from h5ad attaches obsm embeddings as reductions", {
  skip_if_no_py_anndata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 20, n_genes = 30, obsm = list(
    X_umap = matrix(rnorm(40), nrow = 20, ncol = 2),
    X_pca  = matrix(rnorm(200), nrow = 20, ncol = 10)
  ))

  obj <- load_data(rna_path = tmp, project_name = "H5AD_Emb_Test", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_true("X_umap" %in% names(obj@reductions))
  expect_true("X_pca"  %in% names(obj@reductions))
  expect_equal(ncol(Embeddings(obj, "X_umap")), 2L)
  expect_equal(ncol(Embeddings(obj, "X_pca")), 10L)
})

test_that("load_data from h5ad without obs/obsm still works", {
  skip_if_no_py_anndata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))

  .write_test_h5ad(tmp, n_cells = 10, n_genes = 12)

  obj <- load_data(rna_path = tmp, project_name = "H5AD_Bare_Test", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 10L)
  expect_true("orig.ident" %in% colnames(obj@meta.data))
})

test_that("Seurat auto-generated columns are not overwritten by obs", {
  skip_if_no_py_anndata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))

  .write_test_h5ad(tmp, n_cells = 10, n_genes = 12, obs = list(
    orig.ident = rep("SHOULD_NOT_OVERWRITE", 10),
    custom_col = 1:10
  ))

  obj <- load_data(rna_path = tmp, project_name = "ProtectTest", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(unique(as.character(obj$orig.ident)), "ProtectTest")
  expect_true("custom_col" %in% colnames(obj@meta.data))
})

test_that("load_data from h5mu builds RNA and ADT assays and keeps obsm", {
  skip_if_no_py_mudata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5mu"); on.exit(unlink(tmp))
  .write_test_h5mu(tmp, n_cells = 12L, n_genes = 20L, n_adt = 5L,
                   obsm = list(X_umap = matrix(rnorm(24), nrow = 12, ncol = 2)))

  obj <- load_data(rna_path = tmp, project_name = "H5MU_Test", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 12L)
  expect_true("ADT" %in% names(obj@assays))
  expect_equal(nrow(obj[["ADT"]]), 5L)
  expect_true("X_umap" %in% names(obj@reductions))

  expect_false(any(c("attrname", "axes", "dim", "is_view") %in%
                   names(obj@reductions)))
})

test_that("load_data recovers an ADT matrix stored in h5ad obsm", {
  skip_if_no_py_anndata()
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad"); on.exit(unlink(tmp))
  .write_test_h5ad(tmp, n_cells = 12, n_genes = 20,
                   obsm = list(protein = matrix(rpois(60, 20), nrow = 12, ncol = 5)))

  obj <- load_data(rna_path = tmp, project_name = "ADT_Obsm_Test", verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_true("ADT" %in% names(obj@assays))
  expect_equal(nrow(obj[["ADT"]]), 5L)
  expect_equal(ncol(obj[["ADT"]]), 12L)
})
