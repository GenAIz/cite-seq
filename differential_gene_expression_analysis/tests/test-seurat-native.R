test_that("Seurat wilcox FindAllMarkers recovers the cluster DE genes", {

  obj <- make_mini_dge_seurat_clustered()

  result <- perform_dge_analysis(
    seurat_obj  = obj,
    test.use    = "wilcox",
    cluster_col = "CellType",
    verbose     = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("gene" %in% colnames(result))

  expect_true(any(result$gene %in% sprintf("GENE%04d", 21:40)))
})

test_that("FindMarkers with ident.1/ident.2 resolves Idents from condition_col", {

  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj    = obj,
    test.use      = "wilcox",
    cluster_col   = "CellType",
    condition_col = "Condition",
    ident.1       = "Treated",
    ident.2       = "Control",
    verbose       = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("gene" %in% colnames(result))

  expect_true(any(result$gene %in% sprintf("GENE%04d", 1:20)))
})

test_that("FindMarkers with ident.1 alone compares the group against the rest", {
  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj    = obj,
    test.use      = "wilcox",
    cluster_col   = "CellType",
    condition_col = "Condition",
    ident.1       = "Treated",
    verbose       = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true(any(result$gene %in% sprintf("GENE%04d", 1:20)))
})

test_that("FindMarkers ident.1/ident.2 falls back to cluster_col for cluster levels", {

  obj <- make_mini_dge_seurat_clustered()

  result <- perform_dge_analysis(
    seurat_obj  = obj,
    test.use    = "wilcox",
    cluster_col = "CellType",
    ident.1     = "Cluster1",
    ident.2     = "Cluster2",
    verbose     = FALSE
  )

  expect_gt(nrow(result), 0)
  expect_true(any(result$gene %in% sprintf("GENE%04d", 21:40)))
})

test_that("deliberately pre-set Idents containing the requested levels are respected", {

  obj <- make_mini_dge_seurat()
  Idents(obj) <- "Condition"

  result <- perform_dge_analysis(
    seurat_obj  = obj,
    test.use    = "wilcox",
    cluster_col = "CellType",

    ident.1     = "Treated",
    ident.2     = "Control",
    verbose     = FALSE
  )

  expect_gt(nrow(result), 0)
  expect_true(any(result$gene %in% sprintf("GENE%04d", 1:20)))
})

test_that("null test.use defaults to MAST", {
  skip_if_not_installed("MAST")

  obj <- make_mini_dge_seurat()

  expect_message(
    result <- perform_dge_analysis(
      seurat_obj    = obj,
      test.use      = NULL,
      cluster_col   = "CellType",
      condition_col = "Condition",
      ident.1       = "Treated",
      ident.2       = "Control",
      verbose       = TRUE
    ),
    regexp = "test.use = MAST"
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})
