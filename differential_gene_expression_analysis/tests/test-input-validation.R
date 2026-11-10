test_that("invalid test.use throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "nonexistent",
                         cluster_col = "CellType", verbose = FALSE),
    regexp = "Invalid test.use"
  )
})

test_that("bulk method without pseudobulking throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "limma",
                         pseudobulking = FALSE, cluster_col = "CellType",
                         verbose = FALSE),
    regexp = "pseudobulking"
  )
})

test_that("missing cluster column throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "wilcox",
                         cluster_col = "nonexistent_col", verbose = FALSE),
    regexp = "Missing|not found"
  )
})

test_that("missing condition column throws an error for pseudobulk", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "limma",
                         pseudobulking = TRUE,
                         cluster_col = "CellType",
                         sample_col = "SampleID",
                         condition_col = "nonexistent_condition",
                         verbose = FALSE),
    regexp = "Missing|not found"
  )
})

test_that("invalid design_strategy throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "limma",
                         pseudobulking = TRUE,
                         design_strategy = "nonexistent",
                         cluster_col = "CellType",
                         sample_col = "SampleID",
                         condition_col = "Condition",
                         verbose = FALSE),
    regexp = "design_strategy|arg"
  )
})

test_that("invalid ident.1 throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "limma",
                         pseudobulking = TRUE,
                         cluster_col = "CellType",
                         sample_col = "SampleID",
                         condition_col = "Condition",
                         ident.1 = "NONEXISTENT_GROUP",
                         ident.2 = "Control",
                         verbose = FALSE),
    regexp = "not found|ident.1"
  )
})

test_that("invalid reference_condition throws an error", {
  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(seurat_obj = obj, test.use = "limma",
                         pseudobulking = TRUE,
                         cluster_col = "CellType",
                         sample_col = "SampleID",
                         condition_col = "Condition",
                         reference_condition = "NONEXISTENT",
                         verbose = FALSE),
    regexp = "not found|reference_condition"
  )
})
