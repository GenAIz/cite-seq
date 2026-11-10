test_that("non-Seurat input throws an error", {
  expect_error(
    cell_type_annotation(seurat_object = data.frame(x = 1)),
    regexp = "Seurat"
  )
})

test_that("invalid method throws an error", {
  obj <- make_annotatable_seurat()
  expect_error(
    cell_type_annotation(seurat_object = obj, method = "invalid_method"),
    regexp = "method|arg"
  )
})

test_that("missing assay throws an error", {
  obj <- make_annotatable_seurat()
  expect_error(
    cell_type_annotation(seurat_object = obj, assay = "NONEXISTENT"),
    regexp = "not found|Assay"
  )
})

test_that("missing cluster column throws an error", {
  obj <- make_annotatable_seurat()
  expect_error(
    cell_type_annotation(seurat_object = obj,
                         cluster_column = "nonexistent_clusters"),
    regexp = "not found|Cluster"
  )
})

test_that("invalid singler_label_type throws an error", {
  obj <- make_annotatable_seurat()
  expect_error(
    cell_type_annotation(seurat_object = obj,
                         method = "SingleR",
                         singler_label_type = c("invalid")),
    regexp = "label_type|main|fine"
  )
})
