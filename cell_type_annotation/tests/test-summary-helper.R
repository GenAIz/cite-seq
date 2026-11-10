test_that("summary prints without error for SingleR-annotated object", {
  obj <- make_annotatable_seurat()
  obj$SingleR_main <- rep(c("T cells", "B cells"), length.out = ncol(obj))
  obj$SingleR_fine <- rep(c("CD4+", "CD8+", "Naive B", "Memory B"), length.out = ncol(obj))

  expect_message(
    .print_annotation_summary(obj, method = "SingleR", prefix = NULL),
    regexp = "SingleR main"
  )
})

test_that("summary prints without error for STCAT-annotated object", {
  obj <- make_annotatable_seurat()
  obj$Prediction <- rep(c("CD4_T", "CD8_T"), length.out = ncol(obj))

  expect_message(
    .print_annotation_summary(obj, method = "STCAT", prefix = NULL),
    regexp = "STCAT predictions"
  )
})

test_that("summary handles prefixed columns correctly", {
  obj <- make_annotatable_seurat()
  obj$run1_SingleR_main <- rep("T cells", ncol(obj))
  obj$run1_Prediction   <- rep("CD4_T", ncol(obj))

  msgs <- capture_messages(
    .print_annotation_summary(obj, method = "both", prefix = "run1")
  )
  expect_true(any(grepl("SingleR main", msgs)))
  expect_true(any(grepl("STCAT predictions", msgs)))
})

test_that("summary handles missing annotation columns gracefully", {
  obj <- make_annotatable_seurat()

  expect_message(
    .print_annotation_summary(obj, method = "both", prefix = NULL),
    regexp = "==="
  )
})
