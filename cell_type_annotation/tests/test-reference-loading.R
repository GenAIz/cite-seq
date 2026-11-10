test_that("unknown reference name throws an informative error", {
  expect_error(
    .load_reference_dataset("nonexistent_reference"),
    regexp = "Unknown reference|monaco|hpca"
  )
})

test_that("custom SCE object is returned as-is", {
  skip_if_not_installed("SingleCellExperiment")
  suppressPackageStartupMessages(library(SingleCellExperiment))

  mat <- matrix(rpois(100, 5), nrow = 10, ncol = 10,
                dimnames = list(paste0("G", 1:10), paste0("C", 1:10)))
  sce <- SingleCellExperiment(assays = list(logcounts = mat))
  colData(sce)$label.main <- rep(c("TypeA", "TypeB"), each = 5)

  result <- .load_reference_dataset(sce)
  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(ncol(result), 10L)
})

test_that("monaco reference loads successfully", {
  skip_if_not_installed("celldex")
  skip_on_ci()

  ref <- .load_reference_dataset("monaco")
  expect_s4_class(ref, "SummarizedExperiment")
  expect_gt(ncol(ref), 0L)
})
