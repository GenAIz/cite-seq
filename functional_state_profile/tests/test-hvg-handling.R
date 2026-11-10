test_that("use_hvg = FALSE uses all genes", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true("dominant_state" %in% colnames(result@meta.data))
})

test_that("use_hvg = TRUE without HVGs falls back with warning", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  obj[["RNA"]]@meta.data <- data.frame(row.names = rownames(obj[["RNA"]]))

  expect_warning(
    classify_functional_states_ucell(
      seurat_obj       = obj,
      marker_json_path = TEST_MARKERS_JSON,
      use_hvg          = TRUE,
      maxRank          = 500,
      ncores           = 1,
      verbose          = FALSE
    ),
    regexp = "HVG|FindVariableFeatures|Falling back"
  )
})

test_that("use_hvg = 'auto' works silently when HVGs exist", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()
  suppressWarnings({

    obj <- FindVariableFeatures(obj, nfeatures = nrow(obj), verbose = FALSE)
  })

  expect_no_warning(
    result <- classify_functional_states_ucell(
      seurat_obj       = obj,
      marker_json_path = TEST_MARKERS_JSON,
      use_hvg          = "auto",
      maxRank          = 500,
      ncores           = 1,
      verbose          = FALSE
    )
  )

  expect_s4_class(result, "Seurat")
})

test_that("use_hvg = 'auto' works silently when HVGs do not exist", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  expect_no_warning(
    result <- classify_functional_states_ucell(
      seurat_obj       = obj,
      marker_json_path = TEST_MARKERS_JSON,
      use_hvg          = "auto",
      maxRank          = 500,
      ncores           = 1,
      verbose          = FALSE
    )
  )

  expect_s4_class(result, "Seurat")
})
