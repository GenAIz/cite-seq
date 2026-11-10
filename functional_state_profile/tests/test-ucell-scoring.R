test_that("classify_functional_states_ucell adds expected metadata columns", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  expect_s4_class(result, "Seurat")

  expect_true("exhaustion" %in% colnames(result@meta.data))
  expect_true("naive-like" %in% colnames(result@meta.data))
  expect_true("cytotoxicity" %in% colnames(result@meta.data))
  expect_true("dominant_state" %in% colnames(result@meta.data))
  expect_true("dominant_state_score" %in% colnames(result@meta.data))
})

test_that("UCell scores are in 0-1 range", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  states <- c("exhaustion", "naive-like", "cytotoxicity")
  for (state in states) {
    scores <- result@meta.data[[state]]
    expect_true(all(scores >= 0 & scores <= 1),
                info = sprintf("Scores for '%s' outside 0-1 range", state))
  }
})

test_that("dominant_state is one of the state names", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  expected_states <- c("exhaustion", "naive-like", "cytotoxicity")
  dominant <- unique(result$dominant_state)
  expect_true(all(dominant %in% expected_states))
})

test_that("dominant_state_score equals max score across states for each cell", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  states <- c("exhaustion", "naive-like", "cytotoxicity")
  score_mat <- as.matrix(result@meta.data[, states])
  expected_max <- apply(score_mat, 1, max)
  expect_equal(result$dominant_state_score, expected_max, tolerance = 1e-6)
})

test_that("use_core_only = TRUE uses fewer markers", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()

  result_all <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  result_core <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = TRUE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  expect_true("dominant_state" %in% colnames(result_core@meta.data))
  expect_true("exhaustion" %in% colnames(result_core@meta.data))

  for (state in c("exhaustion", "naive-like", "cytotoxicity")) {
    expect_false(
      isTRUE(all.equal(result_all@meta.data[[state]],
                       result_core@meta.data[[state]])),
      info = sprintf("use_core_only had no effect on '%s' scores", state))
  }
})

test_that("cell count is preserved after classification", {
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")

  obj <- make_fsp_seurat()
  n_before <- ncol(obj)

  result <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = TEST_MARKERS_JSON,
    use_core_only    = FALSE,
    use_hvg          = FALSE,
    maxRank          = 500,
    ncores           = 1,
    verbose          = FALSE
  )

  expect_equal(ncol(result), n_before)
})
