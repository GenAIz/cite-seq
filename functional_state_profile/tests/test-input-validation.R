test_that("non-Seurat input throws an error", {
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = data.frame(x = 1),
      marker_json_path = TEST_MARKERS_JSON
    ),
    regexp = "Seurat"
  )
})

test_that("missing marker JSON file throws an error", {
  obj <- make_fsp_seurat()
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = obj,
      marker_json_path = "/nonexistent/path/markers.json"
    ),
    regexp = "not found"
  )
})

test_that("invalid assay throws an error", {
  obj <- make_fsp_seurat()
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = obj,
      marker_json_path = TEST_MARKERS_JSON,
      assay = "NONEXISTENT"
    ),
    regexp = "not found|Assay"
  )
})

test_that("invalid use_hvg value throws an error", {
  obj <- make_fsp_seurat()
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = obj,
      marker_json_path = TEST_MARKERS_JSON,
      use_hvg = 42
    ),
    regexp = "use_hvg"
  )
})

test_that("non-positive maxRank throws an error", {
  obj <- make_fsp_seurat()
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = obj,
      marker_json_path = TEST_MARKERS_JSON,
      maxRank = -1
    ),
    regexp = "maxRank|positive"
  )
})

test_that("zero maxRank throws an error", {
  obj <- make_fsp_seurat()
  expect_error(
    classify_functional_states_ucell(
      seurat_obj = obj,
      marker_json_path = TEST_MARKERS_JSON,
      maxRank = 0
    ),
    regexp = "maxRank|positive"
  )
})
