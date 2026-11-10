test_that("validate_classification passes for valid scores", {
  obj <- make_fsp_seurat()
  obj$exhaustion   <- runif(ncol(obj), 0, 0.8)
  obj$cytotoxicity <- runif(ncol(obj), 0, 0.6)
  obj$dominant_state <- ifelse(obj$exhaustion > obj$cytotoxicity,
                                "exhaustion", "cytotoxicity")

  result <- validate_classification(obj, states = c("exhaustion", "cytotoxicity"))
  expect_true(result)
})

test_that("validate_classification warns for out-of-range scores", {
  obj <- make_fsp_seurat()
  obj$exhaustion   <- runif(ncol(obj), -0.1, 1.1)
  obj$cytotoxicity <- runif(ncol(obj), 0, 0.6)
  obj$dominant_state <- "exhaustion"

  result <- validate_classification(obj, states = c("exhaustion", "cytotoxicity"))
  expect_false(result)
})

test_that("validate_classification errors on missing state columns", {
  obj <- make_fsp_seurat()
  expect_error(
    validate_classification(obj, states = c("nonexistent_state")),
    regexp = "not found"
  )
})

test_that("validate_classification auto-detects states from dominant_state", {
  obj <- make_fsp_seurat()
  obj$exhaustion   <- runif(ncol(obj), 0, 0.5)
  obj$cytotoxicity <- runif(ncol(obj), 0, 0.5)
  obj$dominant_state <- rep(c("exhaustion", "cytotoxicity"), length.out = ncol(obj))

  result <- validate_classification(obj)
  expect_true(is.logical(result))
})

test_that("validate_classification errors when no dominant_state and no states arg", {
  obj <- make_fsp_seurat()

  expect_error(
    validate_classification(obj),
    regexp = "dominant_state|detect"
  )
})
