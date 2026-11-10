skip_assign <- function() {
  if (!exists(".assign_dominant_states", mode = "function"))
    skip("classify_functional_state.R not sourced")
}

mat3 <- function() {
  m <- matrix(c(0.90, 0.10, 0.05,
                0.20, 0.25, 0.05,
                0.02, 0.01, 0.03),
              nrow = 3, byrow = TRUE,
              dimnames = list(paste0("c", 1:3), c("alpha", "beta", "gamma")))
  m
}

test_that("which.max labelling and margin are correct with min_score = NULL", {
  skip_assign()
  a <- .assign_dominant_states(mat3(), min_score = NULL)
  expect_identical(a$state, c("alpha", "beta", "gamma"))
  expect_equal(a$score,  c(0.90, 0.25, 0.03))
  expect_equal(a$margin, c(0.90 - 0.10, 0.25 - 0.20, 0.03 - 0.02))
})

test_that("min_score relabels weak winners as 'unassigned', keeping score and margin", {
  skip_assign()
  a <- .assign_dominant_states(mat3(), min_score = 0.15)
  expect_identical(a$state, c("alpha", "beta", "unassigned"))

  expect_equal(a$score[3],  0.03)
  expect_equal(a$margin[3], 0.01)
})

test_that("min_score = 0 labels every cell (score >= 0 always)", {
  skip_assign()
  a <- .assign_dominant_states(mat3(), min_score = 0)
  expect_false(any(a$state == "unassigned"))
})

test_that("a cell exactly AT the threshold keeps its label (strict less-than)", {
  skip_assign()
  a <- .assign_dominant_states(mat3(), min_score = 0.25)
  expect_identical(a$state[2], "beta")
})

test_that("ties go to the first column, matching historical which.max behaviour", {
  skip_assign()
  m <- matrix(c(0.5, 0.5), nrow = 1,
              dimnames = list("c1", c("first", "second")))
  a <- .assign_dominant_states(m)
  expect_identical(a$state, "first")
  expect_equal(a$margin, 0)
})

test_that("single-signature matrix works and margin equals the score", {
  skip_assign()
  m <- matrix(c(0.4, 0.1), ncol = 1, dimnames = list(c("c1", "c2"), "only"))
  a <- .assign_dominant_states(m, min_score = 0.2)
  expect_identical(a$state, c("only", "unassigned"))
  expect_equal(a$margin, c(0.4, 0.1))
})

test_that("a state literally named 'unassigned' is rejected when min_score is set", {
  skip_assign()
  m <- matrix(0.5, nrow = 1, ncol = 2,
              dimnames = list("c1", c("unassigned", "other")))
  expect_error(.assign_dominant_states(m, min_score = 0.1), "unassigned")

  expect_silent(.assign_dominant_states(m, min_score = NULL))
})

skip_fsp_ms <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")
  if (!exists("classify_functional_states_ucell", mode = "function"))
    skip("classify_functional_state.R not sourced")
}

test_that("classifier validates min_score input", {
  skip_fsp_ms()
  obj <- make_fsp_seurat()
  for (bad in list(-0.1, 1.5, "high", c(0.1, 0.2), NA_real_)) {
    expect_error(
      classify_functional_states_ucell(obj, TEST_MARKERS_JSON,
                                       min_score = bad, verbose = FALSE),
      "min_score")
  }
})

test_that("min_score = NULL preserves historical labels and adds the margin column", {
  skip_fsp_ms()
  obj <- make_fsp_seurat()
  out <- classify_functional_states_ucell(obj, TEST_MARKERS_JSON,
                                          use_hvg = FALSE, verbose = FALSE)
  md <- out@meta.data
  expect_true("dominant_state_margin" %in% colnames(md))
  expect_false(any(md$dominant_state == "unassigned"))

  states <- c("exhaustion", "naive-like", "cytotoxicity")
  sm <- as.matrix(md[, states])
  top2 <- t(apply(sm, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
  expect_equal(unname(md$dominant_state_margin), unname(top2[, 1] - top2[, 2]),
               tolerance = 1e-12)
  expect_true(all(md$dominant_state_margin >= 0))
})

test_that("an unreachably high min_score sends every cell to 'unassigned'", {
  skip_fsp_ms()
  obj <- make_fsp_seurat()
  out <- classify_functional_states_ucell(obj, TEST_MARKERS_JSON,
                                          min_score = 0.999999,
                                          use_hvg = FALSE, verbose = FALSE)
  expect_true(all(out$dominant_state == "unassigned"))

  expect_true(all(out$dominant_state_score >= 0))
})

test_that("validate_classification tolerates 'unassigned' labels", {
  skip_fsp_ms()
  obj <- make_fsp_seurat()

  ref <- classify_functional_states_ucell(obj, TEST_MARKERS_JSON,
                                          use_hvg = FALSE, verbose = FALSE)
  sc  <- ref$dominant_state_score
  if (min(sc) == max(sc)) skip("degenerate fixture: all dominant scores equal")
  thr <- (min(sc) + max(sc)) / 2

  out <- classify_functional_states_ucell(obj, TEST_MARKERS_JSON,
                                          min_score = thr,
                                          use_hvg = FALSE, verbose = FALSE)
  expect_true(any(out$dominant_state == "unassigned"))
  expect_true(any(out$dominant_state != "unassigned"))

  expect_no_error(capture.output(validate_classification(out)))
})
