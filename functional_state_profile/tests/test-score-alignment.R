skip_fsp <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("UCell")
  if (!exists("classify_functional_states_ucell", mode = "function"))
    skip("classify_functional_state.R not sourced")
}

write_markers <- function(states) {
  f <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(states, auto_unbox = TRUE, pretty = TRUE), f)
  f
}

test_that("score columns are bound to the state they belong to", {
  skip_fsp()
  obj <- make_fsp_seurat()

  markers <- list(
    zeta  = list(core = c("GZMB", "PRF1", "NKG7")),
    alpha = list(core = c("PDCD1", "LAG3", "HAVCR2")),
    mu    = list(core = c("CCR7", "LEF1", "TCF7"))
  )
  f <- write_markers(markers); on.exit(unlink(f))

  out <- classify_functional_states_ucell(obj, f, use_hvg = FALSE, verbose = FALSE)

  expect_true(all(c("zeta", "alpha", "mu") %in% colnames(out@meta.data)))
  expect_true(all(c("dominant_state", "dominant_state_score") %in% colnames(out@meta.data)))

  expect_equal(as.character(out@meta.data$dominant_state[1:15]),
               rep("alpha", 15))

  expect_equal(as.character(out@meta.data$dominant_state[16:30]),
               rep("zeta", 15))
})

test_that("dominant_state_score matches the named winning column", {
  skip_fsp()
  obj <- make_fsp_seurat()
  markers <- list(
    zeta  = list(core = c("GZMB", "PRF1", "NKG7")),
    alpha = list(core = c("PDCD1", "LAG3", "HAVCR2"))
  )
  f <- write_markers(markers); on.exit(unlink(f))
  out <- classify_functional_states_ucell(obj, f, use_hvg = FALSE, verbose = FALSE)

  md <- out@meta.data
  winning <- vapply(seq_len(nrow(md)),
                    function(i) md[[as.character(md$dominant_state[i])]][i],
                    numeric(1))
  expect_equal(unname(winning), unname(md$dominant_state_score), tolerance = 1e-12)
})

test_that("a signature with no genes in the universe is a hard error", {
  skip_fsp()
  obj <- make_fsp_seurat()
  markers <- list(
    real   = list(core = c("GZMB", "PRF1", "NKG7")),
    absent = list(core = c("NOT_A_GENE_1", "NOT_A_GENE_2"))
  )
  f <- write_markers(markers); on.exit(unlink(f))

  expect_error(
    classify_functional_states_ucell(obj, f, use_hvg = FALSE, verbose = FALSE),
    "No marker genes survive"
  )
})
