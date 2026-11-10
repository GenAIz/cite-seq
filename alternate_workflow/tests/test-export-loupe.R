skip_loupe <- function() {
  if (!exists(".loupe_reencode_barcodes", mode = "function") ||
      !exists(".loupe_cluster_columns",   mode = "function"))
    skip("export_loupe.R not sourced")
}

test_that(".loupe_reencode_barcodes produces valid, unique 10x barcodes", {
  skip_loupe()
  bc <- c("AAACCTGAGAAAGTGG-1-gPlexA1",
          "AAACCTGAGAAAGTGG-1-gPlexA2",
          "TTTGGTTTCACTGGGC-1-gPlexA1")
  out <- .loupe_reencode_barcodes(bc)
  expect_equal(length(out), length(bc))
  expect_false(anyDuplicated(out) > 0)
  expect_true(all(grepl("^[ACGT]{16}-[0-9]+$", out)))

  expect_false(out[1] == out[2])
})

test_that(".loupe_reencode_barcodes returns NULL when re-encoding collides", {
  skip_loupe()
  bc <- c("AAACCTGAGAAAGTGG-1", "AAACCTGAGAAAGTGG-1")
  expect_null(.loupe_reencode_barcodes(bc))
})

test_that(".loupe_cluster_columns keeps categoricals, drops IDs and continuous", {
  skip_loupe()
  set.seed(1)
  meta <- data.frame(
    cluster  = factor(sample(1:3, 300, replace = TRUE)),
    celltype = sample(c("CD4", "CD8", "NK"), 300, replace = TRUE),
    clone_id = paste0("clone", 1:300),
    pct_mt   = runif(300, 0, 20),
    n_group  = sample(1:5, 300, replace = TRUE),
    stringsAsFactors = FALSE
  )
  cols <- .loupe_cluster_columns(meta, max_groups = 200L)
  expect_true(all(c("cluster", "celltype", "n_group") %in% cols))
  expect_false("clone_id" %in% cols)
  expect_false("pct_mt"   %in% cols)
})
