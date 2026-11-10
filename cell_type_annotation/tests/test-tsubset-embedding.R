skip_tse <- function() {
  if (!exists(".compute_tcell_subset_embedding", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
  skip_if_not_installed("Seurat")
}

.make_tse_so <- function(n_a = 12, n_b = 12, batches = NULL, seed = 7) {
  set.seed(seed)
  n <- n_a + n_b
  m <- matrix(rpois(100 * n, 1), nrow = 100,
              dimnames = list(paste0("g", 1:100), paste0("c", seq_len(n))))
  m[1:40,  seq_len(n_a)]        <- m[1:40,  seq_len(n_a)]        + rpois(40 * n_a, 8)
  m[41:80, n_a + seq_len(n_b)]  <- m[41:80, n_a + seq_len(n_b)]  + rpois(40 * n_b, 8)
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(m, "sparseMatrix")))
  so <- suppressWarnings(Seurat::NormalizeData(so, verbose = FALSE))
  if (!is.null(batches)) so@meta.data$batch <- batches
  so
}

test_that("subset re-embedding returns a well-formed matrix for the gated cells", {
  skip_tse()
  so    <- .make_tse_so()
  cells <- colnames(so)[1:20]
  emb   <- suppressWarnings(.compute_tcell_subset_embedding(
    so, tcell_barcodes = cells, assay = "RNA",
    batch_column = NULL, dims = 10L, verbose = FALSE))
  expect_true(is.matrix(emb))
  expect_identical(sort(rownames(emb)), sort(cells))
  expect_lte(ncol(emb), 10L)
  expect_gte(ncol(emb), 2L)
  expect_true(all(is.finite(emb)))
})

test_that("npcs is capped below the cell count for tiny subsets", {
  skip_tse()
  so    <- .make_tse_so(n_a = 4, n_b = 4)
  cells <- colnames(so)
  emb   <- suppressWarnings(.compute_tcell_subset_embedding(
    so, tcell_barcodes = cells, assay = "RNA",
    batch_column = NULL, dims = 30L, verbose = FALSE))
  expect_true(is.matrix(emb))
  expect_lt(ncol(emb), length(cells))
})

test_that("the subset PCA separates the two expression programs", {
  skip_tse()
  so    <- .make_tse_so(n_a = 12, n_b = 12)
  cells <- colnames(so)
  emb   <- suppressWarnings(.compute_tcell_subset_embedding(
    so, tcell_barcodes = cells, assay = "RNA",
    batch_column = NULL, dims = 10L, verbose = FALSE))
  pc1_a <- emb[paste0("c", 1:12), 1]
  pc1_b <- emb[paste0("c", 13:24), 1]

  expect_true(max(pc1_a) < min(pc1_b) || max(pc1_b) < min(pc1_a))
})

test_that("multi-batch without harmony warns that the embedding is uncorrected", {
  skip_tse()
  skip_if(requireNamespace("harmony", quietly = TRUE),
          "harmony installed - uncorrected-warning path not reachable")
  so    <- .make_tse_so(batches = rep(c("run1", "run2"), each = 12))
  expect_warning(
    emb <- .compute_tcell_subset_embedding(
      so, tcell_barcodes = colnames(so), assay = "RNA",
      batch_column = "batch", dims = 10L, verbose = FALSE),
    "NOT batch-corrected"
  )
  expect_true(is.matrix(emb))
})

test_that("multi-batch with harmony returns a corrected embedding silently", {
  skip_tse()
  skip_if_not_installed("harmony")
  so  <- .make_tse_so(batches = rep(c("run1", "run2"), each = 12))
  emb <- .compute_tcell_subset_embedding(
    so, tcell_barcodes = colnames(so), assay = "RNA",
    batch_column = "batch", dims = 10L, verbose = FALSE)
  expect_true(is.matrix(emb))
  expect_identical(sort(rownames(emb)), sort(colnames(so)))
})

test_that("single-level batch column skips harmony without warning", {
  skip_tse()
  so  <- .make_tse_so(batches = rep("run1", 24))
  expect_no_warning(
    emb <- .compute_tcell_subset_embedding(
      so, tcell_barcodes = colnames(so), assay = "RNA",
      batch_column = "batch", dims = 10L, verbose = FALSE))
  expect_true(is.matrix(emb))
})

test_that("failure returns NULL with a warning (fallback contract)", {
  skip_tse()
  so <- .make_tse_so()
  expect_warning(
    emb <- .compute_tcell_subset_embedding(
      so, tcell_barcodes = colnames(so), assay = "NOPE",
      batch_column = NULL, dims = 10L, verbose = FALSE),
    "re-embedding failed"
  )
  expect_null(emb)
})
