skip_aba <- function() {
  if (!exists(".kbet_rejection", mode = "function"))
    skip("process_seurat_pipeline.R not sourced (Seurat unavailable at runner load)")
}

test_that("kBET rejection is high for separated batches, low for mixed", {
  skip_aba()
  set.seed(1)
  n <- 400

  emb_sep <- rbind(matrix(rnorm(n), ncol = 2),
                   matrix(rnorm(n, mean = 30), ncol = 2))

  emb_mix <- rbind(matrix(rnorm(n), ncol = 2),
                   matrix(rnorm(n), ncol = 2))
  batch <- rep(c("b1", "b2"), each = n / 2)

  r_sep <- .kbet_rejection(emb_sep, batch, k = 15L, n_sample = 400L)
  r_mix <- .kbet_rejection(emb_mix, batch, k = 15L, n_sample = 400L)
  expect_gt(r_sep, 0.9)
  expect_lt(r_mix, 0.3)
  expect_gt(r_sep, r_mix)
})

test_that("kBET is deterministic for a seed and NA for degenerate input", {
  skip_aba()
  set.seed(2)
  emb <- matrix(rnorm(600), ncol = 3)
  batch <- sample(c("x", "y"), 200, replace = TRUE)
  expect_identical(.kbet_rejection(emb, batch, n_sample = 100L, seed = 7L),
                   .kbet_rejection(emb, batch, n_sample = 100L, seed = 7L))

  expect_true(is.na(.kbet_rejection(emb, rep("x", 200))))

  expect_true(is.na(.kbet_rejection(emb[1:10, ], batch[1:10], k = 30L)))

  b2 <- batch; b2[1:50] <- NA
  expect_true(is.finite(.kbet_rejection(emb, b2, n_sample = 100L)))
})

.make_adt_so <- function(n = 120, n_markers = 12, batch_shift = 0) {
  set.seed(5)
  m <- matrix(rnorm(n_markers * n), nrow = n_markers,
              dimnames = list(paste0("Ab", seq_len(n_markers)),
                              paste0("c", seq_len(n))))
  if (batch_shift != 0) m[, seq_len(n / 2)] <- m[, seq_len(n / 2)] + batch_shift
  counts <- matrix(rpois(5 * n, 2), nrow = 5,
                   dimnames = list(paste0("g", 1:5), colnames(m)))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(counts, "sparseMatrix")))
  adt <- suppressWarnings(Seurat::CreateAssayObject(counts = as(exp(m), "sparseMatrix")))
  so[["ADT"]] <- adt
  so <- Seurat::SetAssayData(so, assay = "ADT", layer = "data",
                             new.data = as(m, "sparseMatrix"))
  so@meta.data$batch <- rep(c("b1", "b2"), each = n / 2)
  so
}

test_that("provisional ADT PCA returns well-formed scores", {
  skip_aba()
  skip_if_not_installed("Seurat")
  so <- .make_adt_so()
  qc <- .adt_qc_pca_scores(so, dims = 5L, n_sample = 100L)
  expect_false(is.null(qc))
  expect_true(is.matrix(qc$emb))
  expect_lte(ncol(qc$emb), 5L)
  expect_identical(nrow(qc$emb), length(qc$idx))
  expect_true(all(is.finite(qc$emb)))
})

test_that("provisional PCA + kBET catches a batch-shifted ADT panel", {
  skip_aba()
  skip_if_not_installed("Seurat")
  so_shift <- .make_adt_so(batch_shift = 8)
  so_clean <- .make_adt_so(batch_shift = 0)
  k_shift <- {
    qc <- .adt_qc_pca_scores(so_shift, n_sample = 120L)
    .kbet_rejection(qc$emb, so_shift@meta.data$batch[qc$idx], k = 15L, n_sample = 120L)
  }
  k_clean <- {
    qc <- .adt_qc_pca_scores(so_clean, n_sample = 120L)
    .kbet_rejection(qc$emb, so_clean@meta.data$batch[qc$idx], k = 15L, n_sample = 120L)
  }
  expect_gt(k_shift, 0.8)
  expect_lt(k_clean, 0.3)
})

test_that("guardrail adopts only on sufficient improvement", {
  skip_aba()
  expect_true(.adt_guardrail_adopt(0.90, 0.40, 0.05))
  expect_true(.adt_guardrail_adopt(0.50, 0.45, 0.05))
  expect_false(.adt_guardrail_adopt(0.50, 0.47, 0.05))
  expect_false(.adt_guardrail_adopt(0.50, 0.60, 0.05))

  expect_false(.adt_guardrail_adopt(NA, 0.3, 0.05))
  expect_false(.adt_guardrail_adopt(0.9, NA, 0.05))
  expect_false(.adt_guardrail_adopt(NA, NA, 0.05))
})

test_that("process_and_sketch exposes the new ADT batch parameters", {
  skip_aba()
  fm <- formals(process_and_sketch)
  expect_identical(fm$adt_correction_method, "harmony")
  expect_identical(eval(fm$adt_kbet_threshold), 0.6)
  expect_identical(eval(fm$adt_guardrail_min_improvement), 0.05)
})

test_that(".detect_adt_batch reports the kBET fields", {
  skip_aba()
  skip_if_not_installed("Seurat")
  so <- .make_adt_so(batch_shift = 8)
  res <- .detect_adt_batch(so, batch_column = "batch", n_batches = 2L,
                           wasserstein_threshold = 0.3, verbose = FALSE,
                           kbet_threshold = 0.6)
  expect_true(all(c("kbet_rejection", "kbet_threshold", "kbet_flagged",
                    "correction_recommended") %in% names(res)))
  expect_true(is.finite(res$kbet_rejection))
  expect_gt(res$kbet_rejection, 0.6)
  expect_true(res$kbet_flagged)
  expect_true(res$correction_recommended)
})

test_that(".detect_adt_batch does not flag a clean panel via kBET", {
  skip_aba()
  skip_if_not_installed("Seurat")
  so <- .make_adt_so(batch_shift = 0)
  res <- .detect_adt_batch(so, batch_column = "batch", n_batches = 2L,
                           wasserstein_threshold = 0.3, verbose = FALSE,
                           kbet_threshold = 0.6)
  expect_false(isTRUE(res$kbet_flagged))
})

test_that("harmony corrector errors clearly when harmony is unavailable", {
  skip_aba()
  skip_if(requireNamespace("harmony", quietly = TRUE),
          "harmony installed - dependency error path not reachable")
  so <- .make_adt_so()
  expect_error(
    .correct_adt_batch_harmony(so, batch_column = "batch",
                               theta = 2, niter = 10, verbose = FALSE),
    "harmony"
  )
})

.aba_pipeline_src <- function() {
  cand <- c("script/process_seurat_pipeline.R",
            "../script/process_seurat_pipeline.R")
  path <- cand[file.exists(cand)][1]
  if (is.na(path)) return(NA_character_)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("the ADTnorm guardrail snapshots and restores the whole ADT assay", {
  src <- .aba_pipeline_src()
  skip_if(is.na(src), "process_seurat_pipeline.R not found from test wd")
  expect_true(grepl("adt_assay_backup <- seurat_object\\[\\[\"ADT\"\\]\\]", src),
              info = "stage 4b must snapshot the ADT assay object")
  expect_true(grepl("seurat_object\\[\\[\"ADT\"\\]\\] <- adt_assay_backup", src),
              info = "the revert must restore the whole assay, not just a layer")

  expect_false(grepl("adt_data_backup", src, fixed = TRUE))
})

test_that("forced ADT harmony records a failure when the batch column is absent", {
  src <- .aba_pipeline_src()
  skip_if(is.na(src), "process_seurat_pipeline.R not found from test wd")

  expect_true(grepl("ADT Harmony correction SKIPPED", src, fixed = TRUE))
  expect_true(grepl("not found in metadata", src, fixed = TRUE))
})
