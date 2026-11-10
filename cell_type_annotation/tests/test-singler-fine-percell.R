skip_fpc <- function() {
  if (!exists(".run_singler_fine_percell", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")

  ok <- tryCatch({
    ref <- .make_fpc_ref()
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(logcounts = .fpc_program_matrix(c("A", "B"))))
    p <- SingleR::SingleR(test = sce, ref = ref,
                          labels = SummarizedExperiment::colData(ref)$label.fine,
                          BPPARAM = BiocParallel::SerialParam())
    TRUE
  }, error = function(e) FALSE)
  if (!ok) skip("SingleR not functional on tiny inputs in this environment")
}

.fpc_program_matrix <- function(programs, noise = 0) {
  n <- length(programs)
  m <- matrix(0.1, nrow = 20, ncol = n,
              dimnames = list(paste0("g", 1:20),
                              paste0("cell", seq_len(n))))
  for (j in seq_len(n)) {
    rows <- if (programs[j] == "A") 1:10 else 11:20
    m[rows, j] <- 5 + noise * stats::runif(10)
  }
  m
}

.make_fpc_ref <- function() {
  ref_mat <- cbind(
    .fpc_program_matrix(rep("A", 3)),
    .fpc_program_matrix(rep("B", 3))
  )
  colnames(ref_mat) <- paste0("ref", 1:6)
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(logcounts = ref_mat),
    colData = S4Vectors::DataFrame(
      label.main = c(rep("TypeA", 3), rep("TypeB", 3)),
      label.fine = c(rep("TypeA fine", 3), rep("TypeB fine", 3))
    )
  )
}

.make_fpc_so <- function(n_a = 6, n_b = 4, cluster_fine = "TypeA fine",
                         prefix = NULL) {
  mat <- cbind(.fpc_program_matrix(rep("A", n_a)),
               .fpc_program_matrix(rep("B", n_b)))
  colnames(mat) <- paste0("c", seq_len(n_a + n_b))
  so <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = as(expm1(mat), "sparseMatrix")))

  so <- Seurat::SetAssayData(so, layer = "data", new.data = as(mat, "sparseMatrix"))
  fine_col <- if (!is.null(prefix)) paste0(prefix, "_SingleR_fine") else "SingleR_fine"
  so@meta.data[[fine_col]]                   <- cluster_fine
  so@meta.data[[paste0(fine_col, "_score")]] <- 0.9
  so@meta.data[[paste0(fine_col, "_delta")]] <- 0.5
  so
}

test_that("per-cell fine labels resolve heterogeneity within a blanketed cluster", {
  skip_fpc()
  so  <- .make_fpc_so(n_a = 6, n_b = 4, cluster_fine = "TypeA fine")
  res <- .run_singler_fine_percell(so, ref_data = .make_fpc_ref(), assay = "RNA",
                                   prefix = NULL, chunk_size = 100L, verbose = FALSE)
  md <- res@meta.data
  expect_identical(unname(md$SingleR_fine[1:6]),  rep("TypeA fine", 6))
  expect_identical(unname(md$SingleR_fine[7:10]), rep("TypeB fine", 4))

  expect_gt(length(unique(md$SingleR_fine)), 1)
})

test_that("cluster-level fine columns are preserved as SingleR_fine_cluster*", {
  skip_fpc()
  so  <- .make_fpc_so(cluster_fine = "TypeA fine")
  res <- .run_singler_fine_percell(so, ref_data = .make_fpc_ref(), assay = "RNA",
                                   prefix = NULL, chunk_size = 100L, verbose = FALSE)
  md <- res@meta.data
  expect_true(all(c("SingleR_fine_cluster", "SingleR_fine_cluster_score",
                    "SingleR_fine_cluster_delta") %in% colnames(md)))
  expect_identical(unname(md$SingleR_fine_cluster), rep("TypeA fine", 10))
  expect_identical(unname(md$SingleR_fine_cluster_score), rep(0.9, 10))

  expect_true(all(is.finite(md$SingleR_fine_score)))
  expect_true(is.numeric(md$SingleR_fine_delta))
})

test_that("chunking does not change the result", {
  skip_fpc()
  so <- .make_fpc_so(n_a = 5, n_b = 5)
  r_big   <- .run_singler_fine_percell(so, .make_fpc_ref(), assay = "RNA",
                                       prefix = NULL, chunk_size = 1000L, verbose = FALSE)
  r_small <- .run_singler_fine_percell(so, .make_fpc_ref(), assay = "RNA",
                                       prefix = NULL, chunk_size = 3L, verbose = FALSE)
  expect_identical(r_big@meta.data$SingleR_fine,       r_small@meta.data$SingleR_fine)
  expect_equal(r_big@meta.data$SingleR_fine_score,     r_small@meta.data$SingleR_fine_score)
  expect_equal(r_big@meta.data$SingleR_fine_delta,     r_small@meta.data$SingleR_fine_delta)
})

test_that("prefix routes all columns", {
  skip_fpc()
  so  <- .make_fpc_so(cluster_fine = "TypeA fine", prefix = "run2")
  res <- .run_singler_fine_percell(so, .make_fpc_ref(), assay = "RNA",
                                   prefix = "run2", chunk_size = 100L, verbose = FALSE)
  md <- res@meta.data
  expect_true(all(c("run2_SingleR_fine", "run2_SingleR_fine_cluster") %in% colnames(md)))
  expect_false("SingleR_fine" %in% colnames(md))
})

test_that("reference without label.fine warns and leaves the object unchanged", {
  skip_fpc()
  ref <- .make_fpc_ref()
  SummarizedExperiment::colData(ref)$label.fine <- NULL
  so <- .make_fpc_so(cluster_fine = "TypeA fine")
  expect_warning(
    res <- .run_singler_fine_percell(so, ref, assay = "RNA",
                                     prefix = NULL, chunk_size = 100L, verbose = FALSE),
    "no 'label.fine'"
  )
  expect_identical(unname(res@meta.data$SingleR_fine), rep("TypeA fine", 10))
  expect_false("SingleR_fine_cluster" %in% colnames(res@meta.data))
})

test_that("per-cell fine labels feed the T cell gate", {
  skip_fpc()
  if (!exists(".is_tcell_label", mode = "function")) skip("gate helper not sourced")

  main <- c("NK cells", "NK cells")
  fine_percell <- c("MAIT cells", "Natural killer cells")
  mask <- .is_tcell_label(main) | .is_tcell_label(fine_percell)
  expect_identical(mask, c(TRUE, FALSE))
})
