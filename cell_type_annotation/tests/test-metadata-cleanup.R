skip_mc <- function() {
  if (!exists(".drop_working_columns", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
  skip_if_not_installed("Seurat")
}

.make_mc_so <- function(prefix = NULL) {
  m <- matrix(rpois(50, 2), nrow = 5,
              dimnames = list(paste0("g", 1:5), paste0("c", 1:10)))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(m, "sparseMatrix")))
  p <- function(x) if (!is.null(prefix)) paste0(prefix, "_", x) else x
  so@meta.data[[p("SingleR_main")]]                <- "T cells"
  so@meta.data[[p("SingleR_main_score")]]          <- 0.8
  so@meta.data[[p("SingleR_main_delta")]]          <- 0.05
  so@meta.data[[p("SingleR_fine")]]                <- "Naive CD4 T cells"
  so@meta.data[[p("SingleR_fine_score")]]          <- 0.7
  so@meta.data[[p("SingleR_fine_delta")]]          <- 0.04
  so@meta.data[[p("SingleR_fine_cluster")]]        <- "Th2 cells"
  so@meta.data[[p("SingleR_fine_cluster_score")]]  <- 0.5
  so@meta.data[[p("SingleR_fine_cluster_delta")]]  <- 0.1
  so@meta.data[[p("Cluster")]]                     <- "3"
  so@meta.data[[p("Prediction")]]                  <- "CD4 Tn"
  so@meta.data[[p("Uncertainty_score")]]           <- 0.2
  so@meta.data[[p("vocab_rescue")]]                <- NA_character_
  so@meta.data[[p("vocab_rescue_score")]]          <- NA_real_
  so@meta.data$seurat_clusters                     <- factor("3")
  so@meta.data$RNA_clusters                        <- "3"
  so
}

test_that("working columns are dropped; client-facing columns survive", {
  skip_mc()
  res <- .drop_working_columns(.make_mc_so(), verbose = FALSE)
  cols <- colnames(res@meta.data)
  for (gone in c("SingleR_fine_cluster", "SingleR_fine_cluster_score",
                 "SingleR_fine_cluster_delta", "Cluster", "vocab_rescue_score",
                 "seurat_clusters", "SingleR_main_delta", "SingleR_fine_delta"))
    expect_false(gone %in% cols, info = gone)
  for (kept in c("SingleR_main", "SingleR_main_score", "SingleR_fine",
                 "SingleR_fine_score", "Prediction",
                 "Uncertainty_score", "vocab_rescue", "RNA_clusters"))
    expect_true(kept %in% cols, info = kept)
})

test_that("prefixed working columns are dropped too", {
  skip_mc()
  res <- .drop_working_columns(.make_mc_so(prefix = "run2"), prefix = "run2",
                               verbose = FALSE)
  cols <- colnames(res@meta.data)
  expect_false(any(c("run2_SingleR_fine_cluster", "run2_Cluster",
                     "run2_vocab_rescue_score", "seurat_clusters") %in% cols))
  expect_true(all(c("run2_SingleR_fine", "run2_Prediction", "run2_vocab_rescue") %in% cols))
})

test_that("the configured cluster_column is never dropped, even seurat_clusters", {
  skip_mc()

  res <- .drop_working_columns(.make_mc_so(), cluster_column = "seurat_clusters",
                               verbose = FALSE)
  cols <- colnames(res@meta.data)
  expect_true("seurat_clusters" %in% cols)

  for (gone in c("SingleR_fine_cluster", "Cluster", "vocab_rescue_score",
                 "SingleR_main_delta", "SingleR_fine_delta"))
    expect_false(gone %in% cols, info = gone)

  res2 <- .drop_working_columns(.make_mc_so(), cluster_column = "RNA_clusters",
                                verbose = FALSE)
  expect_false("seurat_clusters" %in% colnames(res2@meta.data))
  expect_true("RNA_clusters" %in% colnames(res2@meta.data))
})

test_that("object without any working columns passes through unchanged", {
  skip_mc()
  so <- .make_mc_so()
  for (dc in c("SingleR_fine_cluster", "SingleR_fine_cluster_score",
               "SingleR_fine_cluster_delta", "Cluster", "vocab_rescue_score",
               "SingleR_main_delta", "SingleR_fine_delta"))
    so@meta.data[[dc]] <- NULL
  so@meta.data$seurat_clusters <- NULL
  before <- colnames(so@meta.data)
  res <- .drop_working_columns(so, verbose = FALSE)
  expect_identical(colnames(res@meta.data), before)
})
