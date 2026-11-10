test_that("SingleR annotation adds expected metadata columns", {
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("celldex")

  obj <- make_annotatable_seurat()

  suppressPackageStartupMessages({
    library(SingleCellExperiment)
  })
  common_genes <- intersect(rownames(obj), c("CD3E", "CD4", "CD8A", "CD14",
                                              "CD19", "MS4A1", "LYZ", "NKG7"))
  n_genes <- max(length(common_genes), 5)

  n_pad     <- max(0L, n_genes - length(common_genes))
  all_genes <- union(common_genes,
                     paste0("GENE", sprintf("%04d", seq_len(n_pad))))
  mock_mat <- matrix(rpois(length(all_genes) * 20, 5),
                     nrow = length(all_genes), ncol = 20,
                     dimnames = list(all_genes, paste0("Ref", 1:20)))
  mock_ref <- SingleCellExperiment(assays = list(logcounts = mock_mat))
  colData(mock_ref)$label.main <- rep(c("T cells", "B cells", "Monocytes", "NK cells"), each = 5)
  colData(mock_ref)$label.fine <- rep(paste0("Sub", 1:10), each = 2)

  result <- .run_singler_annotation(
    seurat_object  = obj,
    ref_data       = mock_ref,
    cluster_column = "RNA_clusters",
    assay          = "RNA",
    label_types    = c("main", "fine"),
    workers        = 1,
    prefix         = NULL,
    verbose        = FALSE
  )

  expect_s4_class(result, "Seurat")
  md <- result@meta.data

  for (col in c("SingleR_main", "SingleR_fine",
                "SingleR_main_score", "SingleR_fine_score",
                "SingleR_main_delta", "SingleR_fine_delta")) {
    expect_true(col %in% colnames(md),
                info = sprintf("column '%s' missing", col))
  }
  expect_true(is.numeric(md$SingleR_main_score) &&
              all(is.finite(md$SingleR_main_score)))
  expect_true(is.numeric(md$SingleR_main_delta))

  for (col in c("SingleR_main", "SingleR_main_score", "SingleR_main_delta")) {
    n_per_cluster <- tapply(md[[col]], md$RNA_clusters,
                            function(x) length(unique(x)))
    expect_true(all(n_per_cluster == 1L),
                info = sprintf("'%s' varies within a cluster", col))
  }

  expect_equal(ncol(result), ncol(obj))
})

test_that("SingleR annotation with prefix names columns correctly", {
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")

  obj <- make_annotatable_seurat()

  suppressPackageStartupMessages(library(SingleCellExperiment))
  genes <- rownames(obj)[1:20]
  mock_mat <- matrix(rpois(length(genes) * 10, 5),
                     nrow = length(genes), ncol = 10,
                     dimnames = list(genes, paste0("R", 1:10)))
  mock_ref <- SingleCellExperiment(assays = list(logcounts = mock_mat))
  colData(mock_ref)$label.main <- rep(c("A", "B"), each = 5)

  result <- .run_singler_annotation(
    seurat_object  = obj,
    ref_data       = mock_ref,
    cluster_column = "RNA_clusters",
    assay          = "RNA",
    label_types    = "main",
    workers        = 1,
    prefix         = "myprefix",
    verbose        = FALSE
  )

  expect_true("myprefix_SingleR_main" %in% colnames(result@meta.data))
  expect_false("SingleR_main" %in% colnames(result@meta.data))
})

test_that("SingleR fails gracefully when zero common genes exist", {
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")

  obj <- make_annotatable_seurat()

  suppressPackageStartupMessages(library(SingleCellExperiment))
  fake_genes <- paste0("FAKEGENE_", 1:20)
  mock_mat <- matrix(rpois(20 * 10, 5), nrow = 20, ncol = 10,
                     dimnames = list(fake_genes, paste0("R", 1:10)))
  mock_ref <- SingleCellExperiment(assays = list(logcounts = mock_mat))
  colData(mock_ref)$label.main <- rep(c("A", "B"), each = 5)

  expect_error(
    .run_singler_annotation(
      seurat_object  = obj,
      ref_data       = mock_ref,
      cluster_column = "RNA_clusters",
      assay          = "RNA",
      label_types    = "main",
      workers        = 1,
      prefix         = NULL,
      verbose        = FALSE
    ),
    regexp = "common genes|No common"
  )
})

test_that("NA in cluster column throws an error", {
  skip_if_not_installed("SingleR")
  skip_if_not_installed("SingleCellExperiment")

  obj <- make_annotatable_seurat()
  obj@meta.data$RNA_clusters[1:5] <- NA

  suppressPackageStartupMessages(library(SingleCellExperiment))
  genes <- rownames(obj)[1:20]
  mock_mat <- matrix(rpois(length(genes) * 10, 5), nrow = length(genes), ncol = 10,
                     dimnames = list(genes, paste0("R", 1:10)))
  mock_ref <- SingleCellExperiment(assays = list(logcounts = mock_mat))
  colData(mock_ref)$label.main <- rep(c("A", "B"), each = 5)

  expect_error(
    .run_singler_annotation(
      seurat_object  = obj,
      ref_data       = mock_ref,
      cluster_column = "RNA_clusters",
      assay          = "RNA",
      label_types    = "main",
      workers        = 1,
      prefix         = NULL,
      verbose        = FALSE
    ),
    regexp = "NA"
  )
})

test_that("cluster_column defaults to RNA_clusters", {
  skip_if_not_installed("Seurat")
  fmls <- formals(cell_type_annotation)
  expect_equal(as.character(fmls$cluster_column), "RNA_clusters")

  obj <- make_annotatable_seurat()
  expect_true("RNA_clusters" %in% colnames(obj@meta.data))
  expect_false("seurat_clusters" %in% colnames(obj@meta.data))
})
