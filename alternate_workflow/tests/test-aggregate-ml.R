skip_agg <- function() {
  skip_if_not_installed("Seurat")
  if (!exists("aggregate_seurat_for_ml_final", mode = "function"))
    skip("aggregate_obj_ML_processing.R not sourced (Seurat unavailable at runner load)")
}

make_agg_obj <- function() {
  obj <- make_cite_seurat()
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
    DefaultAssay(obj) <- "ADT"
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
    DefaultAssay(obj) <- "RNA"
  })
  obj
}

test_that("aggregate_seurat_for_ml_final returns results for a processed object", {
  skip_agg()
  obj <- make_agg_obj()
  res <- aggregate_seurat_for_ml_final(
    seurat_obj      = obj,
    participant_col = "Participant",
    timepoint_col   = "Timepoint",
    celltype_col    = "CellType",
    assays          = c("RNA", "ADT"),
    output_xlsx     = NULL
  )
  expect_false(is.null(res))
})

test_that("aggregate errors when a grouping column is missing", {
  skip_agg()
  obj <- make_agg_obj()
  expect_error(
    aggregate_seurat_for_ml_final(
      seurat_obj = obj, participant_col = "NoSuchCol",
      timepoint_col = "Timepoint", celltype_col = "CellType",
      assays = "RNA", output_xlsx = NULL),
    regexp = "missing"
  )
})

test_that("aggregate works when RNA counts/data are BPCells-backed", {
  skip_agg()
  skip_if_not_installed("BPCells")

  obj <- make_cite_seurat()
  obj <- .to_disk_backed(obj, cache_dir = tempfile("bpc_"), assays = "RNA", verbose = FALSE)
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
    DefaultAssay(obj) <- "ADT"
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
    DefaultAssay(obj) <- "RNA"
  })
  res <- aggregate_seurat_for_ml_final(
    seurat_obj = obj, participant_col = "Participant",
    timepoint_col = "Timepoint", celltype_col = "CellType",
    assays = c("RNA", "ADT"), output_xlsx = NULL)
  expect_false(is.null(res))
})

test_that("aggregate writes an xlsx when output_xlsx is provided", {
  skip_agg()
  skip_if_not_installed("openxlsx")
  obj <- make_agg_obj()
  out <- tempfile(fileext = ".xlsx")
  on.exit(unlink(out), add = TRUE)
  aggregate_seurat_for_ml_final(
    seurat_obj = obj, participant_col = "Participant",
    timepoint_col = "Timepoint", celltype_col = "CellType",
    assays = c("RNA", "ADT"), output_xlsx = out)
  expect_true(file.exists(out))
})

test_that("aggregate handles an assay with a single feature", {
  skip_agg()
  obj <- make_agg_obj()

  suppressWarnings({
    one <- subset(obj[["ADT"]], features = rownames(obj[["ADT"]])[1])
    obj[["ADT"]] <- one
  })
  expect_equal(nrow(obj[["ADT"]]), 1L)

  res <- aggregate_seurat_for_ml_final(
    seurat_obj      = obj,
    participant_col = "Participant",
    timepoint_col   = "Timepoint",
    celltype_col    = "CellType",
    assays          = "ADT",
    output_xlsx     = NULL
  )
  expect_false(is.null(res))

  n_groups <- nrow(unique(obj@meta.data[, c("Participant", "Timepoint", "CellType")]))
  tbl <- if (is.data.frame(res)) res else res[[1]]
  expect_equal(nrow(tbl), n_groups)
})

test_that("aggregate tolerates VariableFeatures not present in the matrix", {
  skip_agg()
  obj <- make_agg_obj()

  VariableFeatures(obj) <- c(VariableFeatures(obj)[1:5], "GeneThatDoesNotExist")
  expect_no_error(
    aggregate_seurat_for_ml_final(
      seurat_obj      = obj,
      participant_col = "Participant",
      timepoint_col   = "Timepoint",
      celltype_col    = "CellType",
      assays          = "RNA",
      output_xlsx     = NULL
    )
  )
})
