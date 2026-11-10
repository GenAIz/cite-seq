skip_ann <- function() {
  skip_if_not_installed("Seurat")
  if (!exists(".load_reference_dataset", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
}
skip_meta <- function() {
  skip_if_not_installed("Seurat")
  if (!exists("restore_original_metadata", mode = "function"))
    skip("add_metadata.R not sourced (Seurat unavailable at runner load)")
}

test_that(".load_reference_dataset rejects an unknown reference by name", {
  skip_ann()
  expect_error(.load_reference_dataset("not_a_real_reference"),
               "Unknown reference")
})

test_that(".load_reference_dataset error names the valid options", {
  skip_ann()
  err <- tryCatch(.load_reference_dataset("nope"), error = function(e) conditionMessage(e))
  for (ref in c("monaco", "hpca", "blueprint_encode", "dice", "immgen"))
    expect_true(grepl(ref, err, fixed = TRUE), info = ref)
})

test_that(".load_reference_dataset passes a SingleCellExperiment straight through", {
  skip_ann()
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = matrix(1:6, nrow = 3,
                         dimnames = list(paste0("g", 1:3), paste0("c", 1:2)))))
  expect_identical(.load_reference_dataset(sce), sce)
})

test_that("reference lookup is case-insensitive", {
  skip_ann()

  res <- tryCatch(.load_reference_dataset("MONACO"), error = function(e) conditionMessage(e))
  if (is.character(res)) {
    expect_false(grepl("Unknown reference", res))
  } else {
    succeed("reference resolved and returned an object")
  }
})

test_that("restore_original_metadata errors when no backup was taken", {
  skip_meta()
  obj <- make_cite_seurat()
  expect_error(restore_original_metadata(obj), "No backup metadata found")
})

test_that("restore_original_metadata puts the original columns back", {
  skip_meta()
  obj <- make_cite_seurat()
  original <- obj@meta.data
  obj@misc$original_metadata <- original
  obj$extra_column_added_later <- "x"
  expect_true("extra_column_added_later" %in% colnames(obj@meta.data))

  restored <- restore_original_metadata(obj)
  expect_false("extra_column_added_later" %in% colnames(restored@meta.data))
  expect_equal(colnames(restored@meta.data), colnames(original))
  expect_equal(nrow(restored@meta.data), ncol(obj))
})

test_that("view_metadata_join_stats is a no-op without stats", {
  skip_meta()
  obj <- make_cite_seurat()
  expect_message(res <- view_metadata_join_stats(obj), "No metadata join statistics")
  expect_null(res)
})

test_that("view_metadata_join_stats prints a recorded join", {
  skip_meta()
  obj <- make_cite_seurat()
  obj@misc$metadata_join_stats <- list(
    metadata_file = "/tmp/meta.xlsx", object_key = "orig.ident",
    metadata_key = "SampleID", timestamp = Sys.time(),
    n_cells_total = ncol(obj), n_cells_matched = ncol(obj),
    match_rate = 1, n_columns_added = 3)
  expect_output(view_metadata_join_stats(obj), "Metadata Join Statistics")
  expect_output(view_metadata_join_stats(obj), "meta.xlsx")
})

skip_add <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("readr")
  if (!exists("add_metadata", mode = "function"))
    skip("add_metadata.R not sourced")
}

.meta_fixture <- function(rows) {
  obj <- make_cite_seurat()
  f <- tempfile(fileext = ".csv")
  write.csv(rows, f, row.names = FALSE)
  list(obj = obj, file = f)
}

test_that("duplicate metadata keys do not multiply cells", {
  skip_add()
  fx <- .meta_fixture(data.frame(
    Participant = c("P1", "P1", "P2"),
    Cohort      = c("A", "SHOULD_BE_DROPPED", "B"),
    stringsAsFactors = FALSE))
  n_cells <- ncol(fx$obj)
  out <- suppressWarnings(add_metadata(
    seurat_object = fx$obj, metadata_file = fx$file,
    object_key = "Participant", metadata_key = "Participant", verbose = FALSE))

  expect_equal(ncol(out), n_cells)
  expect_equal(nrow(out@meta.data), n_cells)
  expect_true(all(out@meta.data$Cohort[out@meta.data$Participant == "P1"] == "A"))
})

test_that("overwrite_existing = FALSE preserves the object's column", {
  skip_add()
  fx <- .meta_fixture(data.frame(
    Participant = c("P1", "P2"),
    Timepoint   = c("FROM_FILE", "FROM_FILE"),
    stringsAsFactors = FALSE))
  out <- suppressWarnings(add_metadata(
    seurat_object = fx$obj, metadata_file = fx$file,
    object_key = "Participant", metadata_key = "Participant",
    overwrite_existing = FALSE, verbose = FALSE))

  expect_true(all(out@meta.data$Timepoint %in% c("T0", "T1")))
  expect_false(any(grepl("\\.(x|y)$", colnames(out@meta.data))))
})

test_that("overwrite_existing = TRUE replaces the object's column", {
  skip_add()
  fx <- .meta_fixture(data.frame(
    Participant = c("P1", "P2"),
    Timepoint   = c("FROM_FILE", "FROM_FILE"),
    stringsAsFactors = FALSE))
  out <- suppressWarnings(add_metadata(
    seurat_object = fx$obj, metadata_file = fx$file,
    object_key = "Participant", metadata_key = "Participant",
    overwrite_existing = TRUE, verbose = FALSE))
  expect_true(all(out@meta.data$Timepoint == "FROM_FILE"))
  expect_false(any(grepl("\\.(x|y)$", colnames(out@meta.data))))
})
