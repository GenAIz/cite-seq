test_that("basic join adds new columns to metadata", {
  obj  <- make_mini_seurat()
  tmp  <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object = obj,
    metadata_file = tmp,
    object_key    = "orig.ident",
    metadata_key  = "Sample",
    verbose       = FALSE
  )

  expect_true("Treatment" %in% colnames(result@meta.data))
  expect_true("Score"     %in% colnames(result@meta.data))
})

test_that("cell barcodes are preserved after join", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object = obj,
    metadata_file = tmp,
    object_key    = "orig.ident",
    metadata_key  = "Sample",
    verbose       = FALSE
  )

  expect_identical(colnames(result), colnames(obj))
  expect_equal(nrow(result@meta.data), ncol(obj))
})

test_that("cells with no metadata match receive NA in joined columns", {
  obj <- make_mini_seurat()

  tmp <- tempfile(fileext = ".csv")
  df  <- data.frame(Sample = paste0("Sample", 1:3), Extra = "X",
                    stringsAsFactors = FALSE)
  write.csv(df, tmp, row.names = FALSE)
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object = obj,
    metadata_file = tmp,
    object_key    = "orig.ident",
    metadata_key  = "Sample",
    verbose       = FALSE
  )

  unmatched <- result@meta.data$Extra[result@meta.data$orig.ident %in% paste0("Sample", 4:5)]
  expect_true(all(is.na(unmatched)))
})

test_that("metadata_matched column is correct when add_match_status = TRUE", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp, samples = paste0("Sample", 1:3))
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object  = obj,
    metadata_file  = tmp,
    object_key     = "orig.ident",
    metadata_key   = "Sample",
    add_match_status = TRUE,
    verbose        = FALSE
  )

  expect_true("metadata_matched" %in% colnames(result@meta.data))
  n_matched   <- sum(result$metadata_matched, na.rm = TRUE)
  n_unmatched <- sum(!result$metadata_matched, na.rm = TRUE)
  expect_gt(n_matched,   0L)
  expect_gt(n_unmatched, 0L)
})

test_that("missing metadata file throws an informative error", {
  obj <- make_mini_seurat()
  expect_error(
    add_metadata(
      seurat_object = obj,
      metadata_file = "/nonexistent/path/meta.csv",
      object_key    = "orig.ident",
      metadata_key  = "Sample",
      verbose       = FALSE
    ),
    regexp = "not found|does not exist"
  )
})

test_that("missing object_key throws an informative error", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  expect_error(
    add_metadata(
      seurat_object = obj,
      metadata_file = tmp,
      object_key    = "nonexistent_column",
      metadata_key  = "Sample",
      verbose       = FALSE
    ),
    regexp = "not found|nonexistent_column"
  )
})

test_that("join statistics are stored in misc$metadata_join_stats", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object = obj,
    metadata_file = tmp,
    object_key    = "orig.ident",
    metadata_key  = "Sample",
    verbose       = FALSE
  )

  stats <- result@misc$metadata_join_stats
  expect_false(is.null(stats))
  expect_true(all(c("n_cells_total", "match_rate", "columns_added") %in% names(stats)))
  expect_equal(stats$n_cells_total, ncol(obj))
})

test_that("backup_metadata = TRUE stores original in misc$original_metadata", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  result <- add_metadata(
    seurat_object   = obj,
    metadata_file   = tmp,
    object_key      = "orig.ident",
    metadata_key    = "Sample",
    backup_metadata = TRUE,
    verbose         = FALSE
  )

  expect_false(is.null(result@misc$original_metadata))
  expect_equal(nrow(result@misc$original_metadata), ncol(obj))
})

test_that("Excel file requires readxl and fails with useful message if absent", {

  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".xlsx")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  assign("requireNamespace", function(package, ...) {
    if (identical(package, "readxl")) return(FALSE)
    base::requireNamespace(package, ...)
  }, envir = globalenv())
  on.exit(rm(list = "requireNamespace", envir = globalenv()), add = TRUE)

  expect_error(
    add_metadata(
      seurat_object = obj,
      metadata_file = tmp,
      object_key    = "orig.ident",
      metadata_key  = "Sample",
      verbose       = FALSE
    ),
    regexp = "readxl.+required for Excel"
  )
})
