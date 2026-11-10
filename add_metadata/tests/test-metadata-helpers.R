test_that("view_metadata_join_stats shows stats after a join", {
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

  expect_output(
    stats <- view_metadata_join_stats(result),
    regexp = "Metadata Join Statistics"
  )
  expect_false(is.null(stats))
  expect_equal(stats$n_cells_total, ncol(obj))
})

test_that("view_metadata_join_stats returns NULL when no stats exist", {
  obj <- make_mini_seurat()
  expect_message(
    result <- view_metadata_join_stats(obj),
    regexp = "No metadata join statistics"
  )
  expect_null(result)
})

test_that("restore_original_metadata reverts to backup", {
  obj <- make_mini_seurat()
  tmp <- tempfile(fileext = ".csv")
  write_meta_csv(tmp)
  on.exit(unlink(tmp))

  original_cols <- ncol(obj@meta.data)

  result <- add_metadata(
    seurat_object   = obj,
    metadata_file   = tmp,
    object_key      = "orig.ident",
    metadata_key    = "Sample",
    backup_metadata = TRUE,
    verbose         = FALSE
  )

  expect_gt(ncol(result@meta.data), original_cols)

  restored <- restore_original_metadata(result)
  expect_equal(ncol(restored@meta.data), original_cols)
})

test_that("restore_original_metadata errors without backup", {
  obj <- make_mini_seurat()
  expect_error(
    restore_original_metadata(obj),
    regexp = "No backup|backup_metadata"
  )
})
