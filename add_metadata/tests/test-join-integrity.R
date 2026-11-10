skip_join <- function() {
  skip_if_not_installed("Seurat")
  if (!exists("add_metadata", mode = "function"))
    skip("add_metadata.R not sourced")
}

write_df <- function(df) {
  f <- tempfile(fileext = ".csv")
  write.csv(df, f, row.names = FALSE)
  f
}

test_that("duplicate join keys do not multiply cells", {
  skip_join()
  obj <- make_mini_seurat()
  n   <- ncol(obj)
  f <- write_df(data.frame(
    Sample    = c("Sample1", "Sample1", "Sample2"),
    Treatment = c("Control", "SHOULD_BE_DROPPED", "Treated"),
    stringsAsFactors = FALSE))
  on.exit(unlink(f))

  out <- suppressWarnings(add_metadata(
    seurat_object = obj, metadata_file = f,
    object_key = "orig.ident", metadata_key = "Sample", verbose = FALSE))

  expect_equal(ncol(out), n)
  expect_equal(nrow(out@meta.data), n)
  expect_true(all(out@meta.data$Treatment[out@meta.data$orig.ident == "Sample1"] == "Control"))
})

test_that("overwrite_existing = FALSE preserves the object's column", {
  skip_join()
  obj <- make_mini_seurat()
  obj$Treatment <- "FROM_OBJECT"
  f <- write_df(data.frame(
    Sample    = paste0("Sample", 1:5),
    Treatment = "FROM_FILE",
    stringsAsFactors = FALSE))
  on.exit(unlink(f))

  out <- suppressWarnings(add_metadata(
    seurat_object = obj, metadata_file = f,
    object_key = "orig.ident", metadata_key = "Sample",
    overwrite_existing = FALSE, verbose = FALSE))

  expect_true(all(out@meta.data$Treatment == "FROM_OBJECT"))
  expect_false(any(grepl("\\.(x|y)$", colnames(out@meta.data))))
})

test_that("overwrite_existing = TRUE replaces the object's column", {
  skip_join()
  obj <- make_mini_seurat()
  obj$Treatment <- "FROM_OBJECT"
  f <- write_df(data.frame(
    Sample    = paste0("Sample", 1:5),
    Treatment = "FROM_FILE",
    stringsAsFactors = FALSE))
  on.exit(unlink(f))

  out <- suppressWarnings(add_metadata(
    seurat_object = obj, metadata_file = f,
    object_key = "orig.ident", metadata_key = "Sample",
    overwrite_existing = TRUE, verbose = FALSE))

  expect_true(all(out@meta.data$Treatment == "FROM_FILE"))
  expect_false(any(grepl("\\.(x|y)$", colnames(out@meta.data))))
})
