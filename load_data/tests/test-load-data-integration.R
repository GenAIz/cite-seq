test_that("load_data from RDS returns object as-is", {
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))

  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds))

  counts <- matrix(sample(0:10, 100, replace = TRUE), nrow = 10,
                   dimnames = list(paste0("G", 1:10), paste0("C", 1:10)))
  tiny   <- CreateSeuratObject(counts = as(counts, "sparseMatrix"))
  saveRDS(tiny, tmp_rds)

  obj <- load_data(rna_path = tmp_rds, verbose = FALSE)
  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 10L)
})

test_that("load_data with non-existent file throws an error", {
  expect_error(
    load_data(rna_path = "/no/such/file.rds", verbose = FALSE),
    regexp = "not found|does not exist|cannot open|no such file"
  )
})

test_that("load_data rejects raw formats with an explanatory error", {
  for (ext in c(".csv", ".csv.gz", ".tar.gz", ".h5")) {
    tmp <- tempfile(fileext = ext)
    writeLines("dummy", tmp)
    expect_error(load_data(rna_path = tmp, verbose = FALSE),
                 regexp = "PRE-PROCESSED|Unsupported input format")
    unlink(tmp)
  }
})

test_that("load_data with unsupported extension throws an error", {
  tmp <- tempfile(fileext = ".xyz")
  writeLines("dummy", tmp)
  on.exit(unlink(tmp))

  expect_error(
    load_data(rna_path = tmp, verbose = FALSE),
    regexp = "Unsupported|format"
  )
})

test_that("load_data rejects an .rds that is not a Seurat object", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(as(matrix(1:20, nrow = 4), "sparseMatrix"), tmp)

  expect_error(load_data(rna_path = tmp, verbose = FALSE),
               regexp = "not a Seurat object|pre-processed Seurat")
})

test_that("adt_path is rejected", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(as(matrix(1:20, nrow = 4), "sparseMatrix"), tmp)

  expect_error(load_data(rna_path = tmp, adt_path = tmp, verbose = FALSE),
               regexp = "adt_path is retired")
})

test_that("a directory of .rds files is no longer accepted", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  saveRDS(1:5, file.path(d, "a.rds"))

  expect_error(load_data(rna_path = d, verbose = FALSE),
               regexp = "Unsupported input format|PRE-PROCESSED")
})
