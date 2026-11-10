.hg_main <- local({
  cands <- c("main.R", file.path("..", "main.R"))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) NA_character_ else normalizePath(hit)
})

.hg_helper <- local({
  if (is.na(.hg_main)) return(NULL)
  exprs <- tryCatch(parse(.hg_main), error = function(e) NULL)
  if (is.null(exprs)) return(NULL)
  env <- new.env()
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(as.character(e[[2]]), ".h5_backend_needed")) {
      eval(e, env)
      return(get(".h5_backend_needed", envir = env))
    }
  }
  NULL
})

skip_hg <- function() {
  if (is.na(.hg_main)) skip("main.R not found from test wd")
  if (is.null(.hg_helper)) skip(".h5_backend_needed not found in main.R")
}

test_that("guard triggers for .h5ad and .h5mu files, by extension, any case", {
  skip_hg()
  expect_true(.hg_helper("data/input.h5ad"))
  expect_true(.hg_helper("data/input.H5AD"))
  expect_true(.hg_helper("data/multimodal.h5mu"))
  expect_true(.hg_helper("/mnt/in/data/Yoshida_et_al_2022_airway 1.h5ad"))
})

test_that("guard triggers for a directory containing h5ad/h5mu files", {
  skip_hg()
  d <- tempfile("hg_dir_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  expect_false(.hg_helper(d))
  file.create(file.path(d, "sample.h5ad"))
  expect_true(.hg_helper(d))
})

test_that("guard does NOT trigger for rds/h5 inputs (BPCells path, no python)", {
  skip_hg()
  expect_false(.hg_helper("data/object.rds"))
  expect_false(.hg_helper("data/matrix.h5"))
  expect_false(.hg_helper("data/archive.tar.gz"))
  d <- tempfile("hg_rds_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  file.create(file.path(d, "object.rds"))
  expect_false(.hg_helper(d))
})

test_that("main.R wires the guard to sf_fail with the stable error code", {
  skip_hg()
  src <- paste(readLines(.hg_main), collapse = "\n")
  expect_true(grepl("H5AD_PYTHON_BACKEND_UNAVAILABLE", src, fixed = TRUE))

  expect_true(grepl("py_run_string", src, fixed = TRUE))
  expect_true(grepl("import h5py, anndata", src, fixed = TRUE))
})
