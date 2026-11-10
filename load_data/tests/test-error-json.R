.ej_sf_root <- local({
  cands <- c("main.R", file.path("..", "main.R"),
             file.path(dirname(dirname(FIXTURES)), "main.R"))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("cannot locate main.R from test wd")
  dirname(normalizePath(hit))
})

test_that("sf_write_error_json writes a complete, parseable error.json", {
  skip_if_not_installed("jsonlite")
  if (!exists("sf_write_error_json", mode = "function")) {
    source(file.path(.ej_sf_root, "error_codes.R"))
  }

  out_dir <- tempfile("ej_out_")
  old_env <- Sys.getenv("SF_OUTPUT_PATH", unset = NA)
  Sys.setenv(SF_OUTPUT_PATH = out_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("SF_OUTPUT_PATH")
    else Sys.setenv(SF_OUTPUT_PATH = old_env)
    unlink(out_dir, recursive = TRUE)
  })

  sf_error_init("test_function")
  sf_step("unit-test")
  path <- sf_write_error_json("SOME_CODE", "Something went wrong.",
                              details = c("detail one", "detail two"))

  expect_true(file.exists(path))

  expect_true(file.exists(file.path(out_dir, "output", "error.json")))
  parsed <- jsonlite::fromJSON(path)
  expect_identical(parsed$status, "error")
  expect_identical(parsed$smart_function, "test_function")
  expect_identical(parsed$error_code, "SOME_CODE")
  expect_identical(parsed$details, c("detail one", "detail two"))
  expect_identical(parsed$step, "unit-test")
})

.run_ld_main <- function(cfg, in_dir, out_dir) {
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), cfg_path)

  old <- setwd(.ej_sf_root); on.exit(setwd(old))
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", out_dir)),
    stdout = TRUE, stderr = TRUE
  ))
  list(status = attr(res, "status") %||% 0L,
       output = paste(res, collapse = "\n"),
       out_dir = out_dir,
       error_json = file.path(out_dir, "error.json"))
}

.read_error_json <- function(r) {
  expect_true(file.exists(r$error_json),
              info = paste("no error.json written; output was:\n", r$output))
  jsonlite::fromJSON(r$error_json)
}

skip_ej_subprocess <- function() {
  skip_if_not_installed("optparse")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("Seurat")
}

test_that("a config that matches config.example.json passes validation and runs (regression)", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ldej_ok_in_"); out_dir <- tempfile("ldej_ok_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  obj <- make_tiny_seurat()
  saveRDS(obj, file.path(in_dir, "data", "preprocessed_object.rds"))

  r <- .run_ld_main(list(load_data_path = "data/preprocessed_object.rds",
                         project_name   = "SmokeTest",
                         verbose        = FALSE),
                    in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)

  expect_false(file.exists(r$error_json))

  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^loaded_data_.*\\.rds$"), 1)
})

test_that("a missing load_data_path produces CONFIG_VALIDATION_FAILED naming the real key", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ldej_nokey_in_"); out_dir <- tempfile("ldej_nokey_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  r <- .run_ld_main(list(project_name = "X"), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$smart_function, "load_data")

  expect_true(any(grepl("load_data_path", ej$details)))
  expect_false(any(grepl("rna_path", ej$details)))
})

test_that("an unsupported input format produces UNSUPPORTED_FORMAT", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ldej_fmt_in_"); out_dir <- tempfile("ldej_fmt_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  writeLines("gene,cell1\nCD3E,4", file.path(in_dir, "data", "raw_counts.csv"))
  r <- .run_ld_main(list(load_data_path = "data/raw_counts.csv"),
                    in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "UNSUPPORTED_FORMAT")
  expect_identical(ej$step, "load-data")
})

test_that("an .rds that is not a Seurat object produces INPUT_NOT_SEURAT", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ldej_mat_in_"); out_dir <- tempfile("ldej_mat_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(matrix(1:4, 2), file.path(in_dir, "data", "bare_matrix.rds"))
  r <- .run_ld_main(list(load_data_path = "data/bare_matrix.rds"),
                    in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "INPUT_NOT_SEURAT")
})
