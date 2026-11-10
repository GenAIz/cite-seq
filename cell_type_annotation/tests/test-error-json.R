.ej_sf_root <- local({

  cands <- c("main.R", file.path("..", "main.R"))
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
                              details = "one detail")

  expect_true(file.exists(path))

  expect_true(file.exists(file.path(out_dir, "output", "error.json")))
  parsed <- jsonlite::fromJSON(path)
  expect_identical(parsed$status, "error")
  expect_identical(parsed$error_code, "SOME_CODE")
  expect_identical(parsed$step, "unit-test")
})

.run_cta_main <- function(cfg, in_dir, out_dir) {
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
  skip_if_not_installed("SingleR")
  skip_if_not_installed("celldex")
}

test_that("an invalid reference produces CONFIG_VALIDATION_FAILED listing the real names", {
  skip_ej_subprocess()
  in_dir <- tempfile("ctaej_ref_in_"); out_dir <- tempfile("ctaej_ref_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "input.rds"))

  r <- .run_cta_main(list(input_rds = "data/input.rds", reference = "zeisel"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$smart_function, "cell_type_annotation")
  expect_true(any(grepl("blueprint_encode", ej$details)))
})

test_that("the loader's real reference names pass validation (regression)", {
  skip_ej_subprocess()
  in_dir <- tempfile("ctaej_be_in_"); out_dir <- tempfile("ctaej_be_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)

  saveRDS(list(x = 1), file.path(in_dir, "data", "input.rds"))

  r <- .run_cta_main(list(input_rds = "data/input.rds",
                          reference = "blueprint_encode",
                          method = "singler"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)

  expect_identical(ej$error_code, "INPUT_NOT_SEURAT")
  expect_identical(ej$step, "annotation")
})

test_that("a missing cluster column produces CLUSTER_COLUMN_INVALID", {
  skip_ej_subprocess()
  in_dir <- tempfile("ctaej_cl_in_"); out_dir <- tempfile("ctaej_cl_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  obj <- make_annotatable_seurat()
  obj@meta.data$RNA_clusters <- NULL
  saveRDS(obj, file.path(in_dir, "data", "input.rds"))

  r <- .run_cta_main(list(input_rds = "data/input.rds",
                          method = "SingleR"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CLUSTER_COLUMN_INVALID")

  expect_true(any(grepl("Available columns", ej$details)))
})

test_that("a missing input file produces CONFIG_VALIDATION_FAILED at validation", {
  skip_ej_subprocess()
  in_dir <- tempfile("ctaej_noin_in_"); out_dir <- tempfile("ctaej_noin_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  r <- .run_cta_main(list(input_rds = "data/absent.rds"), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_true(any(grepl("absent.rds", ej$details)))
})
