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
                              details = "one detail")

  expect_true(file.exists(path))

  expect_true(file.exists(file.path(out_dir, "output", "error.json")))
  parsed <- jsonlite::fromJSON(path)
  expect_identical(parsed$status, "error")
  expect_identical(parsed$error_code, "SOME_CODE")
  expect_identical(parsed$step, "unit-test")
})

.run_am_main <- function(cfg, in_dir, out_dir) {
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

.am_write_inputs <- function(in_dir, samples = paste0("Sample", 1:5)) {
  dir.create(file.path(in_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(make_mini_seurat(), file.path(in_dir, "data", "input.rds"))
  write_meta_csv(file.path(in_dir, "data", "meta.csv"), samples = samples)
}

test_that("a valid run succeeds, writes output, and leaves no error.json", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_ok_in_"); out_dir <- tempfile("amej_ok_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(in_dir, recursive = TRUE); .am_write_inputs(in_dir)

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.csv",
                         object_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_false(file.exists(r$error_json))
  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^metadata_.*\\.rds$"), 1)
})

test_that("the legacy seurat_key spelling is still accepted", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_legacy_in_"); out_dir <- tempfile("amej_legacy_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(in_dir, recursive = TRUE); .am_write_inputs(in_dir)

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.csv",
                         seurat_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_false(file.exists(r$error_json))
})

test_that("a metadata file on the attachment port is discovered by name", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_att_in_"); out_dir <- tempfile("amej_att_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  dir.create(file.path(in_dir, "attachment"), recursive = TRUE)
  saveRDS(make_mini_seurat(), file.path(in_dir, "data", "input.rds"))
  write_meta_csv(file.path(in_dir, "attachment", "meta.csv"))

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "meta.csv",
                         object_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_false(file.exists(r$error_json))
})

test_that("unset metadata_file falls back to the newest attachment-port file", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_att2_in_"); out_dir <- tempfile("amej_att2_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  dir.create(file.path(in_dir, "attachment"), recursive = TRUE)
  saveRDS(make_mini_seurat(), file.path(in_dir, "data", "input.rds"))
  write_meta_csv(file.path(in_dir, "attachment", "uploaded.csv"))

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         object_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_false(file.exists(r$error_json))
})

test_that("missing required config keys produce CONFIG_VALIDATION_FAILED", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_cfg_in_"); out_dir <- tempfile("amej_cfg_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  r <- .run_am_main(list(input_rds = "data/input.rds"), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$smart_function, "add_metadata")
  expect_true(any(grepl("metadata_file", ej$details)))
  expect_true(any(grepl("object_key", ej$details)))
})

test_that("a wrong metadata_key produces METADATA_KEY_NOT_FOUND", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_key_in_"); out_dir <- tempfile("amej_key_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(in_dir, recursive = TRUE); .am_write_inputs(in_dir)

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.csv",
                         object_key = "orig.ident",
                         metadata_key = "NotAColumn",
                         verbose = FALSE), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "METADATA_KEY_NOT_FOUND")
  expect_identical(ej$step, "metadata-join")

  expect_true(any(grepl("Sample", ej$details)))
})

test_that("a wrong object_key produces OBJECT_KEY_NOT_FOUND", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_skey_in_"); out_dir <- tempfile("amej_skey_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(in_dir, recursive = TRUE); .am_write_inputs(in_dir)

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.csv",
                         object_key = "NotAColumn",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "OBJECT_KEY_NOT_FOUND")
})

test_that("an unparseable metadata file produces METADATA_LOAD_FAILED", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_bad_in_"); out_dir <- tempfile("amej_bad_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_mini_seurat(), file.path(in_dir, "data", "input.rds"))

  writeLines("this is not a workbook", file.path(in_dir, "data", "meta.xlsx"))

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.xlsx",
                         object_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "METADATA_LOAD_FAILED")
})

test_that("an .rds metadata file produces CONFIG_VALIDATION_FAILED", {
  skip_ej_subprocess()
  in_dir <- tempfile("amej_rdsm_in_"); out_dir <- tempfile("amej_rdsm_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_mini_seurat(), file.path(in_dir, "data", "input.rds"))
  saveRDS(data.frame(Sample = paste0("Sample", 1:5)),
          file.path(in_dir, "data", "meta.rds"))

  r <- .run_am_main(list(input_rds = "data/input.rds",
                         metadata_file = "data/meta.rds",
                         object_key = "orig.ident",
                         metadata_key = "Sample",
                         verbose = FALSE), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_true(any(grepl("unsupported format", ej$details)))
})
