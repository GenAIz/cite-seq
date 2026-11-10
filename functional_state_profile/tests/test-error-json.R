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
  expect_identical(parsed$message, "Something went wrong.")
  expect_identical(parsed$details, c("detail one", "detail two"))
  expect_identical(parsed$step, "unit-test")
  expect_match(parsed$timestamp, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}")
})

test_that("sf_write_error_json creates the output directory if absent", {
  skip_if_not_installed("jsonlite")
  if (!exists("sf_write_error_json", mode = "function")) {
    source(file.path(.ej_sf_root, "error_codes.R"))
  }
  out_dir <- file.path(tempfile("ej_deep_"), "nested", "out")
  old_env <- Sys.getenv("SF_OUTPUT_PATH", unset = NA)
  Sys.setenv(SF_OUTPUT_PATH = out_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("SF_OUTPUT_PATH")
    else Sys.setenv(SF_OUTPUT_PATH = old_env)
    unlink(dirname(dirname(out_dir)), recursive = TRUE)
  })

  path <- sf_write_error_json("X", "y")
  expect_true(file.exists(path))
})

.run_fsp_main <- function(cfg, in_dir, out_dir) {
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
}

test_that("an invalid config produces error.json with CONFIG_VALIDATION_FAILED", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ej_badcfg_in_"); out_dir <- tempfile("ej_badcfg_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "annotated_seurat.rds"))
  r <- .run_fsp_main(list(input_rds = "data/annotated_seurat.rds",
                          min_score = 2), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$smart_function, "functional_state_profile")
  expect_true(any(grepl("min_score", ej$details)))
  expect_identical(ej$step, "validation")
})

test_that("a missing input file produces error.json at validation", {
  skip_ej_subprocess()
  in_dir  <- tempfile("ej_noin_in_"); out_dir <- tempfile("ej_noin_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  r <- .run_fsp_main(list(input_rds = "data/does_not_exist.rds"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_true(any(grepl("does_not_exist", ej$details)))
})

test_that("a non-Seurat input produces error.json with INPUT_NOT_SEURAT", {
  skip_ej_subprocess()
  skip_if_not_installed("UCell")
  in_dir  <- tempfile("ej_notseu_in_"); out_dir <- tempfile("ej_notseu_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(list(not = "a seurat object"),
          file.path(in_dir, "data", "annotated_seurat.rds"))
  file.copy(TEST_MARKERS_JSON, file.path(in_dir, "test_markers.json"))
  r <- .run_fsp_main(list(input_rds = "data/annotated_seurat.rds",
                          marker_json_path = "test_markers.json"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "INPUT_NOT_SEURAT")
  expect_identical(ej$step, "classification")
})

test_that("markers absent from the data produce error.json with MARKER_STATE_EMPTY", {
  skip_ej_subprocess()
  skip_if_not_installed("UCell")
  in_dir  <- tempfile("ej_empty_in_"); out_dir <- tempfile("ej_empty_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))

  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_fsp_seurat(), file.path(in_dir, "data", "annotated_seurat.rds"))

  writeLines(jsonlite::toJSON(
    list(ghost_state = list(core = c("NOT_A_GENE_1", "NOT_A_GENE_2"))),
    auto_unbox = TRUE), file.path(in_dir, "ghost_markers.json"))
  r <- .run_fsp_main(list(input_rds = "data/annotated_seurat.rds",
                          marker_json_path = "ghost_markers.json"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "MARKER_STATE_EMPTY")
  expect_true(any(grepl("ghost_state", ej$details)))
})
