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

.run_dge_main <- function(cfg, in_dir, out_dir) {
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

.dge_base_cfg <- function(...) {
  modifyList(list(input_rds     = "data/input.rds",
                  cluster_col   = "CellType",
                  sample_col    = "SampleID",
                  condition_col = "Condition",
                  timepoint_col = "Timepoint",
                  verbose       = FALSE),
             list(...))
}

test_that("a valid limma run succeeds end-to-end and leaves no error.json (regression)", {
  skip_ej_subprocess()
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  in_dir <- tempfile("dgej_ok_in_"); out_dir <- tempfile("dgej_ok_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_mini_dge_seurat(), file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(test_use = "limma",
                                   pseudobulking = TRUE,
                                   design_strategy = "stratified",
                                   timepoints_to_analyze = "T1",
                                   ident.1 = "Treated", ident.2 = "Control"),
                     in_dir, out_dir)

  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_false(file.exists(r$error_json))
  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^dge_results_.*\\.csv$"), 1)

  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^dge_passthrough_.*\\.rds$"), 1)
})

test_that("retired output keys fail at VALIDATION with CONFIG_VALIDATION_FAILED", {
  skip_ej_subprocess()
  in_dir <- tempfile("dgej_ret_in_"); out_dir <- tempfile("dgej_ret_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(output_rds = "old_style.rds"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")

  expect_identical(ej$step, "validation")
  expect_true(any(grepl("RETIRED", ej$details)))
})

test_that("pseudobulking without sample_col fails at validation", {
  skip_ej_subprocess()
  in_dir <- tempfile("dgej_scol_in_"); out_dir <- tempfile("dgej_scol_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "input.rds"))

  cfg <- .dge_base_cfg(test_use = "DESeq2", pseudobulking = TRUE)
  cfg$sample_col <- NULL
  r <- .run_dge_main(cfg, in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$step, "validation")
  expect_true(any(grepl("sample_col", ej$details)))
})

test_that("ident_2 without ident_1 fails at validation", {
  skip_ej_subprocess()
  in_dir <- tempfile("dgej_id2_in_"); out_dir <- tempfile("dgej_id2_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(ident_2 = "Control"), in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_identical(ej$step, "validation")
  expect_true(any(grepl("ident_2", ej$details)))
})

test_that("design_strategy values the core rejects now fail validation", {
  skip_ej_subprocess()
  in_dir <- tempfile("dgej_ds_in_"); out_dir <- tempfile("dgej_ds_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(1L, file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(design_strategy = "paired"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "CONFIG_VALIDATION_FAILED")
  expect_true(any(grepl("stratified, interaction", ej$details)))
})

test_that("an invalid timepoint produces INVALID_TIMEPOINTS from the analysis", {
  skip_ej_subprocess()
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  in_dir <- tempfile("dgej_tp_in_"); out_dir <- tempfile("dgej_tp_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_mini_dge_seurat(), file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(test_use = "limma",
                                   pseudobulking = TRUE,
                                   timepoints_to_analyze = "T99"),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "INVALID_TIMEPOINTS")
  expect_identical(ej$step, "dge-analysis")
  expect_true(any(grepl("T99", ej$details)))
})

test_that("an all-skipped run produces NO_DGE_RESULTS instead of an empty CSV", {
  skip_ej_subprocess()
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  in_dir <- tempfile("dgej_empty_in_"); out_dir <- tempfile("dgej_empty_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  saveRDS(make_mini_dge_seurat(), file.path(in_dir, "data", "input.rds"))

  r <- .run_dge_main(.dge_base_cfg(test_use = "limma",
                                   pseudobulking = TRUE,
                                   timepoints_to_analyze = "T1",
                                   min_samples_per_group = 999),
                     in_dir, out_dir)

  expect_identical(r$status, 1L)
  ej <- .read_error_json(r)
  expect_identical(ej$error_code, "NO_DGE_RESULTS")
  expect_true(any(grepl("min_samples_per_group: 999", ej$details)))

  expect_length(list.files(out_dir, recursive = TRUE, pattern = "\\.csv$"), 0)
})
