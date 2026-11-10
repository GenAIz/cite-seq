test_that("valid config passes schema check", {
  cfg <- make_valid_config()
  result <- validate_config_schema(cfg)
  expect_true(result$valid)
})

test_that("missing top-level fields fail schema check", {
  cfg <- list(metadata = list(), paths = list())
  result <- validate_config_schema(cfg)
  expect_false(result$valid)
  expect_match(result$message, "Missing")
})

test_that("each missing field is reported", {
  cfg <- list(metadata = list())
  result <- validate_config_schema(cfg)
  expect_false(result$valid)

  expect_match(result$message, "paths")
  expect_match(result$message, "execution")
})

test_that("existing RNA file passes input validation", {
  tmp <- tempfile(fileext = ".rds")
  writeLines("x", tmp)
  on.exit(unlink(tmp))

  cfg <- make_valid_config()
  cfg$paths$input_data$rna_path <- tmp

  result <- validate_input_files(cfg)
  expect_true(result$valid)
})

test_that("missing RNA file fails input validation", {
  cfg <- make_valid_config()
  cfg$paths$input_data$rna_path <- "/nonexistent/file.rds"

  result <- validate_input_files(cfg)
  expect_false(result$valid)
  expect_match(result$message, "MISSING|not found")
})

test_that("metadata file is checked when add_metadata is enabled", {
  tmp_rna <- tempfile(fileext = ".rds")
  writeLines("x", tmp_rna)
  on.exit(unlink(tmp_rna))

  cfg <- make_valid_config()
  cfg$paths$input_data$rna_path <- tmp_rna
  cfg$add_metadata$enabled <- TRUE
  cfg$paths$input_data$metadata_file <- "/nonexistent/meta.xlsx"

  result <- validate_input_files(cfg)
  expect_false(result$valid)
  expect_match(result$message, "Metadata|MISSING")
})

test_that("metadata file is NOT checked when add_metadata is disabled", {
  tmp_rna <- tempfile(fileext = ".rds")
  writeLines("x", tmp_rna)
  on.exit(unlink(tmp_rna))

  cfg <- make_valid_config()
  cfg$paths$input_data$rna_path <- tmp_rna
  cfg$add_metadata$enabled <- FALSE
  cfg$paths$input_data$metadata_file <- "/nonexistent/meta.xlsx"

  result <- validate_input_files(cfg)
  expect_true(result$valid)
})

test_that("valid parameter ranges pass", {
  cfg <- make_valid_config()
  result <- validate_parameter_ranges(cfg)
  expect_true(result$valid)
})

test_that("negative min_features fails", {
  cfg <- make_valid_config()
  cfg$processing$qc_filters$min_features <- -1
  result <- validate_parameter_ranges(cfg)
  expect_false(result$valid)
})

test_that("max_features < min_features fails", {
  cfg <- make_valid_config()
  cfg$processing$qc_filters$min_features <- 500
  cfg$processing$qc_filters$max_features <- 100
  result <- validate_parameter_ranges(cfg)
  expect_false(result$valid)
})

test_that("max_percent_mt > 100 fails", {
  cfg <- make_valid_config()
  cfg$processing$qc_filters$max_percent_mt <- 150
  result <- validate_parameter_ranges(cfg)
  expect_false(result$valid)
})

test_that("clustering resolution <= 0 fails", {
  cfg <- make_valid_config()
  cfg$processing$clustering$rna_resolution <- 0
  result <- validate_parameter_ranges(cfg)
  expect_false(result$valid)
})

test_that("n_cores < 1 fails", {
  cfg <- make_valid_config()
  cfg$execution$n_cores <- 0
  result <- validate_parameter_ranges(cfg)
  expect_false(result$valid)
})

test_that("writable directories pass", {
  tmp <- tempfile("orch_dir_test_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE))

  cfg <- make_valid_config(tmp)
  cfg$paths$output_dir <- file.path(tmp, "out")

  old_env <- Sys.getenv("SF_LOG_DIR", unset = NA)
  Sys.setenv(SF_LOG_DIR = file.path(tmp, "logs"))
  on.exit(if (is.na(old_env)) Sys.unsetenv("SF_LOG_DIR") else Sys.setenv(SF_LOG_DIR = old_env), add = TRUE)

  result <- ensure_output_directories(cfg)
  expect_true(result$valid)
})

test_that("full valid config passes validate_workflow_config", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  cfg <- jsonlite::fromJSON(config_path)

  old_env <- Sys.getenv("SF_LOG_DIR", unset = NA)
  Sys.setenv(SF_LOG_DIR = cfg$paths$logs_dir)
  on.exit(if (is.na(old_env)) Sys.unsetenv("SF_LOG_DIR") else Sys.setenv(SF_LOG_DIR = old_env), add = TRUE)

  result <- validate_workflow_config(cfg)
  expect_true(result$valid)
})

test_that("config missing required fields fails validate_workflow_config", {
  cfg <- list(metadata = list())

  old_env <- Sys.getenv("SF_LOG_DIR", unset = NA)
  Sys.setenv(SF_LOG_DIR = tempdir())
  on.exit(if (is.na(old_env)) Sys.unsetenv("SF_LOG_DIR") else Sys.setenv(SF_LOG_DIR = old_env), add = TRUE)

  result <- validate_workflow_config(cfg)
  expect_false(result$valid)
})
