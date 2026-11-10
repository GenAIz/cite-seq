test_that("initialize_logger creates a log file", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  expect_true(.logger$initialized)
  expect_true(file.exists(.logger$log_file))

  content <- readLines(.logger$log_file)
  expect_true(any(grepl("T cell Workflow Log", content)))
})

test_that("initialize_logger with NULL config uses default log file", {
  initialize_logger(NULL)
  expect_true(.logger$initialized)
  expect_equal(.logger$log_file, "tcell_workflow_default.log")
  on.exit(unlink("tcell_workflow_default.log"))
})

test_that("initialize_logger with nonexistent config uses default", {
  initialize_logger("/nonexistent/config.json")
  expect_true(.logger$initialized)
  expect_equal(.logger$log_file, "tcell_workflow_default.log")
  on.exit(unlink("tcell_workflow_default.log"))
})

test_that("log_milestone writes to log file", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  log_milestone("Test Milestone", "some details")

  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "Test Milestone")
  expect_match(content, "some details")
})

test_that("log_warning writes warning entry", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  log_warning("TestCat", "something went wrong")

  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "TestCat")
  expect_match(content, "something went wrong")
})

test_that("log_error writes error entry", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  log_error("ErrorCat", "fatal problem")

  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "ERROR")
  expect_match(content, "fatal problem")
})

test_that("log_data_summary handles Seurat objects", {
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))

  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)

  counts <- matrix(sample(0:5, 50, replace = TRUE), nrow = 5, ncol = 10,
                   dimnames = list(paste0("G", 1:5), paste0("C", 1:10)))
  obj <- CreateSeuratObject(counts = as(counts, "sparseMatrix"))

  log_data_summary(obj, "Test summary")
  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "Cells: 10")
  expect_match(content, "Features: 5")
})

test_that("log_config_snapshot writes JSON content to log", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  log_config_snapshot(config_path)

  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "CONFIG SNAPSHOT")
  expect_match(content, "TestProject")
})

test_that("finalize_log writes footer with elapsed time", {
  config_path <- write_valid_config()
  on.exit(unlink(dirname(config_path), recursive = TRUE))

  initialize_logger(config_path)
  finalize_log()

  content <- paste(readLines(.logger$log_file), collapse = "\n")
  expect_match(content, "Completed")
  expect_match(content, "Total time")
})

test_that("logging functions are safe before initialization", {

  .logger$initialized <<- FALSE

  expect_silent(log_milestone("x"))
  expect_silent(log_warning("x", "y"))
  expect_silent(log_error("x", "y"))
  expect_silent(finalize_log())
})
