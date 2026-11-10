suppressPackageStartupMessages({
  library(testthat)
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!testthat::is_testing()) {

.script_path <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) normalizePath(sub("^--file=", "", fa[1]))
  else tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
})
this_dir <- if (!is.na(.script_path)) dirname(.script_path) else normalizePath(getwd())
sf_root <- if (basename(this_dir) == "tests") dirname(this_dir) else this_dir
tests_dir <- file.path(sf_root, "tests")

source(file.path(sf_root, "load_data.R"))

cat("\n======================================================\n")
cat(" Load Data Smart Function - Test Suite\n")
cat("======================================================\n\n")

result <- test_dir(
  tests_dir,
  reporter = "progress",
  stop_on_failure = FALSE,
  load_helpers = TRUE
)

cat("\n")
df <- as.data.frame(result)
n_fail <- sum(df$failed, na.rm = TRUE)
n_error <- sum(df$error, na.rm = TRUE)
n_pass <- sum(df$passed, na.rm = TRUE)
n_skip <- sum(df$skipped, na.rm = TRUE)

cat(sprintf("Results: %d passed | %d failed | %d errors | %d skipped\n",
            n_pass, n_fail, n_error, n_skip))

if (n_fail > 0 || n_error > 0) quit(status = 1) else quit(status = 0)
}
