suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("error_codes.R")
source("portability.R")
source("disk_backed.R")
sf_error_init("load_data")
options(error = sf_uncaught_error_handler)

#' Map a raw load_data() error message onto a stable error code
#'
#' @param raw_msg Condition message of the error thrown by load_data().
#' @return A list with `code` (stable UPPER_SNAKE error code for programmatic
#'   handling) and `message` (user-facing text for the front-end).
#'   Unrecognized messages fall through to `PIPELINE_FAILED`.
.ld_error_for <- function(raw_msg) {
  #' TRUE when `raw_msg` matches `pattern`.
  m <- function(pattern) grepl(pattern, raw_msg)
  if (m("Unsupported input format")) {
    list(code = "UNSUPPORTED_FORMAT",
         message = paste("The input file format is not supported. This module",
                         "loads ONE pre-processed object: .rds, .h5ad or .h5mu."))
  } else if (m("not a Seurat object")) {
    list(code = "INPUT_NOT_SEURAT",
         message = paste("The input .rds does not contain a Seurat object -",
                         "it holds unprocessed data that has not been through",
                         "QC, normalisation or clustering."))
  } else if (m("square X matrix")) {
    list(code = "H5AD_SQUARE_MATRIX",
         message = paste("The .h5ad file has a square cells-by-genes matrix,",
                         "which the reader cannot process safely. Subset the",
                         "input so cells != genes, or convert to 10x format."))
  } else if (m("MuData|h5mu|ADT/PROTEIN modality")) {
    list(code = "H5MU_PARSE_FAILED",
         message = "The .h5mu file could not be parsed into a Seurat object.")
  } else if (m("h5ad|anndata|reticulate|python")) {
    list(code = "H5AD_PARSE_FAILED",
         message = "The .h5ad file could not be parsed into a Seurat object.")
  } else if (m("Unsupported species")) {
    list(code = "UNSUPPORTED_SPECIES",
         message = "Gene ID conversion supports only 'human' and 'mouse'.")
  } else if (m("Required package")) {
    list(code = "GENE_ID_CONVERSION_FAILED",
         message = paste("Gene ID conversion could not run because a required",
                         "annotation package is missing from the image."))
  } else if (m("cannot open|No such file")) {
    list(code = "INPUT_LOAD_FAILED",
         message = "The input file could not be read.")
  } else {
    list(code = "PIPELINE_FAILED",
         message = "Loading the data stopped due to an unexpected error.")
  }
}

cat("============================================================\n")
cat("  CITE-seq Load Data Module\n")
cat("============================================================\n\n")

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to JSON configuration file", metavar = "FILE")
)

args <- parse_args(OptionParser(option_list = option_list))

sf_step("configuration")

if (is.null(args$config)) {
  sf_fail("CONFIG_MISSING",
          "No configuration file was provided to the module.",
          details = "Usage: Rscript main.R --config configs/config.json")
}

if (!file.exists(args$config)) {
  sf_fail("CONFIG_MISSING",
          "The configuration file could not be found.",
          details = sprintf("path: %s", args$config))
}

config <- tryCatch(
  fromJSON(args$config),
  error = function(e) {
    sf_fail("CONFIG_PARSE_FAILED",
            "The configuration file is not valid JSON.",
            details = c(sprintf("path: %s", args$config), e$message))
  }
)

base_input  <- Sys.getenv("SF_INPUT_PATH",  unset = "run/in")
base_output <- Sys.getenv("SF_OUTPUT_PATH", unset = "run/out")

#' Validate the module configuration
#'
#' Structural validation only: checks required keys, allowed values and
#' configured paths, collecting every problem instead of stopping at the
#' first one.
#'
#' @param cfg Parsed configuration list from `configs/config.json`.
#' @param base_input The `$SF_INPUT_PATH` directory, used to resolve
#'   configured file paths.
#' @return Character vector of human-readable validation errors; empty when
#'   the configuration is valid.
validate_config <- function(cfg, base_input) {
  errors <- character(0)

  #' TRUE when a config value is absent or an empty/whitespace string.
  is_missing <- function(x) is.null(x) || (is.character(x) && !nzchar(trimws(x)))

  #' Record a validation error when a required key is missing or empty.
  check_required <- function(field, label = field) {
    if (is_missing(cfg[[field]]))
      errors <<- c(errors, sprintf("  [REQUIRED] '%s' is missing or empty.", label))
  }

  #' Record a validation error when a key's value is outside `allowed`.
  check_one_of <- function(field, allowed, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && !(val %in% allowed))
      errors <<- c(errors, sprintf(
        "  [INVALID]  '%s' = '%s'. Allowed values: %s.",
        label, val, paste(allowed, collapse = ", ")
      ))
  }

  #' Record a validation error when a key is not a positive integer.
  check_positive_int <- function(field, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && (!is.numeric(val) || val < 1 || val != as.integer(val)))
      errors <<- c(errors, sprintf(
        "  [INVALID]  '%s' must be a positive integer (got: %s).", label, val
      ))
  }

  #' Record a validation error when a configured file path does not exist.
  check_file_exists <- function(field, base = base_input, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && nzchar(val)) {
      path <- file.path(base, val)
      if (!file.exists(path))
        errors <<- c(errors, sprintf(
          "  [MISSING]  '%s' path does not exist: %s", label, path
        ))
    }
  }

  #' Record a validation error when a configured path exists as neither a directory nor a file.
  check_dir_exists <- function(field, base = base_input, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && nzchar(val)) {
      path <- file.path(base, val)
      if (!dir.exists(path) && !file.exists(path))
        errors <<- c(errors, sprintf(
          "  [MISSING]  '%s' path does not exist: %s", label, path
        ))
    }
  }

  check_required("load_data_path",
                 "load_data_path (path to the pre-processed input object)")
  check_dir_exists("load_data_path")

  check_one_of("species_name", c("human", "mouse"))

  errors
}

cat("Validating configuration...\n")
sf_step("validation")
validation_errors <- validate_config(config, base_input)

if (length(validation_errors) > 0) {
  sf_fail("CONFIG_VALIDATION_FAILED",
          "The configuration has problems that must be fixed before the module can run.",
          details = trimws(validation_errors))
}

cat("  Configuration OK.\n\n")

source("load_data.R")

output_port_dir <- file.path(base_output, "output")
dir.create(output_port_dir, recursive = TRUE, showWarnings = FALSE)

load_path <- file.path(base_input, config$load_data_path)

#' Does this input need the Python HDF5 backend?
#'
#' `.h5ad`/`.h5mu` reading depends on Python imported in-process through
#' reticulate. A broken linkage previously degraded silently: a CITE-seq h5ad
#' loaded as RNA-only and the run reported success with the protein data
#' dropped. For HDF5-based inputs the backend is therefore probed up front,
#' the same way it is used, and the run stops loudly if it is broken.
#'
#' @param path Configured input path.
#' @return `TRUE` when `path` is an `.h5ad`/`.h5mu` file or a directory
#'   containing one.
.h5_backend_needed <- function(path) {
  if (grepl("\\.h5(ad|mu)$", path, ignore.case = TRUE)) return(TRUE)
  isTRUE(dir.exists(path)) &&
    length(list.files(path, pattern = "\\.h5(ad|mu)$", ignore.case = TRUE)) > 0
}

if (.h5_backend_needed(load_path)) {
  sf_step("python-backend-check")
  backend_err <- tryCatch({
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("R package 'reticulate' is not installed")
    reticulate::py_run_string("import h5py, anndata")
    NULL
  }, error = function(e) conditionMessage(e))

  if (!is.null(backend_err)) {
    sf_fail("H5AD_PYTHON_BACKEND_UNAVAILABLE",
            paste("The Python backend needed to read .h5ad/.h5mu files is not",
                  "working in this environment, so the file cannot be read",
                  "completely. Loading anyway could silently drop the protein",
                  "(ADT) data, so the run was stopped."),
            details = c(
              paste("python:",
                    Sys.getenv("RETICULATE_PYTHON", unset = "reticulate default")),
              backend_err))
  }
  cat("  Python h5ad backend OK.\n")
}

cat("Step 1: Loading CITE-seq data...\n")
sf_step("load-data")
.t_phase <- Sys.time()
if (file.exists(load_path)) {
  cat(sprintf(paste0(
    "  Input file size: %.2f GB - reading a multi-GB object can take several\n",
    "  minutes with no further output; if the log ENDS here, the container\n",
    "  most likely ran out of memory while inflating the object\n"),
    file.size(load_path) / 1024^3))
}

seurat_obj <- tryCatch({
  load_data(
    rna_path     = load_path,
    project_name = config$project_name %||% "CITEseqProject",
    species_name = config$species_name %||% "human",
    verbose      = config$verbose      %||% TRUE
  )
}, error = function(e) {
  err <- .ld_error_for(e$message)
  sf_fail(err$code, err$message,
          details = c(sprintf("input: %s", load_path), e$message))
})

cat(sprintf("  Loaded: %d cells x %d features in %.1f min\n\n",
            ncol(seurat_obj), nrow(seurat_obj),
            as.numeric(difftime(Sys.time(), .t_phase, units = "mins"))))

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_base   <- sub("\\.rds$", "", config$output_name %||% "loaded_data", ignore.case = TRUE)
output_name   <- sprintf("%s_%s.rds", output_base, run_timestamp)
output_path   <- file.path(output_port_dir, output_name)

cat(sprintf("Saving Seurat object to: %s (multi-GB objects take a few minutes)\n",
            output_path))
sf_step("save-output")
.t_phase <- Sys.time()
tryCatch(
  sf_write_backed(seurat_obj, output_path),
  error = function(e) {
    sf_fail("OUTPUT_WRITE_FAILED",
            "The output Seurat object could not be written.",
            details = c(sprintf("path: %s", output_path), e$message))
  }
)
file_size_mb <- file.size(output_path) / 1024^2
cat(sprintf("Saved (%.2f MB) in %.1f min\n\n", file_size_mb,
            as.numeric(difftime(Sys.time(), .t_phase, units = "mins"))))

status_file <- Sys.getenv("SF_STATUS_FILE", unset = "")
if (nzchar(status_file)) {
  writeLines("SUCCESS", status_file)
}

cat("============================================================\n")
cat("  Load Data complete.\n")
cat("============================================================\n")
