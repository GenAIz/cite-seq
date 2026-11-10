suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("error_codes.R")
source("portability.R")
source("disk_backed.R")
sf_error_init("add_metadata")
options(error = sf_uncaught_error_handler)

#' Map a raw add_metadata() error message onto a stable error code
#'
#' @param raw_msg Condition message of the error thrown by add_metadata().
#' @return A list with `code` (stable UPPER_SNAKE error code for programmatic
#'   handling) and `message` (user-facing text for the front-end).
#'   Unrecognized messages fall through to `PIPELINE_FAILED`.
.am_error_for <- function(raw_msg) {
  #' TRUE when `raw_msg` matches `pattern`.
  m <- function(pattern) grepl(pattern, raw_msg)
  if (m("Metadata file not found")) {
    list(code = "METADATA_FILE_MISSING",
         message = "The metadata file could not be found.")
  } else if (m("not found in the object metadata")) {
    list(code = "OBJECT_KEY_NOT_FOUND",
         message = paste("The configured object_key column does not exist in",
                         "the object's metadata. Available columns are",
                         "listed in the details."))
  } else if (m("not found in metadata file")) {
    list(code = "METADATA_KEY_NOT_FOUND",
         message = paste("The configured metadata_key column does not exist in",
                         "the metadata file. Available columns are listed in",
                         "the details."))
  } else if (m("readxl.+required for Excel")) {
    list(code = "EXCEL_SUPPORT_MISSING",
         message = paste("Excel input requires the readxl package, which is",
                         "not installed in this image. Convert the metadata to",
                         "CSV, or rebuild the image with readxl."))
  } else if (m("Unsupported file format")) {
    list(code = "UNSUPPORTED_METADATA_FORMAT",
         message = paste("The metadata file format is not supported. Use",
                         "xlsx, xls, csv, tsv or txt."))
  } else if (m("Error reading metadata file")) {
    list(code = "METADATA_LOAD_FAILED",
         message = "The metadata file exists but could not be read/parsed.")
  } else if (m("No matching entries found")) {
    list(code = "NO_METADATA_MATCHES",
         message = paste("No cells matched any row of the metadata file.",
                         "Check that object_key and metadata_key hold the",
                         "same identifiers."))
  } else if (m("changed the cell count")) {
    list(code = "METADATA_JOIN_CORRUPT",
         message = paste("The metadata join would change the number of cells,",
                         "which indicates a non-unique join key. The run was",
                         "stopped to protect the object."))
  } else {
    list(code = "PIPELINE_FAILED",
         message = "Adding metadata stopped due to an unexpected error.")
  }
}

cat("============================================================\n")
cat("  CITE-seq Add Metadata Module\n")
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

#' Resolve the metadata file from the config and the attachment port
#'
#' Resolution order: the configured path relative to `base_input`; the newest
#' attachment-port file whose name matches the configured basename; the newest
#' supported metadata file on the attachment port when the key is unset or
#' nothing matched.
#'
#' @param base_input The `$SF_INPUT_PATH` directory.
#' @param configured Value of `metadata_file` from the config, or NULL.
#' @return Path of the resolved metadata file, or NULL when nothing was found.
resolve_metadata_file <- function(base_input, configured) {
  att <- list.files(file.path(base_input, "attachment"),
                    full.names = TRUE, recursive = TRUE)
  newest <- function(paths) paths[order(file.mtime(paths), paths,
                                        decreasing = TRUE)][1]
  if (!is.null(configured) && nzchar(configured)) {
    p <- file.path(base_input, configured)
    if (file.exists(p)) return(p)
    hit <- att[basename(att) == basename(configured)]
    if (length(hit) > 0) return(newest(hit))
  }
  supported <- att[grepl("\\.(xlsx|xls|csv|tsv|txt)$", att, ignore.case = TRUE)]
  if (length(supported) > 0) return(newest(supported))
  NULL
}

#' Validate the module configuration
#'
#' Structural validation only: checks required keys, allowed values and configured paths, collecting every problem instead of stopping at the first one.
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

  if (is_missing(cfg$object_key) && is_missing(cfg$seurat_key))
    errors <- c(errors, "  [REQUIRED] 'object_key' is missing or empty.")
  check_required("metadata_key", "metadata_key (column in metadata file)")
  resolved_md <- resolve_metadata_file(base_input, cfg$metadata_file)
  if (is.null(resolved_md)) {
    errors <- c(errors, if (is_missing(cfg$metadata_file))
      "  [REQUIRED] 'metadata_file' is not set and no metadata file was found on the attachment port."
    else sprintf(
      "  [MISSING]  'metadata_file' path does not exist: %s (and no attachment-port match)",
      file.path(base_input, cfg$metadata_file)))
  } else {
    ext <- tolower(tools::file_ext(resolved_md))
    if (!(ext %in% c("xlsx", "xls", "csv", "tsv", "txt")))
      errors <- c(errors, sprintf(
        "  [INVALID]  'metadata_file' has unsupported format '.%s'. Supported formats: xlsx, xls, csv, tsv, txt.",
        ext))
  }

  check_one_of("col_types", c("text", "guess", "list"))

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

source("add_metadata.R")

output_port_dir <- file.path(base_output, "output")
dir.create(output_port_dir, recursive = TRUE, showWarnings = FALSE)

#' Locate the input Seurat .rds
#'
#' In a chained pipeline run the upstream node's output arrives on the `data`
#' input port with a timestamped filename that no config written in advance
#' can name. Resolution order: (1) the configured `input_rds` when it exists,
#' (2) the newest `.rds` found anywhere under the data port (searched
#' recursively, so the upstream dataset's layout does not matter),
#' (3) a `CONFIG_VALIDATION_FAILED` exit whose details inventory the input
#' path and point at the workflow link.
#'
#' @param base_input The `$SF_INPUT_PATH` directory.
#' @param configured The configured `input_rds` value; may be `NULL`.
#' @return Path of the `.rds` file to load. Does not return on failure.
resolve_input_rds <- function(base_input, configured) {
  if (!is.null(configured) && nzchar(configured)) {
    p <- file.path(base_input, configured)
    if (file.exists(p)) return(p)
  }
  port_dir <- file.path(base_input, "data")

  cands <- list.files(port_dir, pattern = "\\.rds$", full.names = TRUE,
                      recursive = TRUE)
  if (length(cands) > 0) {
    newest <- cands[order(file.mtime(cands), cands, decreasing = TRUE)][1]
    if (!is.null(configured) && nzchar(configured)) {
      cat(sprintf("NOTE: configured input_rds not found (%s); using newest .rds on the data port: %s\n",
                  file.path(base_input, configured), newest))
    } else {
      cat(sprintf("input_rds not set; using newest .rds on the data port: %s\n",
                  newest))
    }
    return(newest)
  }

  inv <- tryCatch(list.files(base_input, recursive = TRUE),
                  error = function(e) character(0))
  if (length(inv) > 50)
    inv <- c(inv[1:50], sprintf("... and %d more files", length(inv) - 50))
  sf_fail("CONFIG_VALIDATION_FAILED",
          "No input Seurat .rds file could be found.",
          details = c(
            if (!is.null(configured) && nzchar(configured))
              sprintf("configured input_rds does not exist: %s",
                      file.path(base_input, configured))
            else "input_rds is not set in the config",
            sprintf("no .rds files found anywhere under the data input port: %s",
                    port_dir),
            if (length(inv) == 0)
              sprintf("the input path %s contains no files at all", base_input)
            else sprintf("input file present: %s", inv),
            paste("if this node is chained from an upstream node, check the",
                  "workflow link: its SOURCE port must be the upstream",
                  "node's 'output' port (every smart function's only output",
                  "port), and the upstream run must have produced a",
                  "non-empty output dataset")))
}

input_rds <- resolve_input_rds(base_input, config$input_rds)

cat(sprintf("Loading Seurat object from: %s\n", input_rds))
sf_step("load-input")
seurat_obj <- sf_load_backed(input_rds)
cat(sprintf("  Loaded: %d cells x %d features\n\n", ncol(seurat_obj), nrow(seurat_obj)))

metadata_path <- resolve_metadata_file(base_input, config$metadata_file)
cat(sprintf("Metadata file: %s\n", metadata_path))

cat("Adding external metadata...\n")
sf_step("metadata-join")

seurat_obj <- tryCatch({
  add_metadata(
    seurat_object      = seurat_obj,
    metadata_file      = metadata_path,
    object_key         = config$object_key %||% config$seurat_key,
    metadata_key       = config$metadata_key,
    sheet              = config$sheet              %||% 1,
    col_types          = config$col_types          %||% "text",
    subset_to_existing = config$subset_to_existing %||% TRUE,
    overwrite_existing = config$overwrite_existing %||% FALSE,
    add_match_status   = config$add_match_status   %||% TRUE,
    backup_metadata    = config$backup_metadata    %||% FALSE,
    na_values          = config$na_values          %||% c("", "NA", "N/A", "na", "n/a"),
    verbose            = config$verbose            %||% TRUE
  )
}, error = function(e) {
  err <- .am_error_for(e$message)
  sf_fail(err$code, err$message, details = e$message)
})

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_base   <- sub("\\.rds$", "", config$output_name %||% "metadata", ignore.case = TRUE)
output_name   <- sprintf("%s_%s.rds", output_base, run_timestamp)
output_path   <- file.path(output_port_dir, output_name)

cat(sprintf("\nSaving Seurat object to: %s\n", output_path))
sf_step("save-output")
tryCatch(
  sf_write_backed(seurat_obj, output_path),
  error = function(e) {
    sf_fail("OUTPUT_WRITE_FAILED",
            "The output Seurat object could not be written.",
            details = c(sprintf("path: %s", output_path), e$message))
  }
)
file_size_mb <- file.size(output_path) / 1024^2
cat(sprintf("Saved (%.2f MB)\n\n", file_size_mb))

status_file <- Sys.getenv("SF_STATUS_FILE", unset = "")
if (nzchar(status_file)) writeLines("SUCCESS", status_file)

cat("============================================================\n")
cat("  Add Metadata complete.\n")
cat("============================================================\n")
