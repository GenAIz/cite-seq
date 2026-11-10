suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)

  library(Seurat)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("error_codes.R")
source("portability.R")
source("disk_backed.R")
sf_error_init("functional_state_profile")
options(error = sf_uncaught_error_handler)

#' Map a raw classify_functional_states_ucell() error message onto a stable error code
#'
#' @param raw_msg Condition message of the error thrown by classify_functional_states_ucell().
#' @return A list with `code` (stable UPPER_SNAKE error code for programmatic
#'   handling) and `message` (user-facing text for the front-end).
#'   Unrecognized messages fall through to `PIPELINE_FAILED`.
.fsp_error_for <- function(raw_msg) {
  #' TRUE when `raw_msg` matches `pattern`.
  m <- function(pattern) grepl(pattern, raw_msg, fixed = FALSE)
  if (m("must be a Seurat object")) {
    list(code = "INPUT_NOT_SEURAT",
         message = "The input .rds file does not contain a Seurat object.")
  } else if (m("marker_json_path file not found")) {
    list(code = "MARKER_JSON_MISSING",
         message = "The marker definitions JSON file could not be found.")
  } else if (m("No marker genes survive")) {
    list(code = "MARKER_STATE_EMPTY",
         message = paste("One or more functional states have no marker genes",
                         "present in the input data. This usually means the",
                         "input uses Ensembl IDs instead of gene symbols, or",
                         "an over-restrictive HVG filter."))
  } else if (m("^Assay '") || m("not found in Seurat object")) {
    list(code = "ASSAY_NOT_FOUND",
         message = "The requested assay does not exist in the input object.")
  } else if (m("UCell")) {
    list(code = "SCORING_FAILED",
         message = "UCell signature scoring failed.")
  } else {
    list(code = "PIPELINE_FAILED",
         message = "The functional state analysis stopped due to an unexpected error.")
  }
}

cat("============================================================\n")
cat("  CITE-seq Functional State Profile Module\n")
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

#' Resolve the marker JSON from the config, the attachment port, or the bundle
#'
#' Resolution order: absolute configured paths as given; the configured path
#' relative to `base_input`; the newest attachment-port file whose name
#' matches the configured basename; the newest `.json` on the attachment
#' port; the configured value as an image-relative path; when unset, the
#' newest `.json` on the attachment port, then the bundled
#' `cell_state_markers.json`.
#'
#' @param base_input The `$SF_INPUT_PATH` directory.
#' @param configured Value of `marker_json_path` from the config, or NULL.
#' @return Path of the marker JSON to load.
resolve_marker_json <- function(base_input, configured) {
  att <- list.files(file.path(base_input, "attachment"),
                    full.names = TRUE, recursive = TRUE)
  newest <- function(paths) paths[order(file.mtime(paths), paths,
                                        decreasing = TRUE)][1]
  jsons <- att[grepl("\\.json$", att, ignore.case = TRUE)]
  if (!is.null(configured) && nzchar(configured)) {
    if (grepl("^(/|[A-Za-z]:)", configured)) return(configured)
    p <- file.path(base_input, configured)
    if (file.exists(p)) return(p)
    hit <- att[basename(att) == basename(configured)]
    if (length(hit) > 0) return(newest(hit))
    if (length(jsons) > 0) return(newest(jsons))
    return(configured)
  }
  if (length(jsons) > 0) return(newest(jsons))
  "cell_state_markers.json"
}

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

  #' Record a validation error when a key is not a positive integer.
  check_positive_int <- function(field, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && (!is.numeric(val) || val < 1 || val != as.integer(val)))
      errors <<- c(errors, sprintf(
        "  [INVALID]  '%s' must be a positive integer (got: %s).", label, val
      ))
  }

  check_positive_int("maxRank")
  check_positive_int("ncores")

  ms <- cfg$min_score
  if (!is.null(ms) && (!is.numeric(ms) || length(ms) != 1 || is.na(ms) ||
                       ms < 0 || ms > 1)) {
    errors <- c(errors, sprintf(
      "  [INVALID]  'min_score' must be a number between 0 and 1 (got: %s).",
      paste(ms, collapse = ", ")
    ))
  }

  if (!is.null(cfg$marker_json_path) && nzchar(cfg$marker_json_path)) {
    resolved <- resolve_marker_json(base_input, cfg$marker_json_path)

    if (!file.exists(resolved))
      errors <- c(errors, sprintf(
        "  [MISSING]  'marker_json_path' does not exist: %s (also not found on the attachment port)",
        cfg$marker_json_path
      ))
  }

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

source("classify_functional_state.R")

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

marker_json <- resolve_marker_json(base_input, config$marker_json_path)

cat(sprintf("Marker definitions: %s\n\n", marker_json))

cat("Running functional state classification (UCell)...\n")
sf_step("classification")

seurat_obj <- tryCatch({
  classify_functional_states_ucell(
    seurat_obj       = seurat_obj,
    marker_json_path = marker_json,
    use_core_only    = config$use_core_only %||% FALSE,
    use_hvg          = config$use_hvg       %||% FALSE,
    assay            = config$assay         %||% "RNA",
    maxRank          = config$maxRank       %||% 1500,
    min_score        = config$min_score,
    ncores           = config$ncores        %||% 1,
    verbose          = config$verbose       %||% TRUE
  )
}, error = function(e) {
  err <- .fsp_error_for(e$message)
  sf_fail(err$code, err$message, details = e$message)
})

if (config$run_validation %||% TRUE) {
  cat("\nValidating UCell scores...\n")
  validate_classification(seurat_obj)
}

run_timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_base       <- sub("\\.rds$", "", config$output_rds %||% "functional_state_object", ignore.case = TRUE)
csv_base       <- sub("\\.csv$", "", config$output_csv %||% "functional_state_scores",  ignore.case = TRUE)
output_rds     <- file.path(output_port_dir, sprintf("%s_%s.rds", rds_base, run_timestamp))
output_csv     <- file.path(output_port_dir, sprintf("%s_%s.csv", csv_base, run_timestamp))

cat(sprintf("\nSaving Seurat object to: %s\n", output_rds))
sf_step("save-output")
tryCatch(
  sf_write_backed(seurat_obj, output_rds),
  error = function(e) {
    sf_fail("OUTPUT_WRITE_FAILED",
            "The output Seurat object could not be written.",
            details = c(sprintf("path: %s", output_rds), e$message))
  }
)

marker_states <- tryCatch({
  mk <- fromJSON(marker_json)
  names(mk)[!startsWith(names(mk), "_")]
}, error = function(e) character(0))

state_cols <- c("dominant_state", "dominant_state_score",
                "dominant_state_margin", marker_states)
available_cols <- intersect(state_cols, colnames(seurat_obj@meta.data))

if (length(available_cols) > 0) {
  score_df <- seurat_obj@meta.data[, available_cols, drop = FALSE]
  score_df <- cbind(cell_barcode = rownames(score_df), score_df)
  cat(sprintf("Exporting scores to: %s\n", output_csv))
  tryCatch(
    write.csv(score_df, file = output_csv, row.names = FALSE),
    error = function(e) {
      sf_fail("OUTPUT_WRITE_FAILED",
              "The per-cell score CSV could not be written.",
              details = c(sprintf("path: %s", output_csv), e$message))
    }
  )
}

status_file <- Sys.getenv("SF_STATUS_FILE", unset = "")
if (nzchar(status_file)) writeLines("SUCCESS", status_file)

cat("\n============================================================\n")
cat("  Functional State Profile complete.\n")
cat("============================================================\n")
