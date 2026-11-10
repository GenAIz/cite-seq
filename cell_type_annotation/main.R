suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("error_codes.R")
source("portability.R")
source("disk_backed.R")
sf_error_init("cell_type_annotation")
options(error = sf_uncaught_error_handler)

#' Map a raw cell_type_annotation() error message onto a stable error code
#'
#' @param raw_msg Condition message of the error thrown by cell_type_annotation().
#' @return A list with `code` (stable UPPER_SNAKE error code for programmatic
#'   handling) and `message` (user-facing text for the front-end).
#'   Unrecognized messages fall through to `PIPELINE_FAILED`.
.cta_error_for <- function(raw_msg) {
  #' TRUE when `raw_msg` matches `pattern`.
  m <- function(pattern) grepl(pattern, raw_msg)
  if (m("must be a Seurat object")) {
    list(code = "INPUT_NOT_SEURAT",
         message = "The input .rds file does not contain a Seurat object.")
  } else if (m("^Assay '.*' not found")) {
    list(code = "ASSAY_NOT_FOUND",
         message = "The requested assay does not exist in the input object.")
  } else if (m("Cluster column .* not found|Cluster column .* contains NA")) {
    list(code = "CLUSTER_COLUMN_INVALID",
         message = paste("The configured cluster column is missing from the",
                         "object or contains NA values. Available and",
                         "cluster-like columns are listed in the details."))
  } else if (m("Unknown reference")) {
    list(code = "UNKNOWN_REFERENCE",
         message = paste("The configured reference dataset name is not",
                         "recognised. Use monaco, hpca, blueprint_encode,",
                         "dice or immgen."))
  } else if (m("Failed to load reference|celldex|ExperimentHub")) {
    list(code = "REFERENCE_LOAD_FAILED",
         message = paste("The reference dataset could not be loaded. If this",
                         "image has no internet access, the celldex cache may",
                         "be missing."))
  } else if (m("No common genes")) {
    list(code = "NO_COMMON_GENES",
         message = paste("The input object and the reference share no gene",
                         "names. This usually means the input still uses",
                         "Ensembl IDs - run gene ID conversion in load_data."))
  } else if (m("singler_label_type|Invalid label type")) {
    list(code = "CONFIG_VALIDATION_FAILED",
         message = "singler_label_type must be 'main', 'fine', or both.")
  } else if (m("reticulate|Python environment|No suitable Python|python binary")) {
    list(code = "PYTHON_UNAVAILABLE",
         message = paste("The Python environment needed for STCAT could not",
                         "be configured in this image."))
  } else if (m("STCAT Python script not found")) {
    list(code = "STCAT_SCRIPT_MISSING",
         message = "The STCAT Python script could not be found.")
  } else if (m("STCAT")) {
    list(code = "STCAT_FAILED",
         message = "STCAT annotation failed; details carry the Python output.")
  } else if (m("T cell gating requires SingleR annotations")) {
    list(code = "SINGLER_LABELS_MISSING",
         message = paste("STCAT is only sent the cells SingleR labels as",
                         "T cells, so SingleR annotations must exist first.",
                         "Neither SingleR_main nor SingleR_fine was found in",
                         "the object - run with method 'both'."))
  } else if (m("Pseudo-bulk aggregation failed|Could not retrieve")) {
    list(code = "ANNOTATION_FAILED",
         message = "Building the per-cluster expression profiles failed.")
  } else if (m("SingleR")) {
    list(code = "ANNOTATION_FAILED",
         message = "SingleR annotation failed.")
  } else {
    list(code = "PIPELINE_FAILED",
         message = "Cell type annotation stopped due to an unexpected error.")
  }
}

cat("============================================================\n")
cat("  CITE-seq Cell Type Annotation Module\n")
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

  check_one_of("method", c("SingleR", "STCAT", "both",
                           "singler", "stcat"))

  check_one_of("reference", c(
    "monaco", "hpca", "blueprint_encode", "dice", "immgen"
  ))

  check_positive_int("singler_workers")
  check_positive_int("singler_fine_chunk_size")

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

source("cell_type_annotation.R")

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

cat("Running cell type annotation...\n")
sf_step("annotation")

.method_map <- c(singler = "SingleR", stcat = "STCAT", both = "both",
                 SingleR = "SingleR", STCAT = "STCAT")
method_norm <- unname(.method_map[config$method %||% "both"])

python_script_path <- if (!is.null(config$stcat_python_script)) {
  file.path(base_input, config$stcat_python_script)
} else {
  "python/run_stcat.py"
}

seurat_obj <- tryCatch({
  cell_type_annotation(
    seurat_object       = seurat_obj,
    method              = method_norm,
    reference           = config$reference       %||% "monaco",
    cluster_column      = config$cluster_column  %||% "RNA_clusters",
    assay               = config$assay           %||% "RNA",
    singler_label_type  = config$singler_label_type %||% c("main", "fine"),
    singler_workers     = config$singler_workers %||% 4,
    singler_fine_chunk_size = config$singler_fine_chunk_size %||% 20000L,
    stcat_python_script = python_script_path,

    stcat_batch_column  = config$stcat_batch_column %||% "orig.ident",
    prefix              = config$prefix          %||% NULL,
    fallback_column     = config$fallback_column %||% "SingleR_fine",
    verbose             = config$verbose         %||% TRUE,
    cleanup_temp        = config$cleanup_temp    %||% TRUE
  )
}, error = function(e) {
  err <- .cta_error_for(e$message)
  sf_fail(err$code, err$message, details = e$message)
})

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_base   <- sub("\\.rds$", "", config$output_name %||% "annotated_object", ignore.case = TRUE)
output_name   <- sprintf("%s_%s.rds", output_base, run_timestamp)
output_path   <- file.path(output_port_dir, output_name)

cat(sprintf("\nSaving annotated Seurat object to: %s\n", output_path))
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
cat("  Cell Type Annotation complete.\n")
cat("============================================================\n")
