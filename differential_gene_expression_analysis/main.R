suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("error_codes.R")
source("portability.R")
source("disk_backed.R")
sf_error_init("differential_gene_expression_analysis")
options(error = sf_uncaught_error_handler)

#' Map a raw perform_dge_analysis() error message onto a stable error code
#'
#' @param raw_msg Condition message of the error thrown by perform_dge_analysis().
#' @return A list with `code` (stable UPPER_SNAKE error code for programmatic
#'   handling) and `message` (user-facing text for the front-end).
#'   Unrecognized messages fall through to `PIPELINE_FAILED`.
.dge_error_for <- function(raw_msg) {
  #' TRUE when `raw_msg` matches `pattern`.
  m <- function(pattern) grepl(pattern, raw_msg)
  if (m("Invalid timepoints_to_analyze")) {
    list(code = "INVALID_TIMEPOINTS",
         message = paste("One or more requested timepoints do not exist in",
                         "the data. Available values are in the details."))
  } else if (m("Invalid conditions_to_analyze|ident\\.[12] = |reference_condition = ")) {
    list(code = "INVALID_CONDITIONS",
         message = paste("A requested condition/comparison group does not",
                         "exist in the data. Available values are in the",
                         "details."))
  } else if (m("column\\(s\\) not found in metadata|No valid cluster values|Cannot aggregate: missing columns")) {
    list(code = "COLUMN_NOT_FOUND",
         message = paste("A configured metadata column (cluster, sample,",
                         "condition or timepoint) is missing from the input",
                         "object. Details list the columns involved."))
  } else if (m("Invalid test.use")) {
    list(code = "CONFIG_VALIDATION_FAILED",
         message = "The configured test method is not supported.")
  } else if (m("require pseudobulking")) {
    list(code = "CONFIG_VALIDATION_FAILED",
         message = "Bulk methods (DESeq2/limma/edgeR) require pseudobulking = true.")
  } else {
    list(code = "PIPELINE_FAILED",
         message = "The DGE analysis stopped due to an unexpected error.")
  }
}

cat("============================================================\n")
cat("  CITE-seq Differential Gene Expression Analysis\n")
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

  #' Record a validation error when a key is negative or non-numeric.
  check_non_negative <- function(field, label = field) {
    val <- cfg[[field]]
    if (!is.null(val) && (!is.numeric(val) || val < 0))
      errors <<- c(errors, sprintf(
        "  [INVALID]  '%s' must be >= 0 (got: %s).", label, val
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

  check_required("cluster_col",   "cluster_col (metadata column with cell type labels)")

  if (isTRUE(cfg$pseudobulking)) {
    if (is_missing(cfg$sample_col))
      errors <- c(errors,
        "  [REQUIRED] 'sample_col' is required when 'pseudobulking' is true.")
    if (is_missing(cfg$condition_col))
      errors <- c(errors,
        "  [REQUIRED] 'condition_col' is required when 'pseudobulking' is true.")
  }

  if (is_missing(cfg$ident_1 %||% cfg[["ident.1"]]) &&
      !is_missing(cfg$ident_2 %||% cfg[["ident.2"]]))
    errors <- c(errors,
      "  [INVALID]  'ident_2' requires 'ident_1' to be set.")

  tu <- cfg$test_use
  tu_allowed <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson",
                  "lr", "mast", "deseq2", "limma", "edger")
  if (!is.null(tu) && !(tolower(tu) %in% tu_allowed))
    errors <- c(errors, sprintf(
      "  [INVALID]  'test_use' = '%s'. Allowed values (any case): %s.",
      tu, paste(c("MAST", "wilcox", "bimod", "roc", "t", "negbinom",
                  "poisson", "LR", "DESeq2", "limma", "edgeR"), collapse = ", ")))

  check_one_of("pseudobulk_mode",
               c("cluster_sample", "sample_only", "celltype_condition"))

  check_one_of("design_strategy", c("stratified", "interaction"))

  check_one_of("p_adjust_method",
               c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr", "none"))

  check_positive_int("min_samples_per_group")
  check_positive_int("min_cells_per_pseudobulk")
  check_non_negative("lfc_threshold")

  if (!is.null(cfg$output_rds) || !is.null(cfg$output_csv))
    errors <- c(errors, paste0(
      "  [RETIRED]  'output_rds'/'output_csv' are retired: this module ",
      "writes ONE results CSV named by 'output_name' (plus a fixed-name ",
      "pass-through copy of the input object)."))

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

source("dge.R")

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
seurat_obj <- sf_load_backed(input_rds, keep_backed = NULL)
cat(sprintf("  Loaded: %d cells x %d features\n\n", ncol(seurat_obj), nrow(seurat_obj)))

.dge_assay <- tryCatch(methods::slot(seurat_obj, "active.assay"), error = function(e) NULL)
if (is.null(.dge_assay) || !length(.dge_assay) || !.dge_assay %in% names(seurat_obj@assays))
  .dge_assay <- names(seurat_obj@assays)[1]
.dge_layers <- tryCatch(.sf_assay_layers(seurat_obj@assays[[.dge_assay]]),
                        error = function(e) list())
if (is.null(.dge_layers[["data"]])) {
  sf_fail("INPUT_MISSING_NORMALIZED_DATA",
          paste("The input object has no log-normalized 'data' layer, which every",
                "differential expression test in this module runs on."),
          details = c(
            sprintf("assay checked: %s", .dge_assay),
            sprintf("layers present: %s",
                    if (length(.dge_layers)) paste(names(.dge_layers), collapse = ", ")
                    else "none"),
            paste("objects ingested from .h5ad/.h5mu carry only the file's X matrix",
                  "as 'counts'; normalize before export (e.g. Seurat::NormalizeData)",
                  "or run this node on an object whose 'data' layer is populated")))
}

cat("Running differential gene expression analysis...\n")
sf_step("dge-analysis")
cat(sprintf("  Method         : %s\n", config$test_use   %||% "MAST"))
cat(sprintf("  Pseudobulking  : %s\n", config$pseudobulking %||% FALSE))
cat(sprintf("  Cluster column : %s\n", config$cluster_col %||% "SingleR_main"))
cat(sprintf("  Condition col  : %s\n\n", config$condition_col %||% "Terminal Group"))

dge_results <- tryCatch({
  perform_dge_analysis(
    seurat_obj              = seurat_obj,
    test.use                = config$test_use             %||% "MAST",
    pseudobulking           = config$pseudobulking         %||% FALSE,
    pseudobulk_mode         = config$pseudobulk_mode      %||% "cluster_sample",
    design_strategy         = config$design_strategy      %||% "stratified",
    cluster_col             = config$cluster_col          %||% "SingleR_main",
    sample_col              = config$sample_col           %||% "Participant ID",
    condition_col           = config$condition_col        %||% "Terminal Group",
    timepoint_col           = config$timepoint_col        %||% "Timepoint (months)",
    participant_col         = config$participant_col      %||% config$sample_col %||% "Participant ID",
    pseudobulk_assay        = config$pseudobulk_assay     %||% "RNA",
    reference_condition     = config$reference_condition  %||% NULL,

    timepoints_to_analyze   = config$timepoints_to_analyze %||% NULL,
    conditions_to_analyze   = config$conditions_to_analyze %||% NULL,
    ident.1                 = config$ident_1 %||% config[["ident.1"]] %||% NULL,
    ident.2                 = config$ident_2 %||% config[["ident.2"]] %||% NULL,
    contrasts               = config$contrasts            %||% NULL,
    min_samples_per_group   = config$min_samples_per_group %||% 2,
    min_cells_per_pseudobulk = config$min_cells_per_pseudobulk %||% 10,
    lfc_threshold           = config$lfc_threshold        %||% 0,
    p_adjust_method         = config$p_adjust_method      %||% "BH",
    verbose                 = config$verbose              %||% TRUE
  )
}, error = function(e) {
  err <- .dge_error_for(e$message)
  sf_fail(err$code, err$message, details = e$message)
})

cat(sprintf("  DGE complete: %d rows in results table.\n\n", nrow(dge_results)))

if (nrow(dge_results) == 0) {
  sf_fail("NO_DGE_RESULTS",
          paste("The analysis ran but produced no results: every",
                "cluster/contrast was skipped. Check group sizes against",
                "min_samples_per_group and the condition/timepoint filters."),
          details = c(
            sprintf("min_samples_per_group: %s", config$min_samples_per_group %||% 2),
            sprintf("test_use: %s", config$test_use %||% "MAST"),
            sprintf("design_strategy: %s", config$design_strategy %||% "stratified")))
}

run_timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
csv_base       <- sub("\\.csv$", "", config$output_name %||% "dge_results", ignore.case = TRUE)
output_csv     <- file.path(output_port_dir, sprintf("%s_%s.csv", csv_base, run_timestamp))

cat(sprintf("Saving results to: %s\n", output_csv))
sf_step("save-output")
tryCatch(
  write.csv(dge_results, file = output_csv, row.names = FALSE),
  error = function(e) {
    sf_fail("OUTPUT_WRITE_FAILED",
            "The results CSV could not be written.",
            details = c(sprintf("path: %s", output_csv), e$message))
  }
)

passthrough_rds <- file.path(output_port_dir,
                             sprintf("dge_passthrough_%s.rds", run_timestamp))
cat(sprintf("Passing input object through to: %s\n", passthrough_rds))
if (!sf_copy_bundle(input_rds, passthrough_rds)) {
  sf_fail("OUTPUT_WRITE_FAILED",
          "The pass-through copy of the input object could not be written.",
          details = c(sprintf("from: %s", input_rds),
                      sprintf("to: %s", passthrough_rds)))
}

status_file <- Sys.getenv("SF_STATUS_FILE", unset = "")
if (nzchar(status_file)) writeLines("SUCCESS", status_file)

cat("\n============================================================\n")
cat("  DGE Analysis complete.\n")
cat("============================================================\n")
