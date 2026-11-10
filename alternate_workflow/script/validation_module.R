validate_config_schema <- function(config) {
  required_fields <- c("metadata", "paths", "execution", "load_data", "processing", "aggregation")
  missing <- setdiff(required_fields, names(config))
  
  if (length(missing) > 0) {
    return(list(valid = FALSE, message = paste("Missing required fields:", paste(missing, collapse = ", "))))
  }
  
  return(list(valid = TRUE, message = "Schema valid"))
}

validate_input_files <- function(config) {
  errors <- character(0)

  check_readable <- function(path, label) {
    if (is.null(path) || !isTRUE(nzchar(path))) return()
    if (!file.exists(path)) {
      errors <<- c(errors, sprintf("  [MISSING]  %s not found: %s", label, path))
    } else if (file.access(path, mode = 4) != 0) {
      errors <<- c(errors, sprintf("  [DENIED]   No read permission on %s: %s", label, path))
    }
  }

  check_readable(config$paths$input_data$rna_path, "RNA file")
  check_readable(config$paths$input_data$adt_path,  "ADT file")

  if (isTRUE(config$add_metadata$enabled))
    check_readable(config$paths$input_data$metadata_file, "Metadata file")

  if (length(errors) > 0)
    return(list(valid = FALSE, message = paste(c("Input file problems:", errors), collapse = "\n")))

  return(list(valid = TRUE, message = "All input files readable"))
}

validate_metadata_columns <- function(config) {
  if (!isTRUE(config$add_metadata$enabled)) {
    return(list(valid = TRUE, message = "Metadata validation skipped"))
  }
  
  if (is.null(config$paths$input_data$metadata_file)) {
    return(list(valid = FALSE, message = "Metadata file path not specified"))
  }
  
  return(list(valid = TRUE, message = "Metadata columns configured"))
}

validate_parameter_ranges <- function(config) {

  qc <- config$processing$qc_filters
  if (!is.null(qc$min_features) && qc$min_features < 0) {
    return(list(valid = FALSE, message = "min_features must be >= 0"))
  }
  if (!is.null(qc$max_features) && !is.null(qc$min_features) && qc$max_features < qc$min_features) {
    return(list(valid = FALSE, message = "max_features must be >= min_features"))
  }
  if (!is.null(qc$max_percent_mt) && (qc$max_percent_mt < 0 || qc$max_percent_mt > 100)) {
    return(list(valid = FALSE, message = "max_percent_mt must be between 0 and 100"))
  }
  
  if (!is.null(config$processing$clustering$rna_resolution) && config$processing$clustering$rna_resolution <= 0) {
    return(list(valid = FALSE, message = "rna_resolution must be > 0"))
  }
  
  if (!is.null(config$execution$n_cores) && config$execution$n_cores < 1) {
    return(list(valid = FALSE, message = "n_cores must be >= 1"))
  }
  
  return(list(valid = TRUE, message = "Parameter ranges valid"))
}

ensure_writable_directory <- function(dir_path, label) {
  if (is.null(dir_path) || !isTRUE(nzchar(dir_path)))
    return(sprintf("  [MISSING]  %s path not specified", label))

  if (!dir.exists(dir_path)) {
    created <- tryCatch({
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      TRUE
    }, error = function(e) FALSE)
    if (!created)
      return(sprintf("  [DENIED]   Cannot create %s: %s", label, dir_path))
  }

  test_file <- file.path(dir_path, ".write_test")
  result <- tryCatch({
    write("", file = test_file)
    unlink(test_file)
    NULL
  }, error = function(e) e$message)

  if (!is.null(result))
    return(sprintf("  [DENIED]   No write permission on %s: %s", label, dir_path))

  return(NULL)
}

ensure_output_directories <- function(config) {
  errors <- character(0)

  out_err <- ensure_writable_directory(config$paths$output_dir, "output directory")
  if (!is.null(out_err)) errors <- c(errors, out_err)

  logs_dir <- Sys.getenv("SF_LOG_DIR", unset = "")
  if (!nzchar(logs_dir)) logs_dir <- config$paths$logs_dir %||% "logs"
  log_err <- ensure_writable_directory(logs_dir, "logs directory")
  if (!is.null(log_err)) errors <- c(errors, log_err)

  if (length(errors) > 0)
    return(list(valid = FALSE, message = paste(c("Directory permission problems:", errors), collapse = "\n")))

  return(list(valid = TRUE, message = sprintf(
    "Output directories ready:\n  output: %s\n  logs:   %s",
    config$paths$output_dir, logs_dir
  )))
}

validate_workflow_config <- function(config) {
  checks <- list(
    schema    = validate_config_schema(config),
    inputs    = validate_input_files(config),
    metadata  = validate_metadata_columns(config),
    params    = validate_parameter_ranges(config),
    dirs      = ensure_output_directories(config)
  )

  failures <- Filter(function(r) !r$valid, checks)

  if (length(failures) > 0) {
    combined <- paste(sapply(failures, `[[`, "message"), collapse = "\n\n")
    return(list(valid = FALSE, message = combined))
  }

  return(list(valid = TRUE, message = "All validation checks passed"))
}

validate_step_requirements <- function(config) {
  steps <- config$execution$steps
  if (is.null(steps) || length(steps) == 0)
    return(list(valid = TRUE, message = "No steps to check"))

  have <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))

  hard <- list(
    load_data            = c("Seurat", "SeuratObject"),
    add_metadata         = c("Seurat"),
    process_seurat       = c("Seurat", "SeuratObject"),
    cell_type_annotation = c("Seurat", "SingleR", "celldex"),
    aggregate_ml         = c("Seurat", "openxlsx"),
    export_loupe         = c("Seurat", "loupeR")
  )

  soft <- list(
    load_data            = c("BPCells", "reticulate"),
    process_seurat       = c("BPCells", "harmony", "transport", "scDblFinder"),
    cell_type_annotation = c("reticulate")
  )

  errors <- character(0)
  warns  <- character(0)

  for (s in steps) {
    for (p in hard[[s]] %||% character(0))
      if (!have(p))
        errors <- c(errors, sprintf(
          "  [MISSING] step '%s' requires R package '%s'", s, p))
    for (p in soft[[s]] %||% character(0))
      if (!have(p))
        warns <- c(warns, sprintf(
          "  [optional] step '%s' works better with '%s' (not installed)", s, p))
  }

  if ("process_seurat" %in% steps) {
    adt_on <- !identical(config$batch$correct_adt, FALSE) &&
              !identical(config$batch$correct, FALSE)
    if (isTRUE(adt_on) && !have("ADTnorm"))
      warns <- c(warns,
        "  [optional] ADT batch correction is enabled but 'ADTnorm' is not installed")
  }

  if (length(warns) > 0)
    message(paste(c("Step prerequisite warnings:", warns), collapse = "\n"))

  if (length(errors) > 0)
    return(list(valid = FALSE, message = paste(c(
      "Step prerequisites not met (checked before execution to avoid failing hours in):",
      errors,
      "  Rebuild the Docker image, or remove the affected step from execution.steps."),
      collapse = "\n")))

  list(valid = TRUE, message = sprintf(
    "All prerequisites present for %d step(s)", length(steps)))
}
