library(optparse, quietly = TRUE)
library(jsonlite, quietly = TRUE)

cat("+================================================================+\n")
cat("| T cell Workflow Orchestrator v3.0 (Dispatch-Loop) |\n")
cat("| Loading configuration and validating pipeline... |\n")
cat("+================================================================+\n\n")

option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to configuration file (JSON)", metavar = "FILE"),
  make_option(c("--dry-run"), action = "store_true", default = FALSE,
              help = "Validate configuration and print the plan without executing"),
  make_option(c("--validate"), action = "store_true", default = FALSE,
              help = "Validate configuration and input files, then exit"),
  make_option(c("--resume"), type = "character", default = NULL,
              help = "Resume from specified step", metavar = "STEP"),
  make_option(c("--debug"), action = "store_true", default = FALSE,
              help = "Enable debug mode with verbose output")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

source("script/logging_module.R", local = FALSE)
source("script/validation_module.R", local = FALSE)

initialize_logger(args$config)

start_log_tee()
log_milestone("Orchestrator Started", "v3.0 - Dispatch-Loop Architecture")

if (is.null(args$config)) {
  log_error("Configuration", "No config file specified. Use: --config config.json")
  finalize_log()
  quit(status = 1)
}

if (!file.exists(args$config)) {
  log_error("Configuration", sprintf("Config file not found: %s", args$config))
  finalize_log()
  quit(status = 1)
}

tryCatch({
  config <- fromJSON(args$config)
  log_milestone("Configuration Loaded", args$config)
}, error = function(e) {
  log_error("JSON Parse", sprintf("Invalid JSON: %s", e$message))
  finalize_log()
  quit(status = 1)
})

log_milestone("Validation Started", "Running pre-execution checks...")
validation_result <- validate_workflow_config(config)

if (!validation_result$valid) {
  log_error("Validation Failed", validation_result$message)
  finalize_log()
  quit(status = 1)
}

step_req <- validate_step_requirements(config)
if (!step_req$valid) {
  log_error("Step Prerequisites", step_req$message)
  finalize_log()
  quit(status = 1)
}
log_milestone("Step Prerequisites", step_req$message)

log_milestone("Validation Complete", "All checks passed")
log_config_snapshot(args$config)

if (isTRUE(args$validate)) {
  cat("\nConfiguration and input files are valid.\n")
  cat("Config: ", args$config, "\n", sep = "")
  cat("Steps: ", paste(config$execution$steps %||% "(default)", collapse = " -> "), "\n\n", sep = "")
  finalize_log()
  quit(status = 0)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

step_output_path <- function(step_key, default_name) {
  name <- config$output_files[[step_key]] %||% default_name
  base <- sub("\\.rds$", "", name, ignore.case = TRUE)
  file.path(config$paths$output_dir, sprintf("%s_%s.rds", base, run_timestamp))
}

step_resume_path <- function(step_key, default_name) {
  name <- config$output_files[[step_key]] %||% default_name
  base <- sub("\\.rds$", "", name, ignore.case = TRUE)
  out_dir <- config$paths$output_dir
  pattern <- sprintf("^%s_[0-9]{8}_[0-9]{6}\\.rds$", base)
  candidates <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  if (length(candidates) > 0) {
    return(candidates[which.max(file.mtime(candidates))])
  }
  file.path(out_dir, paste0(base, ".rds"))
}

current_object <- NULL
generated_files <- character(0)

load_previous_output <- function(current_step_idx, executed_steps) {

  all_step_outputs <- list(
    load_data = c("load_data", "load_data.rds"),
    add_metadata = c("add_metadata", "add_metadata.rds"),
    process_seurat = c("process_seurat", "process_seurat.rds"),
    cell_type_annotation = c("cell_type_annotation", "cell_type_annotation.rds")

  )

  chained_outputs <- list(
    cell_type_annotation = c("cell_type_annotation", "annotated_seurat.rds")
  )

  canonical <- c("load_data", "add_metadata", "process_seurat",
                 "cell_type_annotation", "aggregate_ml")

  prior_steps <- if (current_step_idx > 1) rev(executed_steps[1:(current_step_idx - 1)]) else character(0)

  current_step <- if (current_step_idx <= length(executed_steps))
    executed_steps[current_step_idx] else NA_character_
  cur_pos <- match(current_step, canonical)
  upstream <- if (!is.na(cur_pos) && cur_pos > 1)
    rev(canonical[seq_len(cur_pos - 1)]) else character(0)

  search_order <- unique(c(prior_steps, upstream))

  for (step_name in search_order) {
    if (!is.null(all_step_outputs[[step_name]])) {
      key <- all_step_outputs[[step_name]][1]
      name <- all_step_outputs[[step_name]][2]
      path <- step_resume_path(key, name)
      if (file.exists(path)) return(path)
    }
    if (!is.null(chained_outputs[[step_name]])) {
      key <- chained_outputs[[step_name]][1]
      name <- chained_outputs[[step_name]][2]
      path <- step_resume_path(key, name)
      if (file.exists(path)) return(path)
    }
  }
  return(NULL)
}

save_step_output <- function(step_key, default_name, chained_name = NULL, is_last = FALSE) {
  if (config$execution$mode == "step-by-step") {
    step_file <- step_output_path(step_key, default_name)
    log_milestone(sprintf("%s Saving", step_key), sprintf("Output: %s", basename(step_file)))
    saveRDS(current_object, file = step_file, compress = TRUE)
    generated_files <<- c(generated_files, basename(step_file))
    log_milestone(sprintf("%s Saved", step_key), sprintf("Output: %s", basename(step_file)))
    current_object <<- NULL
    gc()
  } else if (config$execution$mode == "chained" && !is.null(chained_name)) {
    step_file <- step_output_path(step_key, chained_name)
    log_milestone(sprintf("%s Saving", step_key), sprintf("Output: %s", basename(step_file)))
    saveRDS(current_object, file = step_file, compress = TRUE)
    generated_files <<- c(generated_files, basename(step_file))
    log_milestone(sprintf("%s Saved", step_key), sprintf("Output: %s", basename(step_file)))
  }
}

run_load_data <- function(step_num) {
  log_milestone(
    sprintf("STEP %d: Load CITE-seq Data", step_num),
    sprintf("RNA: %s", basename(config$paths$input_data$rna_path))
  )

  source("script/load_data_updated.R", local = FALSE)

  .h5_backend_needed <- function(path) {
    !is.null(path) && grepl("\\.h5(ad|mu)$", path, ignore.case = TRUE)
  }
  if (.h5_backend_needed(config$paths$input_data$rna_path) ||
      .h5_backend_needed(config$paths$input_data$adt_path %||% NULL)) {
    backend_err <- tryCatch({
      if (!requireNamespace("reticulate", quietly = TRUE))
        stop("R package 'reticulate' is not installed")
      reticulate::py_run_string("import h5py, anndata")
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(backend_err)) {
      stop(sprintf(paste0(
        "The Python backend needed to read .h5ad/.h5mu files is not working ",
        "in this environment (python: %s). Loading anyway could silently drop ",
        "the protein (ADT) data, so the run was stopped.\n %s"),
        Sys.getenv("RETICULATE_PYTHON", unset = "reticulate default"),
        backend_err), call. = FALSE)
    }
    cat(" +- Python h5ad backend OK (in-process h5py/anndata import)\n")
  }

  cat(sprintf(" +- Calling: load_data(rna_path, adt_path, project_name, species_name)\n"))
  current_object <<- load_data(
    rna_path = config$paths$input_data$rna_path,
    adt_path = config$paths$input_data$adt_path %||% NULL,
    project_name = config$load_data$project_name %||% "CellProject",
    species_name = config$load_data$species_name %||% "human",
    top_mouse_genes = config$load_data$top_mouse_genes %||% 100,
    rds_rna_key = config$load_data$rds_rna_key %||% "rna",
    rds_adt_key = config$load_data$rds_adt_key %||% "adt",
    join_layers = config$load_data$join_layers %||% TRUE,
    save_output = FALSE,
    use_disk_backed = config$processing$use_disk_backed %||% TRUE,
    bpcells_dir = file.path(config$paths$output_dir %||% tempdir(), "bpcells_cache"),

    slim_to_counts = "process_seurat" %in% (config$execution$steps %||% default_steps),
    verbose = config$execution$verbose %||% FALSE
  )
  cat(" +- Result: Seurat object loaded\n\n")

  log_data_summary(current_object, "After loading")
  log_milestone(sprintf("Step %d Complete", step_num),
                sprintf("Loaded %d cells | %d features",
                        ncol(current_object), nrow(current_object)))

  save_step_output("load_data", "load_data.rds")
}

run_add_metadata <- function(step_num) {
  if (!isTRUE(config$add_metadata$enabled)) {
    log_milestone(sprintf("STEP %d: Add Metadata", step_num), "Skipped (disabled in config)")
    return(invisible(NULL))
  }

  log_milestone(
    sprintf("STEP %d: Add Metadata", step_num),
    sprintf("Source: %s", basename(config$paths$input_data$metadata_file))
  )

  source("script/add_metadata.R", local = FALSE)

  cat(" +- Calling: add_metadata(seurat_object, metadata_file, ...)\n")
  current_object <<- add_metadata(
    seurat_object = current_object,
    metadata_file = config$paths$input_data$metadata_file,
    object_key = config$add_metadata$merge_by %||% "orig.ident",
    metadata_key = config$add_metadata$merge_by_metadata %||% "Sample",
    sheet = config$add_metadata$sheet %||% 1,
    col_types = config$add_metadata$col_types %||% "text",
    subset_to_existing = config$add_metadata$subset_to_existing %||% TRUE,
    overwrite_existing = config$add_metadata$overwrite_existing %||% FALSE,
    add_match_status = config$add_metadata$add_match_status %||% TRUE,
    backup_metadata = config$add_metadata$backup_metadata %||% FALSE,
    na_values = config$add_metadata$na_values %||% c("", "NA", "N/A", "na", "n/a"),
    verbose = config$execution$verbose %||% TRUE
  )
  cat(" +- Result: Metadata added\n\n")

  log_data_summary(current_object, "After adding metadata")
  log_milestone(sprintf("Step %d Complete", step_num), "Metadata added successfully")

  save_step_output("add_metadata", "add_metadata.rds")
}

run_process_seurat <- function(step_num) {
  log_milestone(
    sprintf("STEP %d: Process Seurat", step_num),
    sprintf("Normalization: %s | Clustering resolution: %s",
            config$processing$normalization$rna_method %||% "LogNormalize",
            config$processing$clustering$rna_resolution %||% 1.5)
  )

  source("script/process_seurat_pipeline.R", local = FALSE)

  cat(" +- Calling: process_and_sketch(seurat_object, ...)\n")
  current_object <<- process_and_sketch(
    seurat_object = current_object,
    normalization_method = config$processing$normalization$rna_method %||% "LogNormalize",
    use_sketch = config$processing$use_sketch %||% NULL,
    sketch_ncells = config$processing$sketch_ncells %||% NULL,
    sketch_method = config$processing$sketch_method %||% "LeverageScore",
    remove_ids = config$processing$remove_ids %||% NULL,
    n_variable_features = config$processing$features$n_features %||% 5000,
    pca_dims_rna = config$processing$reduction$pca_dims_rna %||% 30,
    pca_dims_adt = config$processing$reduction$pca_dims_adt %||% 10,
    clustering_resolution = config$processing$clustering$rna_resolution %||% 1.5,
    wnn_resolution = config$processing$clustering$wnn_resolution %||% 1.5,
    k_param = config$processing$clustering$k_param %||% 30,
    clustering_algorithm = config$processing$clustering$algorithm %||% 1,
    qc_params = list(
      min_features = config$processing$qc_filters$min_features %||% 200,
      max_features = config$processing$qc_filters$max_features %||% 6000,
      max_mt = config$processing$qc_filters$max_percent_mt %||% 20,
      max_rp = config$processing$qc_filters$max_percent_rp %||% 40,
      min_cells = config$processing$qc_filters$min_cells %||% 3,
      min_cells_per_sample = config$processing$qc_filters$min_cells_per_sample %||% 500
    ),
    batch_column = config$batch$column %||% "orig.ident",
    detect_batch = config$batch$detect %||% TRUE,
    correct_batch = config$batch$correct %||% "auto",
    correct_rna = config$batch$correct_rna %||% NULL,
    correct_adt = config$batch$correct_adt %||% NULL,
    harmony_theta = config$batch$harmony_theta %||% 2,
    harmony_niter = config$batch$harmony_niter %||% 10,
    lisi_threshold = config$batch$lisi_threshold %||% 0.7,
    silhouette_threshold = config$batch$silhouette_threshold %||% 0.15,
    silhouette_n_cells = config$batch$silhouette_n_cells %||% 5000L,
    wasserstein_threshold = config$batch$wasserstein_threshold %||% 0.3,
    min_cells_per_batch = config$batch$min_cells_per_batch %||% 50L,
    adt_correction_method = config$batch$adt_correction_method %||% "harmony",
    adt_kbet_threshold = config$batch$adt_kbet_threshold %||% 0.6,
    adt_guardrail_min_improvement = config$batch$adt_guardrail_min_improvement %||% 0.05,
    detect_doublets = config$processing$doublets$detect %||% TRUE,
    remove_doublets = config$processing$doublets$remove %||% TRUE,
    doublet_sample_column = config$processing$doublets$sample_column %||% NULL,
    doublet_use_adt = config$processing$doublets$use_adt %||% TRUE,
    doublet_adt_quantile = config$processing$doublets$adt_quantile %||% 0.90,
    doublet_min_cells = config$processing$doublets$min_cells %||% 50L,
    doublet_seed = config$processing$doublets$seed %||% 42L,
    use_disk_backed = config$processing$use_disk_backed %||% TRUE,
    disk_cache_dir = file.path(config$paths$output_dir %||% tempdir(), "bpcells_cache"),
    adt_keep_pattern = config$processing$adt_filters$keep_pattern %||% NULL,
    adt_drop_pattern = config$processing$adt_filters$drop_pattern %||% "isotype",
    output_dir = config$paths$output_dir %||% NULL,
    output_prefix = tools::file_path_sans_ext(config$output_files$process_seurat %||% "process_seurat"),
    verbose = config$execution$verbose %||% TRUE
  )
  cat(" +- Result: Seurat object processed, clustered, and reduced\n\n")

  dqc <- tryCatch(current_object@misc$doublet_qc, error = function(e) NULL)
  if (!is.null(dqc)) {
    log_milestone("Doublet QC", sprintf(
      "Detected %s of %s cells (%.2f%%) | RNA %s, ADT cross-lineage %s | %s",
      dqc$n_total_doublet, dqc$n_cells,
      100 * (dqc$n_total_doublet %||% 0) / max(1L, dqc$n_cells %||% 1L),
      dqc$n_rna_doublet, dqc$n_adt_flagged,
      if (isTRUE(dqc$removed)) sprintf("removed (cells now %s)", dqc$n_cells_after)
      else "flagged only"))
  }

  log_data_summary(current_object, "After processing")
  log_milestone(sprintf("Step %d Complete", step_num),
                "Seurat processing and clustering complete")

  save_step_output("process_seurat", "process_seurat.rds")
}

run_cell_type_annotation <- function(step_num) {
  log_milestone(
    sprintf("STEP %d: Cell Type Annotation", step_num),
    "SingleR + STCAT annotation methods"
  )

  source("script/cell_type_annotation.R", local = FALSE)

  cta <- config$cell_type_annotation
  cat(sprintf(" +- Calling: cell_type_annotation(method='%s', reference='%s')\n",
              cta$method %||% "both", cta$reference %||% "monaco"))
  current_object <<- cell_type_annotation(
    seurat_object = current_object,
    method = cta$method %||% "both",
    reference = cta$reference %||% "monaco",
    cluster_column = cta$cluster_column %||% "RNA_clusters",
    assay = cta$assay %||% "RNA",
    singler_label_type = cta$singler_label_type %||% c("main", "fine"),
    singler_workers = cta$singler_workers %||% config$execution$n_cores %||% 4,
    stcat_python_script = cta$stcat_python_script %||% NULL,
    stcat_batch_column = cta$stcat_batch_column %||% config$batch$column %||% "orig.ident",
    singler_fine_chunk_size = cta$singler_fine_chunk_size %||% 20000L,
    prefix = cta$prefix %||% NULL,
    fallback_column = cta$fallback_column %||% "SingleR_fine",
    verbose = config$execution$verbose %||% TRUE,
    cleanup_temp = cta$cleanup_temp %||% TRUE
  )
  cat(" +- Result: Cell type annotations added\n\n")

  log_data_summary(current_object, "After cell type annotation")
  log_milestone(sprintf("Step %d Complete", step_num), "Cell type annotation complete")

  save_step_output("cell_type_annotation", "cell_type_annotation.rds",
                   chained_name = "annotated_seurat.rds")
}

run_aggregate_ml <- function(step_num) {
  log_milestone(
    sprintf("STEP %d: Aggregate for ML", step_num),
    sprintf("Participant: %s | CellType: %s",
            config$aggregation$participant_col,
            config$aggregation$celltype_cols[1] %||% "Prediction")
  )

  source("script/aggregate_obj_ML_processing.R", local = FALSE)

  cat(" +- Calling: aggregate_seurat_for_ml_final(...)\n")

  agg_name <- config$output_files$aggregate_ml %||% "aggregate_ml.xlsx"
  agg_base <- sub("\\.(xlsx|rds|csv)$", "", agg_name, ignore.case = TRUE)
  xlsx_file <- file.path(config$paths$output_dir,
                         sprintf("%s_%s.xlsx", agg_base, run_timestamp))
  aggregate_seurat_for_ml_final(
    seurat_obj = current_object,
    participant_col = config$aggregation$participant_col,
    timepoint_col = config$aggregation$timepoint_col %||% "Timepoint",
    celltype_col = config$aggregation$celltype_cols %||% "Prediction",
    columns_to_preserve = config$aggregation$preserve_columns %||% NULL,
    assays = config$aggregation$assays %||% c("RNA", "ADT"),
    output_xlsx = xlsx_file
  )
  generated_files <<- c(generated_files, basename(xlsx_file))
  cat(sprintf(" +- Result: %s\n\n", basename(xlsx_file)))

  log_data_summary(current_object, "After aggregation")
  log_milestone(sprintf("Step %d Complete", step_num), "ML aggregation complete")

  current_object <<- NULL
  gc()
}

run_export_loupe <- function(step_num) {
  log_milestone(
    sprintf("STEP %d: Export Loupe", step_num),
    "Generating .cloupe file(s) for Loupe Browser")

  source("script/export_loupe.R", local = FALSE)

  cfg <- config$export_loupe %||% list()

  src_step <- cfg$source_step %||% "cell_type_annotation"
  src_path <- step_resume_path(src_step, paste0(src_step, ".rds"))
  if (!file.exists(src_path)) {
    log_error("Export Loupe",
              sprintf("Source object for '%s' not found: %s", src_step, src_path))
    stop(sprintf("Loupe source object not found: %s", src_path))
  }
  log_milestone("Export Loupe",
                sprintf("Loading single-cell object: %s", basename(src_path)))
  loupe_obj <- readRDS(src_path)

  base_name <- cfg$output_name %||% config$metadata$project_name %||% "loupe"
  out_name <- sprintf("%s_%s", base_name, run_timestamp)

  cat(" +- Calling: export_seurat_to_loupe(...)\n")
  created <- export_seurat_to_loupe(
    seurat_obj = loupe_obj,
    output_dir = config$paths$output_dir,
    output_name = out_name,
    assays = cfg$assays %||% c("RNA", "ADT"),
    executable_path = cfg$executable_path %||% NULL,
    verbose = config$execution$verbose %||% TRUE)

  if (length(created) > 0) {
    generated_files <<- c(generated_files, basename(created))
    cat(sprintf(" +- Result: %d Loupe file(s) created\n\n", length(created)))
  } else {
    log_warning("Export Loupe", "No .cloupe files were created (see warnings above).")
  }

  rm(loupe_obj); gc()

  current_object <<- NULL
  log_milestone(sprintf("Step %d Complete", step_num), "Loupe export complete")
}

step_dispatch <- list(
  load_data = run_load_data,
  add_metadata = run_add_metadata,
  process_seurat = run_process_seurat,
  cell_type_annotation = run_cell_type_annotation,
  aggregate_ml = run_aggregate_ml,
  export_loupe = run_export_loupe
)

valid_step_names <- names(step_dispatch)

default_steps <- c("load_data", "add_metadata", "process_seurat",
                   "cell_type_annotation", "aggregate_ml")

if (config$execution$mode %in% c("step-by-step", "chained")) {
  steps <- config$execution$steps %||% default_steps
  log_milestone("Execution Mode", config$execution$mode)
} else {
  log_error("Configuration", sprintf("Unknown mode: %s", config$execution$mode))
  finalize_log()
  quit(status = 1)
}

invalid_steps <- setdiff(steps, valid_step_names)
if (length(invalid_steps) > 0) {
  log_error("Configuration", sprintf(
    "Unknown steps: %s. Valid steps: %s",
    paste(invalid_steps, collapse = ", "),
    paste(valid_step_names, collapse = ", ")))
  finalize_log()
  quit(status = 1)
}

if ("load_data" %in% steps) {
  ld_pos <- which(steps == "load_data")
  if (ld_pos != 1) {
    log_error("Configuration", "load_data must be the first step when included")
    finalize_log()
    quit(status = 1)
  }
}

if (!is.null(args$resume)) {
  resume_step <- args$resume
  if (!resume_step %in% valid_step_names) {
    log_error("Resume", sprintf("Unknown step: %s", resume_step))
    finalize_log()
    quit(status = 1)
  }
  resume_index <- which(steps == resume_step)
  if (length(resume_index) == 0) {
    log_error("Resume", sprintf(
      "Step '%s' is not in the configured steps list (%s). Nothing to run.",
      resume_step, paste(steps, collapse = ", ")))
    finalize_log()
    quit(status = 1)
  }
  steps <- steps[resume_index:length(steps)]
  log_milestone("Resume Mode", sprintf("Resuming from: %s", resume_step))
}

log_milestone("Pipeline Steps", paste(steps, collapse = " -> "))

if (args$`dry-run`) {
  log_milestone("Dry-run Mode", "Configuration is valid")
  cat("\nAll validation checks passed!\n")
  cat("Steps: ", paste(steps, collapse = " -> "), "\n")
  cat("Mode: ", config$execution$mode, "\n\n")
  cat("To execute, run:\n")
  cat(" Rscript orchestrator.R --config ", args$config, "\n\n", sep = "")
  finalize_log()
  quit(status = 0)
}

dir.create(config$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

local({
  py <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(py) && file.exists(py)) {
    log_milestone("Python Interpreter", sprintf("RETICULATE_PYTHON -> %s", py))
    return(invisible(NULL))
  }

  env_name <- Sys.getenv("PIPELINE_PY_ENV", "")
  if (!nzchar(env_name)) {
    log_warning("Python Interpreter",
      "RETICULATE_PYTHON and PIPELINE_PY_ENV are unset; reticulate will use its default interpreter.")
    return(invisible(NULL))
  }

  conda_path <- Sys.which("conda")
  if (!nzchar(conda_path))
    conda_path <- Find(file.exists, path.expand(c(
      "/opt/conda/bin/conda",
      "~/miniforge3/bin/conda", "~/mambaforge/bin/conda",
      "~/miniconda3/bin/conda", "~/anaconda3/bin/conda")))
  resolved <- NULL
  if (requireNamespace("reticulate", quietly = TRUE) &&
      !is.null(conda_path) && nzchar(conda_path)) {
    resolved <- tryCatch(reticulate::conda_python(env_name, conda = conda_path),
                         error = function(e) NULL)
    if (is.null(resolved) || !file.exists(resolved)) {
      cand <- file.path(dirname(dirname(conda_path)), "envs", env_name, "bin", "python")
      resolved <- if (file.exists(cand)) cand else NULL
    }
  }
  if (!is.null(resolved) && file.exists(resolved)) {
    Sys.setenv(RETICULATE_PYTHON = resolved)
    log_milestone("Python Interpreter", sprintf("Resolved '%s' -> %s", env_name, resolved))
  } else {
    log_warning("Python Interpreter",
      sprintf("Could not resolve conda env '%s'; reticulate will use its default.", env_name))
  }
})

for (step_idx in seq_along(steps)) {
  step_name <- steps[step_idx]
  step_fn <- step_dispatch[[step_name]]

  if (is.null(current_object) && !step_name %in% c("load_data", "export_loupe")) {
    prior_file <- load_previous_output(step_idx, steps)
    if (is.null(prior_file) || !file.exists(prior_file)) {
      log_error(step_name, "No prior step output found to resume from")
      finalize_log()
      quit(status = 1)
    }
    current_object <- readRDS(prior_file)
    log_milestone("Loaded prior output", sprintf("File: %s", basename(prior_file)))
  }

  tryCatch({
    step_fn(step_idx)
  }, error = function(e) {
    log_error(sprintf("Step %d: %s", step_idx, step_name), as.character(e))
    finalize_log()
    quit(status = 1)
  })
}

log_milestone("PIPELINE EXECUTION COMPLETE", "All steps executed successfully")
finalize_log()

cat("\n+================================================================+\n")
cat("| PIPELINE COMPLETE |\n")
cat("+================================================================+\n")
cat("Results saved to: ", config$paths$output_dir, "\n\n", sep = "")

if (length(generated_files) > 0) {
  cat("Output files generated this run:\n")
  for (f in generated_files) cat(" - ", f, "\n", sep = "")
}

cat("\n")
quit(status = 0)
