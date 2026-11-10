`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(optparse)
})

if (!interactive() && sys.nframe() == 0) {

  option_list <- list(
    make_option(c("-c", "--config"),
                type = "character",
                default = "process_and_sketch_config.yml",
                help = "Path to YAML configuration file")
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  if (!file.exists(opt$config))
    stop(sprintf("Configuration file not found: %s", opt$config))

  config <- yaml::read_yaml(opt$config)

  input_path <- file.path(config$paths$input_dir, config$paths$input_filename)
  output_path <- file.path(config$paths$output_dir, config$paths$output_filename)

  processed_obj <- process_and_sketch(
    input_path = input_path,
    output_path = output_path,
    normalization_method = config$processing$normalization_method,
    use_sketch = config$processing$use_sketch,
    sketch_ncells = config$processing$sketch_ncells,
    remove_ids = config$remove_ids,
    qc_params = config$qc_params,
    n_variable_features = config$analysis$n_variable_features,
    pca_dims_rna = config$analysis$pca_dims_rna,
    pca_dims_adt = config$analysis$pca_dims_adt,
    clustering_resolution = config$analysis$clustering_resolution,
    wnn_resolution = config$analysis$wnn_resolution,
    k_param = config$analysis$k_param,

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
    verbose = config$verbose %||% TRUE
  )

  cat("Processing complete!\n")
  quit(save = "no")
}

#' Process and Sketch Seurat Object with Multi-modal Integration
#'
#' @description
#' Production function for processing single-cell data with optional protein
#' (ADT) data. Supports Seurat v4/v5, multiple normalization methods, and
#' automatic sketching for large datasets. Batch effect detection and
#' correction are embedded at the correct points in the pipeline so that UMAP,
#' clustering, and WNN are computed only once, from the (optionally corrected)
#' embeddings.
#'
#' Pipeline order
#' --------------
#' 1. Load / validate
#' 2. QC filter
#' 2b. Doublet detection + removal (scDblFinder per capture + ADT cross-lineage)
#' 2c. Sample & feature cleanup (drop small batches, clean ADT panel)
#' 3. Normalize (RNA LogNormalize/SCT + ADT CLR)
#' 4. ADT batch detection (Wasserstein per protein on CLR data)
#' 4a. ADT batch correction (ADTnorm -> re-CLR) - if recommended / forced
#' 5. RNA PCA (full or sketch projection)
#' 5a. RNA batch detection (silhouette on PCA)
#' 5b. RNA batch correction (Harmony) - if recommended / forced
#' 5c. Post-correction validation - re-runs detection on corrected embedding
#' to verify improvement; results stored in batch_qc$rna_post_correction
#' 6. Embeddings, clustering, WNN (uses corrected reductions throughout)
#' 7. Store batch_qc in misc, save, return
#'
#' @param seurat_object Seurat object (or NULL when input_path is used).
#' @param input_path Path to RDS input file (default: NULL).
#' @param output_path Path to save output RDS (default: NULL - no save).
#' @param normalization_method "LogNormalize" or "SCT" (default: "LogNormalize").
#' @param use_sketch Logical - force sketch mode (default: NULL, auto >= 850k cells).
#' @param sketch_ncells Number of cells for sketch (default: NULL -> 40 % of total).
#' @param sketch_method Sketching algorithm (default: "LeverageScore").
#' @param remove_ids Character vector of orig.ident values to remove (default: NULL).
#' @param qc_params Named list: min_features, max_features, max_mt, max_rp.
#' @param n_variable_features Number of variable features for RNA (default: 5000).
#' @param pca_dims_rna Number of RNA PCA dims for downstream steps (default: 30).
#' @param pca_dims_adt Number of ADT PCA dims for downstream steps (default: 20).
#' @param clustering_resolution Resolution for RNA and ADT clusters (default: 2).
#' @param wnn_resolution Resolution for WNN clusters (default: 2).
#' @param k_param k for nearest-neighbour graphs (default: 30).
#' @param clustering_algorithm Seurat clustering algorithm (default: 1 = Louvain).
#'
#' @param detect_doublets Run per-capture doublet detection (default: TRUE).
#' @param remove_doublets Remove detected doublets, as for other QC filters
#' (default: TRUE). When FALSE the cells are flagged and retained.
#' @param doublet_sample_column Metadata column identifying the capture/library.
#' Doublets form within one droplet emulsion, so detection is run separately
#' per capture. Defaults to NULL -> falls back to "orig.ident".
#' @param doublet_use_adt Use the ADT cross-lineage signature as a confirmatory
#' signal on CITE-seq data (default: TRUE).
#' @param doublet_adt_quantile Per-marker upper-tail quantile defining ADT
#' positivity for the cross-lineage check (default: 0.90).
#' @param doublet_min_cells Minimum cells for a capture to be classified
#' (default: 50).
#' @param doublet_seed Random seed for scDblFinder (default: 42).
#' @param batch_column Metadata column holding batch labels (default: "orig.ident").
#' @param detect_batch Run batch-effect detection? (default: TRUE).
#' @param correct_batch "auto" corrects when detection flags it; TRUE forces
#' correction; FALSE disables it entirely (default: "auto").
#' @param harmony_theta Harmony diversity penalty theta (default: 2).
#' @param harmony_niter Maximum Harmony iterations (default: 10).
#' @param lisi_threshold Retained for backward compatibility but UNUSED - LISI was
#' decoupled; silhouette ASW-by-batch is the sole RNA batch metric (default: 0.7).
#' @param silhouette_threshold ASW-by-batch above which the RNA batch flag is
#' set (default: 0.15).
#' @param silhouette_n_cells Maximum number of cells to subsample for silhouette
#' computation (default: 5000). Higher values are more accurate but use O(n^2)
#' memory for the distance matrix.
#' @param adt_correction_method HOW to correct ADT batch effects when detection
#' (or forcing) says correction should run: "harmony" (default - Harmony on
#' the apca embedding consumed by WNN; marker values untouched), "adtnorm"
#' (legacy value-level correction of the ADT data layer), or "none".
#' @param adt_kbet_threshold kBET-style rejection rate on the ADT embedding
#' above which the joint-embedding detector flags a batch effect
#' (default: 0.6; 0 = perfectly mixed, 1 = fully separated).
#' @param adt_guardrail_min_improvement Minimum kBET rejection improvement the
#' correction must achieve to be adopted; otherwise it is reverted and the
#' pipeline continues uncorrected (default: 0.05).
#' @param wasserstein_threshold Per-protein mean Wasserstein distance above which
#' a protein is flagged (default: 0.3). ADT correction is recommended when
#' > 10 % of proteins are flagged.
#' @param min_cells_per_batch Minimum cells per batch for reliable detection
#' (default: 50). A warning is issued if any batch is smaller.
#'
#' @param use_disk_backed Logical - move raw counts to on-disk BPCells
#' matrices before QC so subset + normalization stay within RAM
#' (default: TRUE; silently falls back to in-memory if BPCells absent).
#' @param disk_cache_dir Directory for the on-disk BPCells matrices
#' (default: NULL -> a subdir of tempdir()). Must persist for the life
#' of the object; pass the run output dir for resume safety.
#' @param qc_params Adds \code{min_cells}: drop genes detected in fewer
#' than this many cells before normalization (default: 3).
#' @param verbose Print progress messages (default: TRUE).
#'
#' @return Processed Seurat object. \code{seurat_object@misc$doublet_qc} holds
#' the per-capture doublet summary (cells, doublets, percentage, whether the
#' capture was skipped) and whether removal was applied.
#' \code{seurat_object@misc$batch_qc} holds
#' all detection results and flags.
#'
#' @export
process_and_sketch <- function(
    seurat_object = NULL,
    input_path = NULL,
    output_path = NULL,
    normalization_method = "LogNormalize",
    use_sketch = NULL,
    sketch_ncells = NULL,
    sketch_method = "LeverageScore",
    remove_ids = NULL,
    qc_params = list(
      min_features = 200,
      max_features = 6000,
      max_mt = 20,
      max_rp = 40,
      min_cells = 3,
      min_cells_per_sample = 500
    ),
    n_variable_features = 5000,
    pca_dims_rna = 30,
    pca_dims_adt = 20,
    clustering_resolution = 2,
    wnn_resolution = 2,
    k_param = 30,
    clustering_algorithm = 1,

    batch_column = "orig.ident",
    detect_batch = TRUE,
    correct_batch = "auto",
    correct_rna = NULL,
    correct_adt = NULL,
    harmony_theta = 2,
    harmony_niter = 10,
    lisi_threshold = 0.7,
    silhouette_threshold = 0.15,
    silhouette_n_cells = 5000L,
    wasserstein_threshold = 0.3,
    min_cells_per_batch = 50L,
    adt_correction_method = "harmony",
    adt_kbet_threshold = 0.6,
    adt_guardrail_min_improvement = 0.05,

    detect_doublets = TRUE,
    remove_doublets = TRUE,
    doublet_sample_column = NULL,
    doublet_use_adt = TRUE,
    doublet_adt_quantile = 0.90,
    doublet_min_cells = 50L,
    doublet_seed = 42L,
    use_disk_backed = TRUE,
    disk_cache_dir = NULL,
    adt_keep_pattern = NULL,
    adt_drop_pattern = "isotype",
    output_dir = NULL,
    output_prefix = NULL,
    verbose = TRUE
) {

  options(future.globals.maxSize = 8 * 1024^3)
  if (requireNamespace("future", quietly = TRUE)) {
    try(future::plan("sequential"), silent = TRUE)
  }

  if (verbose) cat("=== Starting Seurat Processing Pipeline ===\n")

  if (!is.null(input_path)) {
    if (!file.exists(input_path))
      stop(sprintf("Input file not found: %s", input_path))
    if (verbose) {
      cat(sprintf("\n--- Loading Data ---\nReading from: %s\n", input_path))
    }
    seurat_object <- readRDS(input_path)
    if (verbose) {
      cat(sprintf("Initial cell count: %d\n", ncol(seurat_object)))
      cat(sprintf("Initial feature count: %d\n", nrow(seurat_object)))
    }
  }

  if (is.null(seurat_object))
    stop("Either 'seurat_object' or 'input_path' must be provided")

  if (!normalization_method %in% c("LogNormalize", "SCT"))
    stop("normalization_method must be 'LogNormalize' or 'SCT'")

  has_adt <- "ADT" %in% names(seurat_object@assays)
  if (verbose)
    cat(sprintf("ADT data detected: %s\n", ifelse(has_adt, "YES", "NO")))

  n_cells <- ncol(seurat_object)
  if (is.null(use_sketch)) {
    use_sketch <- n_cells > 850000
    if (verbose)
      cat(sprintf(
        "Sketch mode: %s (auto; threshold 850,000 cells, current %d)\n",
        ifelse(use_sketch, "ENABLED", "DISABLED"), n_cells))
  } else {
    if (verbose)
      cat(sprintf("Sketch mode: %s (explicitly set)\n",
                  ifelse(use_sketch, "ENABLED", "DISABLED")))
  }

  if (use_sketch && is.null(sketch_ncells)) {
    sketch_ncells <- round(n_cells * 0.4)
    if (verbose)
      cat(sprintf("Sketch size set to 40%% of cells: %d\n", sketch_ncells))
  }

  if (!is.null(remove_ids) && length(remove_ids) > 0) {
    if ("orig.ident" %in% colnames(seurat_object@meta.data)) {
      n_before <- ncol(seurat_object)
      keep_cells <- colnames(seurat_object)[
        !seurat_object$orig.ident %in% remove_ids]
      if (length(keep_cells) < 2)
        stop(sprintf(paste0(
          "remove_ids would remove all (or all but %d) cells. Requested: %s. ",
          "Available orig.ident values: %s"),
          length(keep_cells), paste(remove_ids, collapse = ", "),
          paste(unique(seurat_object$orig.ident), collapse = ", ")))
      seurat_object <- subset(seurat_object, cells = keep_cells)
      n_after <- ncol(seurat_object)
      if (verbose && n_before != n_after)
        cat(sprintf(
          "\n--- Sample Removal ---\nRemoved %d cells from: %s\nCells remaining: %d\n",
          n_before - n_after,
          paste(remove_ids, collapse = ", "),
          n_after))
    } else if (verbose) {
      cat("\nWarning: 'orig.ident' not found in metadata. Skipping sample removal.\n")
    }
  }

  if (isTRUE(use_disk_backed)) {
    if (verbose) cat("\n--- Converting counts to disk-backed (BPCells) ---\n")
    seurat_object <- .to_disk_backed(
      seurat_object, cache_dir = disk_cache_dir, verbose = verbose)
  }

  seurat_object <- tryCatch({
    keep_reducs <- tryCatch(SeuratObject::Reductions(seurat_object),
                             error = function(e) NULL)

    count_layers <- unique(unlist(lapply(
      names(seurat_object@assays),
      function(a) tryCatch(
        SeuratObject::Layers(seurat_object, assay = a, search = "counts"),
        error = function(e) character(0)))))
    if (length(count_layers) == 0) count_layers <- "counts"
    dieted <- DietSeurat(seurat_object, layers = count_layers, dimreducs = keep_reducs)
    if (verbose)
      cat(sprintf("\n--- Slimming object (DietSeurat: keeping layers %s) ---\n",
                  paste(count_layers, collapse = ", ")))
    dieted
  }, error = function(e) {
    if (verbose) cat(sprintf("\nDietSeurat skipped (%s)\n", conditionMessage(e)))
    seurat_object
  })
  gc()

  if (verbose) cat("\n--- Quality Control ---\n")

  seurat_object <- PercentageFeatureSet(
    seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object <- PercentageFeatureSet(
    seurat_object, pattern = "^RP[SL]", col.name = "percent.rp")

  n_before <- ncol(seurat_object)
  seurat_object <- subset(seurat_object,
    subset = nFeature_RNA > qc_params$min_features &
             nFeature_RNA < qc_params$max_features &
             percent.mt < qc_params$max_mt &
             percent.rp < qc_params$max_rp)
  n_after <- ncol(seurat_object)

  if (n_after < 2)
    stop(sprintf(paste0(
      "QC filtering retained %d cell(s) of %d. Thresholds are too strict: ",
      "min_features=%s, max_features=%s, max_mt=%s, max_rp=%s"),
      n_after, n_before, qc_params$min_features, qc_params$max_features,
      qc_params$max_mt, qc_params$max_rp))

  if (verbose)
    cat(sprintf("Cells before QC: %d\nCells after QC: %d (%.1f%% retained)\n",
                n_before, n_after, 100 * n_after / max(1L, n_before)))
  gc()

  seurat_object <- .filter_genes_by_min_cells(
    seurat_object, min_cells = qc_params$min_cells %||% 3, verbose = verbose)
  gc()

  doublet_qc <- NULL
  if (isTRUE(detect_doublets)) {
    if (verbose) cat("\n--- Doublet Detection ---\n")
    seurat_object <- tryCatch(
      .detect_doublets(
        seurat_object,
        sample_column = doublet_sample_column %||% batch_column,
        assay = "RNA",
        adt_assay = "ADT",
        use_adt_signal = isTRUE(doublet_use_adt) && has_adt,
        adt_positive_quantile = doublet_adt_quantile,
        min_cells_per_sample = doublet_min_cells,
        random_seed = doublet_seed,
        verbose = verbose),
      error = function(e) {
        if (verbose)
          cat(sprintf("Doublet detection failed (%s) - continuing without it.\n",
                      conditionMessage(e)))
        seurat_object
      })

    doublet_qc <- seurat_object@misc$doublet_qc

    if (isTRUE(remove_doublets)) {
      seurat_object <- .remove_doublets(seurat_object, verbose = verbose)
    } else if (verbose) {
      cat(" remove_doublets = FALSE - doublets flagged but retained.\n")
    }

    if (!is.null(doublet_qc)) {
      doublet_qc$removed <- isTRUE(remove_doublets)
      doublet_qc$n_cells_after <- ncol(seurat_object)

      dbl_dir <- output_dir %||% (if (!is.null(output_path)) dirname(output_path) else NULL)
      if (!is.null(dbl_dir) && nzchar(dbl_dir) &&
          is.data.frame(doublet_qc$per_sample) && nrow(doublet_qc$per_sample) > 0) {
        tryCatch({
          dir.create(dbl_dir, showWarnings = FALSE, recursive = TRUE)
          pfx <- output_prefix %||% "process_seurat"
          csv <- file.path(dbl_dir, sprintf("%s_doublet_per_sample.csv", pfx))
          utils::write.csv(doublet_qc$per_sample, csv, row.names = FALSE)
          if (verbose) cat(sprintf(" Per-capture doublet summary: %s\n", csv))
        }, error = function(e)
          if (verbose) cat(sprintf(" Could not write doublet CSV (%s)\n",
                                   conditionMessage(e))))
      }
    }
    if (!is.null(doublet_qc)) seurat_object@misc$doublet_qc <- doublet_qc
    gc()
  } else if (verbose) {
    cat("\n--- Doublet Detection: SKIPPED (detect_doublets = FALSE) ---\n")
  }

  seurat_object <- .drop_small_batches(
    seurat_object, batch_column = batch_column,
    min_cells_per_sample = qc_params$min_cells_per_sample %||% 500, verbose = verbose)

  if (has_adt)
    seurat_object <- .clean_adt_features(
      seurat_object, keep_pattern = adt_keep_pattern,
      drop_pattern = adt_drop_pattern, batch_column = batch_column,
      verbose = verbose)
  gc()

  if (verbose) cat(sprintf("\n--- Normalization (%s) ---\n", normalization_method))

  if (normalization_method == "LogNormalize") {
    seurat_object <- NormalizeData(
      seurat_object,
      normalization.method = "LogNormalize",
      verbose = FALSE)
    if (has_adt)
      seurat_object <- .normalize_adt(seurat_object, verbose = verbose)
  } else {
    seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", verbose = FALSE)
    if (has_adt)
      seurat_object <- .normalize_adt(seurat_object, verbose = verbose)
  }
  gc()

  batch_qc <- list(
    batch_column = batch_column,
    n_batches = 1L,
    detection_skipped = TRUE,
    rna = NULL,
    adt = NULL,
    rna_correction_applied = FALSE,
    adt_correction_applied = FALSE,
    overall_correction_recommended = FALSE
  )

  .resolve_correct <- function(x) {
    if (is.logical(x)) return(x[1])
    xs <- tolower(trimws(as.character(x)[1]))
    if (identical(xs, "auto")) NA
    else if (xs %in% c("true", "yes", "1", "t")) TRUE
    else if (xs %in% c("false", "no", "0", "f")) FALSE
    else FALSE
  }
  run_correct_rna <- .resolve_correct(correct_rna %||% correct_batch)
  run_correct_adt <- .resolve_correct(correct_adt %||% correct_batch)

  adt_correction_method <- match.arg(
    tolower(trimws(as.character(adt_correction_method)[1])),
    c("harmony", "adtnorm", "none"))

  run_detect <- isTRUE(detect_batch)

  if (!batch_column %in% colnames(seurat_object@meta.data)) {
    if (verbose)
      cat(sprintf(
        "\nBatch column '%s' not found in metadata - skipping batch detection.\n",
        batch_column))
    run_detect <- FALSE
    if (is.na(run_correct_rna)) run_correct_rna <- FALSE
    if (is.na(run_correct_adt)) run_correct_adt <- FALSE
  } else {
    batch_labels <- seurat_object@meta.data[[batch_column]]
    n_batches <- length(unique(na.omit(batch_labels)))
    batch_qc$n_batches <- n_batches

    if (n_batches < 2L) {
      if (verbose)
        cat(sprintf(
          "\nOnly %d unique batch value found - skipping batch detection.\n",
          n_batches))
      run_detect <- FALSE
      if (is.na(run_correct_rna)) run_correct_rna <- FALSE
      if (is.na(run_correct_adt)) run_correct_adt <- FALSE
    } else {
      min_batch_n <- min(table(batch_labels))
      if (min_batch_n < min_cells_per_batch && verbose)
        cat(sprintf(
          "\nWarning: smallest batch has %d cells (< min_cells_per_batch = %d). ",
          min_batch_n, min_cells_per_batch),
          "Detection metrics may be unreliable.\n")
    }
  }

  adt_qc_dir <- if (!is.null(output_dir) && nzchar(output_dir)) output_dir
                else if (!is.null(output_path) && nzchar(output_path)) dirname(output_path)
                else NULL
  .adt_base_prefix <- if (!is.null(output_prefix) && nzchar(output_prefix)) output_prefix
                      else if (!is.null(output_path) && nzchar(output_path))
                        tools::file_path_sans_ext(basename(output_path))
                      else "process_seurat"
  adt_qc_prefix <- paste0(.adt_base_prefix, "_adt")

  if (run_detect && has_adt) {
    if (verbose) cat("\n--- ADT Batch Effect Detection ---\n")
    adt_qc <- .detect_adt_batch(
      seurat_object = seurat_object,
      batch_column = batch_column,
      n_batches = batch_qc$n_batches,
      wasserstein_threshold = wasserstein_threshold,
      verbose = verbose,
      save_dir = adt_qc_dir,
      save_prefix = adt_qc_prefix,
      kbet_threshold = adt_kbet_threshold)
    batch_qc$adt <- adt_qc
    batch_qc$detection_skipped <- FALSE
  }

  if (isTRUE(run_correct_adt) && has_adt && is.null(batch_qc$adt)) {
    if (verbose) cat("\n--- ADT Batch Detection (forced before correction) ---\n")
    adt_qc <- .detect_adt_batch(
      seurat_object = seurat_object,
      batch_column = batch_column,
      n_batches = batch_qc$n_batches,
      wasserstein_threshold = wasserstein_threshold,
      verbose = verbose,
      save_dir = adt_qc_dir,
      save_prefix = adt_qc_prefix,
      kbet_threshold = adt_kbet_threshold)
    batch_qc$adt <- adt_qc
    batch_qc$detection_skipped <- FALSE
  }

  should_correct_adt <-
    (isTRUE(run_correct_adt) ||
     (isTRUE(is.na(run_correct_adt)) && !is.null(batch_qc$adt) && isTRUE(batch_qc$adt$correction_recommended))) &&
    adt_correction_method != "none"

  if (has_adt && adt_correction_method == "none" &&
      (isTRUE(run_correct_adt) ||
       (!is.null(batch_qc$adt) && isTRUE(batch_qc$adt$correction_recommended)))) {
    if (verbose) cat("\n ADT correction requested/recommended but adt_correction_method = 'none' - skipping.\n")
  }

  batch_qc$adt_correction_method <- adt_correction_method

  pending_adt_harmony <- should_correct_adt && has_adt &&
                         adt_correction_method == "harmony"
  if (pending_adt_harmony && verbose)
    cat("\n ADT correction deferred to embedding stage (Harmony on apca).\n")

  if (should_correct_adt && has_adt && adt_correction_method == "adtnorm") {
    if (verbose) cat("\n--- ADT Batch Correction (ADTnorm) ---\n")

    labels_gr <- seurat_object@meta.data[[batch_column]]
    qc_before <- .adt_qc_pca_scores(seurat_object)
    kbet_before <- if (!is.null(qc_before))
      .kbet_rejection(qc_before$emb, labels_gr[qc_before$idx]) else NA_real_

    adt_assay_backup <- seurat_object[["ADT"]]

    adt_corrected <- tryCatch({
      seurat_object <- .correct_adt_batch(
        seurat_object = seurat_object,
        batch_column = batch_column,
        verbose = verbose)
      TRUE
    }, error = function(e) {
      if (verbose) cat(sprintf(
        paste0(" ADT batch correction FAILED (%s).\n",
               " Continuing with CLR-normalized ADT (uncorrected).\n"),
        conditionMessage(e)))
      FALSE
    })

    if (isTRUE(adt_corrected)) {
      qc_after <- .adt_qc_pca_scores(seurat_object)
      kbet_after <- if (!is.null(qc_after))
        .kbet_rejection(qc_after$emb, labels_gr[qc_after$idx]) else NA_real_
      adopt <- .adt_guardrail_adopt(kbet_before, kbet_after,
                                    adt_guardrail_min_improvement)
      batch_qc$adt_guardrail <- list(
        method = "adtnorm",
        kbet_before = kbet_before,
        kbet_after = kbet_after,
        min_improvement = adt_guardrail_min_improvement,
        adopted = adopt)
      if (!adopt) {
        seurat_object[["ADT"]] <- adt_assay_backup
        adt_corrected <- FALSE
        if (verbose) cat(sprintf(
          paste0(" GUARDRAIL: ADTnorm did not improve batch mixing ",
                 "(kBET %.3f -> %s) - REVERTED to CLR data.\n"),
          kbet_before,
          if (is.finite(kbet_after)) sprintf("%.3f", kbet_after) else "NA"))
      } else if (verbose) {
        cat(sprintf(
          " Guardrail: batch mixing improved (kBET %.3f -> %.3f) - ADTnorm adopted.\n",
          kbet_before, kbet_after))
      }
    }
    rm(adt_assay_backup); gc(FALSE)

    batch_qc$adt_correction_applied <- isTRUE(adt_corrected)
    batch_qc$adt_correction_failed <- !isTRUE(adt_corrected)
    if (isTRUE(adt_corrected) && verbose) cat(" ADT batch correction applied.\n")
  }

  assay_name <- ifelse(normalization_method == "SCT", "SCT", "RNA")

  if (use_sketch) {
    if (verbose) cat("\n--- RNA Analysis (Sketch Mode) ---\n")
    seurat_object <- .compute_rna_pca_sketch(
      seurat_object = seurat_object,
      assay_name = assay_name,
      n_variable_features = n_variable_features,
      sketch_ncells = sketch_ncells,
      sketch_method = sketch_method,
      verbose = verbose)
    rna_pca_name <- "full.pca"
  } else {
    if (verbose) cat("\n--- RNA Analysis (Full Mode) ---\n")
    seurat_object <- .compute_rna_pca_full(
      seurat_object = seurat_object,
      assay_name = assay_name,
      n_variable_features = n_variable_features,
      verbose = verbose)
    rna_pca_name <- "pca"
  }
  gc()

  if (run_detect && batch_qc$n_batches >= 2L) {
    if (verbose) cat("\n--- RNA Batch Effect Detection ---\n")
    rna_qc <- .detect_rna_batch(
      seurat_object = seurat_object,
      rna_reduction = rna_pca_name,
      batch_column = batch_column,
      n_batches = batch_qc$n_batches,
      lisi_threshold = lisi_threshold,
      silhouette_threshold = silhouette_threshold,
      silhouette_n_cells = silhouette_n_cells,
      verbose = verbose)
    batch_qc$rna <- rna_qc
    batch_qc$detection_skipped <- FALSE
  }

  should_correct_rna <-
    isTRUE(run_correct_rna) ||
    (isTRUE(is.na(run_correct_rna)) && !is.null(batch_qc$rna) && isTRUE(batch_qc$rna$correction_recommended))

  rna_reduction_final <- rna_pca_name

  if (should_correct_rna) {
    if (verbose) cat("\n--- RNA Batch Correction (Harmony) ---\n")
    if (use_sketch) {

      corr <- tryCatch({
        obj <- .correct_rna_batch_sketch(
          seurat_object = seurat_object, assay_name = assay_name,
          batch_column = batch_column, theta = harmony_theta,
          niter = harmony_niter, pca_dims = pca_dims_rna, verbose = verbose)
        list(obj = obj, red = "harmony.full")
      }, error = function(e) {
        if (verbose) cat(sprintf(
          " Sketch Harmony/projection failed (%s); falling back to Harmony on full PCA.\n",
          conditionMessage(e)))
        obj <- .correct_rna_batch(
          seurat_object = seurat_object, rna_reduction = rna_pca_name,
          batch_column = batch_column, theta = harmony_theta,
          niter = harmony_niter, verbose = verbose)
        list(obj = obj, red = "harmony")
      })
      seurat_object <- corr$obj
      rna_reduction_final <- corr$red
    } else {
      seurat_object <- .correct_rna_batch(
        seurat_object = seurat_object,
        rna_reduction = rna_pca_name,
        batch_column = batch_column,
        theta = harmony_theta,
        niter = harmony_niter,
        verbose = verbose)
      rna_reduction_final <- "harmony"
    }
    batch_qc$rna_correction_applied <- TRUE
    if (verbose) cat(sprintf(" RNA batch correction applied (reduction: %s).\n",
                             rna_reduction_final))
  }

  batch_qc$overall_correction_recommended <-
    (!is.null(batch_qc$rna) && isTRUE(batch_qc$rna$correction_recommended)) ||
    (!is.null(batch_qc$adt) && isTRUE(batch_qc$adt$correction_recommended))

  if (isTRUE(batch_qc$rna_correction_applied) && batch_qc$n_batches >= 2L) {
    if (verbose) cat("\n--- Post-Correction Validation (RNA) ---\n")
    rna_qc_post <- .detect_rna_batch(
      seurat_object = seurat_object,
      rna_reduction = rna_reduction_final,
      batch_column = batch_column,
      n_batches = batch_qc$n_batches,
      lisi_threshold = lisi_threshold,
      silhouette_threshold = silhouette_threshold,
      silhouette_n_cells = silhouette_n_cells,
      verbose = verbose)
    batch_qc$rna_post_correction <- rna_qc_post
    if (verbose) {
      improved <- !isTRUE(rna_qc_post$correction_recommended)
      cat(sprintf(" Post-correction status: %s\n",
                  ifelse(improved, "IMPROVED (no longer flagged)",
                         "STILL FLAGGED (correction may be insufficient)")))
    }
  }

  if (verbose) cat("\n--- Computing Embeddings, Clusters, and WNN ---\n")
  seurat_object <- .compute_embeddings(
    seurat_object = seurat_object,
    has_adt = has_adt,
    rna_reduction = rna_reduction_final,
    assay_name = assay_name,
    pca_dims_rna = pca_dims_rna,
    pca_dims_adt = pca_dims_adt,
    clustering_resolution = clustering_resolution,
    wnn_resolution = wnn_resolution,
    k_param = k_param,
    clustering_algorithm = clustering_algorithm,
    verbose = verbose,
    adt_harmony = pending_adt_harmony,
    batch_column = batch_column,
    harmony_theta = harmony_theta,
    harmony_niter = harmony_niter,
    adt_guardrail_min_improvement = adt_guardrail_min_improvement)
  gc()

  if (!is.null(seurat_object@misc$adt_harmony_qc)) {
    batch_qc$adt_guardrail <- seurat_object@misc$adt_harmony_qc
    batch_qc$adt_correction_applied <- isTRUE(batch_qc$adt_guardrail$adopted)
    batch_qc$adt_correction_failed <- isTRUE(batch_qc$adt_guardrail$failed)
    seurat_object@misc$adt_harmony_qc <- NULL
  }
  seurat_object@misc$batch_qc <- batch_qc
  if (!is.null(doublet_qc)) seurat_object@misc$doublet_qc <- doublet_qc

  if (verbose) {
    cat("\n=== Pipeline Complete ===\n")
    cat(sprintf("Final cell count: %d\n", ncol(seurat_object)))
    cat(sprintf("Reductions available: %s\n",
                paste(names(seurat_object@reductions), collapse = ", ")))
    cat(sprintf("RNA reduction used for embeddings: %s\n", rna_reduction_final))
    if (has_adt) cat("Multi-modal analysis completed\n")
    if (batch_qc$rna_correction_applied || batch_qc$adt_correction_applied)
      cat(sprintf("Batch correction applied - RNA: %s | ADT: %s\n",
                  batch_qc$rna_correction_applied,
                  batch_qc$adt_correction_applied))
  }

  if (!is.null(output_path)) {
    out_dir <- dirname(output_path)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    if (verbose) cat(sprintf("\n--- Saving Results ---\nWriting to: %s\n", output_path))
    saveRDS(seurat_object, file = output_path)
    if (verbose)
      cat(sprintf("File saved (%.2f GB)\n", file.info(output_path)$size / 1024^3))
  }

  return(seurat_object)
}

.filter_genes_by_min_cells <- function(seurat_object, min_cells = 3, verbose = TRUE) {
  if (is.null(min_cells) || min_cells <= 0 ||
      !"RNA" %in% names(seurat_object@assays))
    return(seurat_object)
  tryCatch({
    rna_counts <- SeuratObject::LayerData(seurat_object, assay = "RNA", layer = "counts")
    cells_per_gene <- rowSums(rna_counts > 0)
    keep_genes <- names(cells_per_gene)[cells_per_gene >= min_cells]
    n_genes_before <- length(cells_per_gene)
    if (length(keep_genes) > 0 && length(keep_genes) < n_genes_before) {
      suppressWarnings(seurat_object[["RNA"]] <- subset(seurat_object[["RNA"]], features = keep_genes))
      if (verbose)
        cat(sprintf("Genes before filter: %d\nGenes after filter: %d (min_cells = %d)\n",
                    n_genes_before, length(keep_genes), min_cells))
    }
    seurat_object
  }, error = function(e) {
    if (verbose) cat(sprintf("Gene filter skipped (%s)\n", conditionMessage(e)))
    seurat_object
  })
}

.drop_small_batches <- function(seurat_object, batch_column,
                                min_cells_per_sample = 500, verbose = TRUE) {
  if (is.null(min_cells_per_sample) || min_cells_per_sample <= 0 ||
      !batch_column %in% colnames(seurat_object@meta.data))
    return(seurat_object)
  bt <- table(seurat_object@meta.data[[batch_column]])
  small <- names(bt)[bt < min_cells_per_sample]
  if (length(small) == 0) return(seurat_object)
  keep_cells <- colnames(seurat_object)[
    !seurat_object@meta.data[[batch_column]] %in% small]
  if (length(keep_cells) >= 1 && length(keep_cells) < ncol(seurat_object)) {
    if (verbose)
      cat(sprintf(
        "\n--- Dropping small %s batches (< %d cells) ---\nRemoved %d batch(es): %s\nCells: %d -> %d\n",
        batch_column, min_cells_per_sample, length(small),
        paste(small, collapse = ", "), ncol(seurat_object), length(keep_cells)))
    seurat_object <- subset(seurat_object, cells = keep_cells)
  } else if (verbose) {
    cat(sprintf(
      "\nSmall-batch removal skipped: would drop all cells (threshold %d).\n",
      min_cells_per_sample))
  }
  seurat_object
}

.adt_lineage_groups <- function() {
  list(
    T_cell = c("CD3"),
    B_cell = c("CD19", "CD20"),
    NK_cell = c("CD56"),
    Myeloid = c("CD14", "CD33"),
    Platelet = c("CD41", "CD42B", "CD61")
  )
}

.adt_match_features <- function(features, wanted) {

  tokenise <- function(x) {
    tk <- strsplit(toupper(gsub("[^A-Za-z0-9]+", "_", x)), "_", fixed = TRUE)
    lapply(tk, function(z) z[nzchar(z)])
  }
  tok_list <- tokenise(features)
  hits <- character(0)
  for (w_tok in tokenise(wanted)) {
    if (!length(w_tok)) next
    w_glued <- paste(w_tok, collapse = "")
    n <- length(w_tok)
    idx <- which(vapply(tok_list, function(tk) {
      if (w_glued %in% tk) return(TRUE)
      if (n == 1L || length(tk) < n) return(FALSE)
      any(vapply(seq_len(length(tk) - n + 1L),
                 function(i) identical(tk[i:(i + n - 1L)], w_tok), logical(1)))
    }, logical(1)))
    hits <- c(hits, features[idx])
  }
  unique(hits)
}

.adt_crosslineage_doublets <- function(seurat_object, assay = "ADT",
                                       positive_quantile = 0.90,
                                       min_groups = 2L, verbose = TRUE) {
  n_cells <- ncol(seurat_object)
  none <- structure(rep(FALSE, n_cells), names = colnames(seurat_object),
                       groups_tested = character(0))

  if (!assay %in% names(seurat_object@assays)) return(none)

  counts <- tryCatch(
    SeuratObject::LayerData(seurat_object, assay = assay, layer = "counts"),
    error = function(e) tryCatch(
      SeuratObject::GetAssayData(seurat_object, assay = assay, slot = "counts"),
      error = function(e2) NULL))
  if (is.null(counts) || nrow(counts) == 0) return(none)

  features <- rownames(counts)
  groups <- .adt_lineage_groups()

  pos_list <- list()
  for (g in names(groups)) {
    feats <- .adt_match_features(features, groups[[g]])
    if (length(feats) == 0) next
    sub <- tryCatch(as.matrix(counts[feats, , drop = FALSE]),
                    error = function(e) NULL)
    if (is.null(sub)) next
    lg <- log1p(sub)

    grp_pos <- rep(FALSE, n_cells)
    for (i in seq_len(nrow(lg))) {
      v <- lg[i, ]
      th <- suppressWarnings(stats::quantile(v, probs = positive_quantile,
                                            na.rm = TRUE, names = FALSE))
      if (!is.finite(th) || th <= 0) next
      grp_pos <- grp_pos | (v > th)
    }
    pos_list[[g]] <- grp_pos
  }

  if (length(pos_list) < min_groups) {
    if (verbose)
      cat(sprintf(
        " ADT cross-lineage check skipped (only %d testable lineage group(s) in panel).\n",
        length(pos_list)))
    return(none)
  }

  pos_mat <- do.call(cbind, pos_list)
  n_pos <- rowSums(pos_mat)
  flagged <- n_pos >= min_groups

  if (verbose)
    cat(sprintf(
      " ADT cross-lineage check: groups tested = %s | flagged %d of %d cells (%.2f%%)\n",
      paste(names(pos_list), collapse = ", "),
      sum(flagged), n_cells, 100 * sum(flagged) / max(1L, n_cells)))

  structure(as.logical(flagged), names = colnames(seurat_object),
            groups_tested = names(pos_list))
}

.detect_doublets <- function(seurat_object,
                             sample_column = NULL,
                             assay = "RNA",
                             adt_assay = "ADT",
                             use_adt_signal = TRUE,
                             adt_positive_quantile = 0.90,
                             min_cells_per_sample = 50L,
                             random_seed = 42L,
                             verbose = TRUE) {

  n_cells <- ncol(seurat_object)
  cells <- colnames(seurat_object)

  score_v <- rep(NA_real_, n_cells); names(score_v) <- cells
  class_v <- rep(NA_character_, n_cells); names(class_v) <- cells

  have_scdbl <- requireNamespace("scDblFinder", quietly = TRUE) &&
                requireNamespace("SingleCellExperiment", quietly = TRUE)

  if (!have_scdbl && verbose)
    cat(" scDblFinder / SingleCellExperiment not installed - RNA doublet detection skipped.\n")

  meta <- seurat_object@meta.data
  if (is.null(sample_column) || !nzchar(sample_column) ||
      !sample_column %in% colnames(meta)) {
    fallback <- if ("orig.ident" %in% colnames(meta)) "orig.ident" else NULL
    if (!is.null(sample_column) && nzchar(sample_column) && verbose)
      cat(sprintf(" Warning: doublet sample column '%s' not found.\n", sample_column))
    sample_column <- fallback
  }

  if (is.null(sample_column)) {
    if (verbose)
      cat(paste0(" Warning: no sample/capture column available - running detection on the\n",
                 " whole object. Calls will be less reliable than per-capture detection.\n"))
    groups <- setNames(rep("all", n_cells), cells)
  } else {
    groups <- setNames(as.character(meta[[sample_column]]), cells)
    groups[is.na(groups)] <- "unassigned"
    if (verbose)
      cat(sprintf(" Per-capture detection using column '%s' (%d capture(s)).\n",
                  sample_column, length(unique(groups))))
  }

  counts_all <- tryCatch(
    SeuratObject::LayerData(seurat_object, assay = assay, layer = "counts"),
    error = function(e) tryCatch(
      SeuratObject::GetAssayData(seurat_object, assay = assay, slot = "counts"),
      error = function(e2) NULL))
  if (is.null(counts_all)) {
    if (verbose) cat(sprintf(" No counts layer on assay '%s' - RNA detection skipped.\n", assay))
    have_scdbl <- FALSE
  }

  per_sample <- list()
  if (have_scdbl) {
    for (g in sort(unique(groups))) {
      idx <- which(groups == g)
      if (length(idx) < min_cells_per_sample) {
        if (verbose)
          cat(sprintf(" %s: %d cells - below min_cells_per_sample (%d), skipped.\n",
                      g, length(idx), min_cells_per_sample))
        per_sample[[g]] <- list(n = length(idx), n_doublet = NA_integer_, skipped = TRUE)
        next
      }
      res <- tryCatch({

        mat <- counts_all[, idx, drop = FALSE]
        if (inherits(mat, "IterableMatrix"))
          mat <- methods::as(mat, "dgCMatrix")
        if (!inherits(mat, "CsparseMatrix"))
          mat <- methods::as(methods::as(mat, "Matrix"), "CsparseMatrix")
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
        set.seed(random_seed)
        sce <- suppressWarnings(scDblFinder::scDblFinder(sce, verbose = FALSE))
        out <- list(
          class = as.character(sce$scDblFinder.class),
          score = as.numeric(sce$scDblFinder.score))
        rm(mat, sce); gc()
        out
      }, error = function(e) {
        if (verbose)
          cat(sprintf(" %s: scDblFinder failed (%s) - capture left unclassified.\n",
                      g, conditionMessage(e)))
        NULL
      })

      if (is.null(res)) {
        per_sample[[g]] <- list(n = length(idx), n_doublet = NA_integer_, skipped = TRUE)
        next
      }
      class_v[idx] <- res$class
      score_v[idx] <- res$score
      nd <- sum(res$class == "doublet", na.rm = TRUE)
      per_sample[[g]] <- list(n = length(idx), n_doublet = nd, skipped = FALSE)
      if (verbose)
        cat(sprintf(" %s: %d cells | %d doublet(s) (%.2f%%)\n",
                    g, length(idx), nd, 100 * nd / length(idx)))
    }
  }

  adt_flag <- rep(FALSE, n_cells); names(adt_flag) <- cells
  if (isTRUE(use_adt_signal) && adt_assay %in% names(seurat_object@assays)) {
    adt_flag <- .adt_crosslineage_doublets(
      seurat_object, assay = adt_assay,
      positive_quantile = adt_positive_quantile, verbose = verbose)
  }

  rna_doublet <- !is.na(class_v) & class_v == "doublet"

  consensus <- rna_doublet | adt_flag
  method <- ifelse(rna_doublet & adt_flag, "RNA+ADT",
               ifelse(rna_doublet, "RNA",
               ifelse(adt_flag, "ADT", NA_character_)))

  seurat_object$doublet_score <- as.numeric(score_v)
  seurat_object$doublet_class_rna <- as.character(class_v)
  seurat_object$adt_crosslineage <- as.logical(adt_flag)
  seurat_object$doublet_class <- ifelse(consensus, "doublet", "singlet")
  seurat_object$doublet_method <- as.character(method)

  summary_df <- data.frame(
    sample = names(per_sample),
    n_cells = vapply(per_sample, function(x) as.integer(x$n), integer(1)),
    n_doublet = vapply(per_sample, function(x) as.integer(x$n_doublet), integer(1)),
    skipped = vapply(per_sample, function(x) isTRUE(x$skipped), logical(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  if (nrow(summary_df) > 0)
    summary_df$pct_doublet <- round(100 * summary_df$n_doublet / summary_df$n_cells, 3)

  seurat_object@misc$doublet_qc <- list(
    sample_column = sample_column,
    per_sample = summary_df,
    n_rna_doublet = sum(rna_doublet),
    n_adt_flagged = sum(adt_flag),
    n_total_doublet = sum(consensus),
    n_cells = n_cells,
    adt_groups_tested = attr(adt_flag, "groups_tested") %||% character(0),
    scdblfinder_available = have_scdbl)

  if (verbose)
    cat(sprintf(
      paste0(" Totals: RNA doublets %d | ADT cross-lineage %d | ",
             "consensus %d of %d cells (%.2f%%)\n"),
      sum(rna_doublet), sum(adt_flag), sum(consensus), n_cells,
      100 * sum(consensus) / max(1L, n_cells)))

  seurat_object
}

.remove_doublets <- function(seurat_object, verbose = TRUE) {
  if (!"doublet_class" %in% colnames(seurat_object@meta.data)) return(seurat_object)
  keep <- colnames(seurat_object)[
    !is.na(seurat_object$doublet_class) & seurat_object$doublet_class == "singlet"]
  n_before <- ncol(seurat_object)
  if (length(keep) < 2) {
    if (verbose)
      cat(sprintf(
        " Doublet removal skipped: would retain only %d cell(s) of %d - check thresholds.\n",
        length(keep), n_before))
    return(seurat_object)
  }
  if (length(keep) == n_before) {
    if (verbose) cat(" No doublets to remove.\n")
    return(seurat_object)
  }
  summ <- seurat_object@misc$doublet_qc
  seurat_object <- subset(seurat_object, cells = keep)

  if (!is.null(summ)) seurat_object@misc$doublet_qc <- summ
  if (verbose)
    cat(sprintf(" Removed %d doublet(s). Cells: %d -> %d (%.1f%% retained)\n",
                n_before - length(keep), n_before, length(keep),
                100 * length(keep) / n_before))
  seurat_object
}

.adt_degenerate_markers <- function(mat, samples = NULL, var_eps = 1e-8,
                                    min_distinct = 10L) {
  if (ncol(mat) == 0) return(character(0))
  use_ms <- requireNamespace("matrixStats", quietly = TRUE)
  bad <- stats::setNames(rep(FALSE, ncol(mat)), colnames(mat))

  gvar <- if (use_ms) matrixStats::colVars(as.matrix(mat)) else apply(mat, 2, stats::var)
  bad <- bad | (colSums(abs(mat), na.rm = TRUE) == 0) | !is.finite(gvar) | (gvar <= var_eps)

  if (!is.null(samples)) {
    for (g in unique(samples)) {
      rows <- which(samples == g)

      if (length(rows) < max(2L, min_distinct)) next
      sub <- mat[rows, , drop = FALSE]

      mads <- if (use_ms) matrixStats::colMads(as.matrix(sub)) else apply(sub, 2, stats::mad)
      mad0 <- !is.finite(mads) | (mads <= .Machine$double.eps)

      ndist <- apply(sub, 2, function(z) length(unique(z[is.finite(z)])))
      lowc <- ndist < min_distinct

      bad <- bad | mad0 | lowc
    }
  }
  names(bad)[bad]
}

.NON_MARKER_TAG_PATTERN <- paste0(
  "^HTO[-_. ]?[0-9]*$",
  "|^HTO[-_. ]",
  "|hash[-_. ]?tag",
  "|hashing",
  "|^HT[0-9]+$",
  "|^BC[0-9]+$",
  "|totalseq.*hashtag",
  "|^Set[0-9]+[-_. ]Sample[0-9]+$",
  "|^Sample[-_. ]?[0-9]+$",
  "|^CMO[-_. ]?[0-9]*$",
  "|^LMO[-_. ]?[0-9]*$",
  "|^MULTI[-_. ]"
)

.detect_nonmarker_tags <- function(features) {
  if (length(features) == 0) return(character(0))
  features[grepl(.NON_MARKER_TAG_PATTERN, features,
                 ignore.case = TRUE, perl = TRUE)]
}

.clean_adt_features <- function(seurat_object, keep_pattern = NULL,
                                drop_pattern = "isotype", batch_column = NULL,
                                verbose = TRUE) {
  if (!"ADT" %in% names(seurat_object@assays)) return(seurat_object)
  adt_feats <- rownames(seurat_object[["ADT"]])
  keep_adt <- adt_feats

  dropped_tags <- .detect_nonmarker_tags(keep_adt)
  if (length(dropped_tags) > 0) {
    if (length(setdiff(keep_adt, dropped_tags)) >= 2) {
      keep_adt <- setdiff(keep_adt, dropped_tags)
      if (verbose)
        cat(sprintf(
          "\n--- Removing %d sample-multiplexing tag(s) from the ADT panel ---\n %s\n",
          length(dropped_tags),
          paste(utils::head(dropped_tags, 20), collapse = ", ")))
    } else {
      if (verbose)
        cat(sprintf(
          "\nMultiplexing-tag removal skipped: would leave only %d ADT feature(s).\n",
          length(setdiff(keep_adt, dropped_tags))))
      dropped_tags <- character(0)
    }
  }

  if (!is.null(keep_pattern) && nzchar(keep_pattern))
    keep_adt <- keep_adt[grepl(keep_pattern, keep_adt, ignore.case = TRUE, perl = TRUE)]
  if (!is.null(drop_pattern) && nzchar(drop_pattern))
    keep_adt <- keep_adt[!grepl(drop_pattern, keep_adt, ignore.case = TRUE, perl = TRUE)]

  dropped_pattern <- setdiff(adt_feats, keep_adt)

  dropped_var <- character(0)
  if (length(keep_adt) >= 2) {
    cnt <- t(as.matrix(
      GetAssayData(seurat_object, assay = "ADT", layer = "counts")[keep_adt, , drop = FALSE]))
    samples <- if (!is.null(batch_column) && batch_column %in% colnames(seurat_object@meta.data))
                 as.character(seurat_object@meta.data[[batch_column]]) else NULL
    dropped_var <- .adt_degenerate_markers(cnt, samples = samples)
    keep_adt <- setdiff(keep_adt, dropped_var)
  }

  dropped_adt <- c(dropped_pattern, dropped_var)
  if (length(dropped_adt) > 0 && length(keep_adt) >= 2) {
    suppressWarnings(seurat_object[["ADT"]] <- subset(seurat_object[["ADT"]], features = keep_adt))
    if (verbose)
      cat(sprintf(
        "\n--- Cleaning ADT panel ---\nADT features: %d -> %d (dropped %d: %d by pattern, %d zero/near-zero-variance)%s\n",
        length(adt_feats), length(keep_adt), length(dropped_adt),
        length(dropped_pattern), length(dropped_var),
        if (length(dropped_adt)) paste0("\n ", paste(utils::head(dropped_adt, 20), collapse = ", "),
                                        if (length(dropped_adt) > 20) ", ..." else "") else ""))
  } else if (length(keep_adt) < 2 && verbose) {
    cat(sprintf(
      "\nADT cleaning skipped: filters would leave only %d feature(s).\n",
      length(keep_adt)))
  }
  seurat_object
}

.normalize_adt <- function(seurat_object, verbose = TRUE) {
  if (!"ADT" %in% names(seurat_object@assays)) return(seurat_object)
  cnt <- GetAssayData(seurat_object, assay = "ADT", layer = "counts")

  v <- tryCatch(utils::head(cnt@x, 10000L), error = function(e) NULL)
  if (is.null(v) || !length(v))
    v <- as.numeric(as.matrix(
      cnt[seq_len(min(nrow(cnt), 32L)), seq_len(min(ncol(cnt), 1000L)), drop = FALSE]))
  v <- v[is.finite(v) & v != 0]
  is_raw <- !length(v) || all(abs(v - round(v)) < 1e-8)

  if (is_raw) {
    if (verbose) cat("Normalizing ADT with CLR...\n")
    seurat_object <- NormalizeData(
      seurat_object, assay = "ADT",
      normalization.method = "CLR", margin = 2, verbose = FALSE)
  } else {
    if (verbose)
      cat(paste0("ADT counts are non-integer (already normalized) - skipping CLR; ",
                 "using the provided values as the ADT 'data' layer.\n"))
    seurat_object <- SetAssayData(
      seurat_object, assay = "ADT", layer = "data", new.data = cnt)
  }
  seurat_object
}

.to_disk_backed <- function(seurat_object, cache_dir = NULL,
                            assays = "RNA", verbose = TRUE) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    if (verbose) cat(" BPCells not installed - keeping counts in memory.\n")
    return(seurat_object)
  }
  if (is.null(cache_dir) || !nzchar(cache_dir))
    cache_dir <- file.path(tempdir(), "bpcells_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  for (assay_name in intersect(assays, names(seurat_object@assays))) {
    seurat_object <- tryCatch({

      counts <- tryCatch(
        SeuratObject::LayerData(seurat_object, assay = assay_name, layer = "counts"),
        error = function(e) tryCatch(
          SeuratObject::GetAssayData(seurat_object, assay = assay_name, slot = "counts"),
          error = function(e2) NULL))
      if (is.null(counts)) {
        if (verbose) cat(sprintf(" %s: no counts layer found - skipping.\n", assay_name))
        seurat_object
      } else if (inherits(counts, "IterableMatrix")) {
        if (verbose) cat(sprintf(" %s counts already on-disk - skipping.\n", assay_name))
        seurat_object
      } else {
        mat_dir <- file.path(cache_dir, sprintf("%s_counts", assay_name))
        if (dir.exists(mat_dir)) unlink(mat_dir, recursive = TRUE)

        BPCells::write_matrix_dir(mat = counts, dir = mat_dir, overwrite = TRUE)
        disk <- BPCells::open_matrix_dir(dir = mat_dir)
        suppressWarnings(seurat_object[[assay_name]] <- SeuratObject::CreateAssay5Object(counts = disk))
        if (verbose) cat(sprintf(" %s counts moved on-disk (BPCells): %s\n",
                                 assay_name, mat_dir))
        gc()
        seurat_object
      }
    }, error = function(e) {
      if (verbose) cat(sprintf(" Disk-backing %s failed (%s) - kept in memory.\n",
                               assay_name, conditionMessage(e)))
      seurat_object
    })
  }
  gc()
  return(seurat_object)
}

.compute_rna_pca_full <- function(seurat_object, assay_name,
                                   n_variable_features, verbose) {
  DefaultAssay(seurat_object) <- assay_name

  if (assay_name != "SCT") {
    if (verbose) cat("Finding variable features (RNA)...\n")
    seurat_object <- FindVariableFeatures(
      seurat_object,
      selection.method = "vst",
      nfeatures = n_variable_features,
      verbose = FALSE)
  }

  if (verbose) cat("Scaling RNA data...\n")
  seurat_object <- ScaleData(
    seurat_object,
    features = VariableFeatures(seurat_object),
    verbose = FALSE)
  gc()

  if (verbose) cat("Running PCA on RNA...\n")
  seurat_object <- RunPCA(
    seurat_object,
    features = VariableFeatures(seurat_object),
    npcs = 50,
    verbose = FALSE)

  return(seurat_object)
}

.compute_rna_pca_sketch <- function(seurat_object, assay_name, n_variable_features,
                                     sketch_ncells, sketch_method, verbose) {
  DefaultAssay(seurat_object) <- assay_name

  if (assay_name != "SCT") {
    if (verbose) cat("Finding variable features (RNA)...\n")
    seurat_object <- FindVariableFeatures(
      seurat_object,
      selection.method = "vst",
      nfeatures = n_variable_features,
      verbose = FALSE)
  }

  if (verbose)
    cat(sprintf("Creating sketch (%s, %d cells)...\n", sketch_method, sketch_ncells))
  seurat_object <- SketchData(
    object = seurat_object,
    ncells = sketch_ncells,
    method = sketch_method,
    sketched.assay = "sketch",
    assay = assay_name)

  DefaultAssay(seurat_object) <- "sketch"

  if (verbose) cat("Scaling sketch data...\n")
  seurat_object <- ScaleData(
    seurat_object,
    features = VariableFeatures(seurat_object),
    verbose = FALSE)

  if (verbose) cat("Running PCA on sketch...\n")
  seurat_object <- RunPCA(
    seurat_object,
    features = VariableFeatures(seurat_object),
    npcs = 50,
    verbose = FALSE)
  gc()

  if (verbose) cat("Projecting PCA to full dataset...\n")
  seurat_object <- ProjectData(
    object = seurat_object,
    assay = assay_name,
    full.reduction = "full.pca",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    dims = 1:50,
    verbose = verbose)

  DefaultAssay(seurat_object) <- assay_name

  return(seurat_object)
}

.compute_embeddings <- function(seurat_object, has_adt, rna_reduction, assay_name,
                                 pca_dims_rna, pca_dims_adt, clustering_resolution,
                                 wnn_resolution, k_param, clustering_algorithm,
                                 verbose,

                                 adt_harmony = FALSE,
                                 batch_column = NULL,
                                 harmony_theta = 2,
                                 harmony_niter = 10,
                                 adt_guardrail_min_improvement = 0.05) {

  DefaultAssay(seurat_object) <- assay_name

  cluster_method <- "matrix"
  if (identical(as.integer(clustering_algorithm), 4L)) {
    leiden_ok <- requireNamespace("reticulate", quietly = TRUE) &&
                 isTRUE(tryCatch(reticulate::py_module_available("leidenalg"),
                                 error = function(e) FALSE))
    if (leiden_ok) {
      cluster_method <- "igraph"
      if (verbose) cat(" Leiden clustering via igraph (sparse graph).\n")
    } else {
      clustering_algorithm <- 1L
      if (verbose)
        cat(" 'leidenalg' not available - falling back to Louvain (algorithm 1).\n")
    }
  }

  if (verbose) cat("Running UMAP on RNA...\n")
  seurat_object <- RunUMAP(
    seurat_object,
    reduction = rna_reduction,
    dims = 1:pca_dims_rna,
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_",
    verbose = FALSE)

  if (verbose) cat("Clustering based on RNA...\n")
  seurat_object <- FindNeighbors(
    seurat_object,
    reduction = rna_reduction,
    dims = 1:pca_dims_rna,
    k.param = k_param,
    verbose = FALSE)

  seurat_object <- FindClusters(
    seurat_object,
    resolution = clustering_resolution,
    cluster.name = "RNA_clusters",
    algorithm = clustering_algorithm,
    method = cluster_method,
    verbose = FALSE)

  if (verbose)
    cat(sprintf("RNA clusters identified: %d\n",
                length(unique(seurat_object$RNA_clusters))))
  gc()

  if (has_adt) {
    if (verbose) cat("\n--- ADT Analysis ---\n")
    DefaultAssay(seurat_object) <- "ADT"

    if (verbose) cat("Scaling ADT data...\n")
    seurat_object <- ScaleData(
      seurat_object,
      features = rownames(seurat_object[["ADT"]]),
      verbose = FALSE)
    gc()

    if (verbose) cat("Running PCA on ADT...\n")
    seurat_object <- RunPCA(
      seurat_object,
      assay = "ADT",
      features = rownames(seurat_object[["ADT"]]),
      reduction.name = "apca",
      reduction.key = "apca_",
      npcs = 30,
      verbose = FALSE)

    adt_reduction <- "apca"
    if (isTRUE(adt_harmony) &&
        (is.null(batch_column) ||
         !batch_column %in% colnames(seurat_object@meta.data))) {

      if (verbose) cat(sprintf(
        paste0(" ADT Harmony correction SKIPPED: batch column '%s' not found ",
               "in metadata - using uncorrected apca for WNN/UMAP.\n"),
        batch_column %||% "NULL"))
      seurat_object@misc$adt_harmony_qc <- list(
        method = "harmony", failed = TRUE, adopted = FALSE,
        reason = sprintf("batch column '%s' not found in metadata",
                         batch_column %||% "NULL"),
        kbet_before = NA_real_, kbet_after = NA_real_,
        reduction_used = adt_reduction)
    }
    if (isTRUE(adt_harmony) && !is.null(batch_column) &&
        batch_column %in% colnames(seurat_object@meta.data)) {
      if (verbose) cat("\n--- ADT Batch Correction (Harmony on apca) ---\n")
      labels_adt <- seurat_object@meta.data[[batch_column]]
      kbet_before <- .kbet_rejection(Embeddings(seurat_object, "apca"), labels_adt)

      harmony_ok <- tryCatch({
        seurat_object <- .correct_adt_batch_harmony(
          seurat_object = seurat_object,
          batch_column = batch_column,
          theta = harmony_theta,
          niter = harmony_niter,
          verbose = FALSE)
        TRUE
      }, error = function(e) {
        if (verbose) cat(sprintf(
          " Harmony on apca FAILED (%s) - continuing with uncorrected apca.\n",
          conditionMessage(e)))
        FALSE
      })

      if (harmony_ok) {
        kbet_after <- .kbet_rejection(Embeddings(seurat_object, "harmony.adt"),
                                      labels_adt)
        adopt <- .adt_guardrail_adopt(kbet_before, kbet_after,
                                      adt_guardrail_min_improvement)
        if (adopt) {
          adt_reduction <- "harmony.adt"
          if (verbose) cat(sprintf(
            " Guardrail: batch mixing improved (kBET %.3f -> %.3f) - using harmony.adt for WNN/UMAP.\n",
            kbet_before, kbet_after))
        } else if (verbose) {
          cat(sprintf(
            paste0(" GUARDRAIL: Harmony did not improve batch mixing ",
                   "(kBET %s -> %s) - REVERTED to uncorrected apca.\n"),
            if (is.finite(kbet_before)) sprintf("%.3f", kbet_before) else "NA",
            if (is.finite(kbet_after)) sprintf("%.3f", kbet_after) else "NA"))
        }
        seurat_object@misc$adt_harmony_qc <- list(
          method = "harmony",
          kbet_before = kbet_before,
          kbet_after = kbet_after,
          min_improvement = adt_guardrail_min_improvement,
          adopted = adopt,
          reduction_used = adt_reduction)
      } else {
        seurat_object@misc$adt_harmony_qc <- list(
          method = "harmony", failed = TRUE, adopted = FALSE,
          kbet_before = kbet_before, kbet_after = NA_real_,
          reduction_used = adt_reduction)
      }
    }

    if (verbose) cat("Running UMAP on ADT...\n")
    seurat_object <- RunUMAP(
      seurat_object,
      reduction = adt_reduction,
      dims = 1:pca_dims_adt,
      reduction.name = "umap.adt",
      reduction.key = "adtUMAP_",
      verbose = FALSE)

    if (verbose) cat("Clustering based on ADT...\n")
    seurat_object <- FindNeighbors(
      seurat_object,
      reduction = adt_reduction,
      dims = 1:pca_dims_adt,
      graph.name = c("adt_nn", "adt_snn"),
      k.param = k_param,
      verbose = FALSE)

    seurat_object <- FindClusters(
      seurat_object,
      graph.name = "adt_snn",
      resolution = clustering_resolution,
      cluster.name = "ADT_clusters",
      algorithm = clustering_algorithm,
      method = cluster_method,
      verbose = FALSE)

    if (verbose)
      cat(sprintf("ADT clusters identified: %d\n",
                  length(unique(seurat_object$ADT_clusters))))

    DefaultAssay(seurat_object) <- assay_name
    gc()

    if (verbose) cat("\n--- WNN Multimodal Integration ---\n")
    seurat_object <- FindMultiModalNeighbors(
      seurat_object,
      reduction.list = list(rna_reduction, adt_reduction),
      dims.list = list(1:pca_dims_rna, 1:pca_dims_adt),
      modality.weight.name = c("RNA.weight", "ADT.weight"),
      k.nn = k_param,
      verbose = verbose)
    gc()

    if (verbose) cat("Clustering based on WNN...\n")
    seurat_object <- FindClusters(
      seurat_object,
      graph.name = "wsnn",
      resolution = wnn_resolution,
      cluster.name = "wsnn_clusters",
      algorithm = clustering_algorithm,
      method = cluster_method,
      verbose = FALSE)

    if (verbose)
      cat(sprintf("WNN clusters identified: %d\n",
                  length(unique(seurat_object$wsnn_clusters))))
    gc()

    if (verbose) cat("Running UMAP on WNN graph...\n")
    seurat_object <- RunUMAP(
      seurat_object,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_",
      verbose = FALSE)
  }

  return(seurat_object)
}

.detect_rna_batch <- function(seurat_object, rna_reduction, batch_column,
                               n_batches, lisi_threshold, silhouette_threshold,
                               silhouette_n_cells = 5000L, verbose) {
  result <- list(
    median_lisi = NA_real_,
    mean_silhouette = NA_real_,
    lisi_flag = FALSE,
    silhouette_flag = FALSE,
    correction_recommended = FALSE
  )

  emb <- Embeddings(seurat_object, reduction = rna_reduction)
  meta <- seurat_object@meta.data

  tryCatch({
    max_cells <- min(nrow(emb), silhouette_n_cells)
    set.seed(42L)
    idx <- sample(nrow(emb), max_cells)
    emb_sub <- emb[idx, , drop = FALSE]
    batch_int <- as.integer(factor(meta[[batch_column]][idx]))
    d_sub <- dist(emb_sub)
    sil <- cluster::silhouette(batch_int, d_sub)
    result$mean_silhouette <- mean(sil[, "sil_width"], na.rm = TRUE)
    result$silhouette_flag <- result$mean_silhouette > silhouette_threshold
    if (verbose)
      cat(sprintf(
        " Silhouette (batch, n = %d subsample): mean = %.3f (threshold = %.3f) - %s\n",
        max_cells, result$mean_silhouette, silhouette_threshold,
        ifelse(result$silhouette_flag, "FLAGGED", "OK")))
  }, error = function(e) {
    if (verbose)
      cat(sprintf(" Silhouette computation failed: %s\n", conditionMessage(e)))
  })

  result$correction_recommended <-
    isTRUE(result$lisi_flag) || isTRUE(result$silhouette_flag)
  return(result)
}

.detect_adt_batch <- function(seurat_object, batch_column, n_batches,
                               wasserstein_threshold, verbose,
                               save_dir = NULL, save_prefix = "adt",
                               kbet_threshold = 0.6) {
  result <- list(
    flagged_proteins = character(0),
    pct_proteins_flagged = 0,
    mean_wasserstein = NA_real_,
    kbet_rejection = NA_real_,
    kbet_threshold = kbet_threshold,
    kbet_flagged = FALSE,
    correction_recommended = FALSE
  )

  adt_sparse <- GetAssayData(seurat_object, assay = "ADT", layer = "data")
  batch_labels <- seurat_object@meta.data[[batch_column]]
  batches <- sort(unique(na.omit(batch_labels)))
  markers <- rownames(adt_sparse)
  batch_pairs <- combn(batches, 2, simplify = FALSE)

  pair_names <- vapply(batch_pairs, function(p) paste(p, collapse = "_vs_"), character(1))

  if (verbose)
    cat(sprintf(
      " Computing per-protein distances: %d markers x %d batch pairs\n",
      length(markers), length(batch_pairs)))

  use_transport <- requireNamespace("transport", quietly = TRUE)
  if (!use_transport && verbose)
    cat(" 'transport' not installed - using KS statistic as distance proxy.\n")

  wass_mat <- matrix(
    NA_real_,
    nrow = length(markers),
    ncol = length(batch_pairs),
    dimnames = list(markers, pair_names))

  for (m in markers) {
    expr_m <- as.numeric(adt_sparse[m, ])
    for (j in seq_along(batch_pairs)) {
      b1 <- batch_pairs[[j]][1]
      b2 <- batch_pairs[[j]][2]
      x1 <- expr_m[!is.na(batch_labels) & batch_labels == b1]
      x2 <- expr_m[!is.na(batch_labels) & batch_labels == b2]
      x1 <- x1[!is.na(x1)]
      x2 <- x2[!is.na(x2)]
      if (length(x1) < 5L || length(x2) < 5L) next
      wass_mat[m, j] <- if (use_transport) {
        transport::wasserstein1d(x1, x2)
      } else {

        suppressWarnings(ks.test(x1, x2)$statistic)
      }
    }
  }

  per_marker_mean <- rowMeans(wass_mat, na.rm = TRUE)
  flagged <- names(per_marker_mean)[
    !is.na(per_marker_mean) & per_marker_mean > wasserstein_threshold]

  result$flagged_proteins <- flagged
  result$pct_proteins_flagged <- length(flagged) / length(markers)
  result$mean_wasserstein <- mean(per_marker_mean, na.rm = TRUE)

  qc <- .adt_qc_pca_scores(seurat_object)
  if (!is.null(qc)) {
    result$kbet_rejection <- .kbet_rejection(qc$emb, batch_labels[qc$idx])
    result$kbet_flagged <- is.finite(result$kbet_rejection) &&
                             result$kbet_rejection > kbet_threshold
  }

  result$correction_recommended <-
    (length(flagged) > 0L && result$pct_proteins_flagged > 0.10) ||
    isTRUE(result$kbet_flagged)

  wass_str <- if (!is.na(result$mean_wasserstein))
    sprintf("%.3f", result$mean_wasserstein) else "NA"
  kbet_str <- if (is.finite(result$kbet_rejection))
    sprintf("%.3f (threshold %.2f)", result$kbet_rejection, kbet_threshold) else "NA"

  if (verbose)
    cat(sprintf(
      " Mean distance: %s | Flagged: %d / %d proteins (%.0f%%) | kBET rejection: %s - %s\n",
      wass_str,
      length(flagged), length(markers),
      100 * result$pct_proteins_flagged,
      kbet_str,
      ifelse(result$correction_recommended, "FLAGGED", "OK")))

  if (verbose && length(flagged)) {
    ord <- order(per_marker_mean[flagged], decreasing = TRUE)
    cat(sprintf(
      " Flagged proteins (mean distance): %s\n",
      paste(sprintf("%s (%.3f)", flagged[ord], per_marker_mean[flagged][ord]),
            collapse = ", ")))
  }

  result$per_protein_distance <- sort(per_marker_mean, decreasing = TRUE)

  if (!is.null(save_dir) && nzchar(save_dir)) {
    ok <- dir.exists(save_dir) ||
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    if (ok) {
      summ <- data.frame(
        protein = markers,
        mean_distance = as.numeric(per_marker_mean),
        n_pairs = rowSums(!is.na(wass_mat)),
        flagged = markers %in% flagged,
        stringsAsFactors = FALSE)
      summ <- summ[order(summ$mean_distance, decreasing = TRUE), ]
      f_summary <- file.path(save_dir, paste0(save_prefix, "_per_protein_distance.csv"))
      f_pairwise <- file.path(save_dir, paste0(save_prefix, "_pairwise_distance.csv"))
      tryCatch({
        utils::write.csv(summ, f_summary, row.names = FALSE)
        utils::write.csv(
          data.frame(protein = markers, wass_mat, check.names = FALSE),
          f_pairwise, row.names = FALSE)
        if (verbose)
          cat(sprintf(" Per-protein distances written to: %s\n", f_summary))
      }, error = function(e) if (verbose)
        cat(sprintf(" Could not write per-protein distance CSV: %s\n", e$message)))
    }
  }

  return(result)
}

.kbet_rejection <- function(emb, batch_labels, k = 30L, n_sample = 3000L,
                            seed = 42L, alpha = 0.05) {
  ok <- !is.na(batch_labels)
  emb <- emb[ok, , drop = FALSE]
  batch_labels <- as.character(batch_labels[ok])
  n <- nrow(emb)
  batches <- unique(batch_labels)
  if (length(batches) < 2L || n < (k + 1L) * 2L) return(NA_real_)

  set.seed(seed)
  idx <- if (n > n_sample) sort(sample.int(n, n_sample)) else seq_len(n)
  emb_s <- emb[idx, , drop = FALSE]
  batch_s <- batch_labels[idx]
  ns <- length(idx)
  k <- min(k, ns - 1L)

  glob <- table(factor(batch_s, levels = batches)) / ns

  d <- as.matrix(stats::dist(emb_s))
  diag(d) <- Inf

  rej <- 0L; tested <- 0L
  expected <- as.numeric(glob) * k
  for (i in seq_len(ns)) {
    nn <- order(d[i, ])[seq_len(k)]
    obs <- table(factor(batch_s[nn], levels = batches))

    stat <- sum((as.numeric(obs) - expected)^2 / pmax(expected, 1e-9))
    p <- stats::pchisq(stat, df = length(batches) - 1L, lower.tail = FALSE)
    tested <- tested + 1L
    if (p < alpha) rej <- rej + 1L
  }
  if (tested == 0L) return(NA_real_)
  rej / tested
}

.adt_qc_pca_scores <- function(seurat_object, dims = 20L, n_sample = 3000L,
                               seed = 42L) {
  adt <- tryCatch(GetAssayData(seurat_object, assay = "ADT", layer = "data"),
                  error = function(e) NULL)
  if (is.null(adt) || nrow(adt) < 3L) return(NULL)

  n <- ncol(adt)
  set.seed(seed)
  idx <- if (n > n_sample) sort(sample.int(n, n_sample)) else seq_len(n)

  m <- as.matrix(adt[, idx, drop = FALSE])
  sds <- apply(m, 1, stats::sd)
  m <- m[is.finite(sds) & sds > 0, , drop = FALSE]
  if (nrow(m) < 3L) return(NULL)

  dims <- min(as.integer(dims), nrow(m) - 1L, length(idx) - 1L)
  pc <- tryCatch(stats::prcomp(t(m), center = TRUE, scale. = TRUE, rank. = dims),
                 error = function(e) NULL)
  if (is.null(pc)) return(NULL)
  list(emb = pc$x, idx = idx)
}

.adt_guardrail_adopt <- function(kbet_before, kbet_after, min_improvement = 0.05) {
  if (!is.finite(kbet_before) || !is.finite(kbet_after)) return(FALSE)

  (kbet_before - kbet_after) >= (min_improvement - 1e-9)
}

.correct_adt_batch_harmony <- function(seurat_object, batch_column,
                                       theta, niter, verbose) {
  if (!requireNamespace("harmony", quietly = TRUE))
    stop(
      "Package 'harmony' is required for ADT batch correction (method 'harmony').\n",
      "Install with: install.packages('harmony')")

  rh <- tryCatch(utils::getS3method("RunHarmony", "Seurat"),
                    error = function(e) harmony::RunHarmony)
  fmls <- names(formals(rh))
  red_name <- if ("reduction.use" %in% fmls) "reduction.use" else "reduction"
  iter_name <- if (utils::packageVersion("harmony") >= "1.0.0") "max_iter" else "max.iter.harmony"

  rh_args <- list(
    object = seurat_object,
    group.by.vars = batch_column,
    reduction.save = "harmony.adt",
    theta = theta,
    verbose = verbose)
  rh_args[[red_name]] <- "apca"
  rh_args[[iter_name]] <- niter
  if ("project.dim" %in% fmls) rh_args$project.dim <- FALSE

  seurat_object <- do.call(harmony::RunHarmony, rh_args)
  return(seurat_object)
}

.correct_rna_batch_sketch <- function(seurat_object, assay_name, batch_column,
                                      theta, niter, pca_dims, verbose) {
  if (!requireNamespace("harmony", quietly = TRUE))
    stop("Package 'harmony' is required for sketch batch correction.")
  if (!"sketch" %in% names(seurat_object@assays) ||
      !"pca" %in% names(seurat_object@reductions))
    stop("Sketch assay or sketch 'pca' reduction not found.")

  DefaultAssay(seurat_object) <- "sketch"
  rh <- tryCatch(utils::getS3method("RunHarmony", "Seurat"),
                    error = function(e) harmony::RunHarmony)
  fmls <- names(formals(rh))
  red_name <- if ("reduction.use" %in% fmls) "reduction.use" else "reduction"
  iter_name <- if (utils::packageVersion("harmony") >= "1.0.0") "max_iter" else "max.iter.harmony"
  rh_args <- list(
    object = seurat_object,
    group.by.vars = batch_column,
    reduction.save = "harmony",
    theta = theta,
    verbose = verbose)
  rh_args[[red_name]] <- "pca"
  rh_args[[iter_name]] <- niter

  seurat_object <- do.call(harmony::RunHarmony, rh_args)

  if (verbose) cat(" Projecting Harmony-corrected sketch embedding to all cells...\n")
  seurat_object <- ProjectData(
    object = seurat_object,
    assay = assay_name,
    sketched.assay = "sketch",
    sketched.reduction = "harmony",
    full.reduction = "harmony.full",
    dims = 1:pca_dims,
    verbose = verbose)

  DefaultAssay(seurat_object) <- assay_name
  return(seurat_object)
}

.correct_rna_batch <- function(seurat_object, rna_reduction, batch_column,
                                theta, niter, verbose) {
  if (!requireNamespace("harmony", quietly = TRUE))
    stop(
      "Package 'harmony' is required for RNA batch correction.\n",
      "Install with: install.packages('harmony')")

  rh <- tryCatch(utils::getS3method("RunHarmony", "Seurat"),
                    error = function(e) harmony::RunHarmony)
  fmls <- names(formals(rh))
  red_name <- if ("reduction.use" %in% fmls) "reduction.use" else "reduction"
  iter_name <- if (utils::packageVersion("harmony") >= "1.0.0") "max_iter" else "max.iter.harmony"

  rh_args <- list(
    object = seurat_object,
    group.by.vars = batch_column,
    reduction.save = "harmony",
    theta = theta,
    verbose = verbose)
  rh_args[[red_name]] <- rna_reduction
  rh_args[[iter_name]] <- niter

  if ("project.dim" %in% fmls) rh_args$project.dim <- FALSE

  seurat_object <- do.call(harmony::RunHarmony, rh_args)

  return(seurat_object)
}

.correct_adt_batch <- function(seurat_object, batch_column, verbose) {
  if (!requireNamespace("ADTnorm", quietly = TRUE))
    stop(
      "Package 'ADTnorm' is required for ADT batch correction.\n",
      "Install with: BiocManager::install('ADTnorm')")

  raw_counts <- t(as.matrix(
    GetAssayData(seurat_object, assay = "ADT", layer = "counts")))

  orig_markers <- colnames(raw_counts)
  safe_markers <- make.unique(gsub("[/\\:]", "-", orig_markers))
  colnames(raw_counts) <- safe_markers

  cell_meta <- data.frame(
    sample = as.character(seurat_object@meta.data[[batch_column]]),
    row.names = colnames(seurat_object),
    stringsAsFactors = FALSE)

  dropped <- .adt_degenerate_markers(raw_counts, samples = cell_meta$sample)
  good <- setdiff(colnames(raw_counts), dropped)
  good_orig <- setNames(orig_markers, safe_markers)[good]
  dropped_orig <- setNames(orig_markers, safe_markers)[dropped]

  if (length(dropped)) {
    if (verbose)
      cat(sprintf(" Removing %d degenerate ADT marker(s) (constant within a batch / all-zero) from the assay: %s\n",
                  length(dropped_orig), paste(dropped_orig, collapse = ", ")))
    if (length(good) < 2L) {
      if (verbose)
        cat(" Too few non-degenerate ADT markers remain; leaving ADT assay unchanged and keeping CLR values.\n")
      return(seurat_object)
    }

    seurat_object[["ADT"]] <- subset(seurat_object[["ADT"]], features = good_orig)
    raw_counts <- raw_counts[, good, drop = FALSE]
  }

  adt_temp <- tempfile("adtnorm_")
  dir.create(adt_temp, showWarnings = FALSE)
  on.exit(unlink(adt_temp, recursive = TRUE), add = TRUE)

  run_adtnorm <- function(mat) {
    res <- NULL
    invisible(utils::capture.output(suppressMessages(
      res <- ADTnorm::ADTnorm(
        cell_x_adt = mat,
        cell_x_feature = cell_meta,
        save_outpath = adt_temp,
        study_name = "batch_correction",
        marker_to_process = NULL))))
    res
  }

  if (verbose)
    cat(sprintf(" Running ADTnorm: %d cells x %d markers...\n",
                nrow(raw_counts), ncol(raw_counts)))

  failed <- character(0)
  norm_mat <- tryCatch(run_adtnorm(raw_counts), error = function(e) NULL)

  if (is.null(norm_mat)) {
    if (verbose)
      cat(" Whole-panel ADTnorm failed; retrying marker-by-marker to isolate unmodellable markers...\n")
    cols <- list()
    for (mk in colnames(raw_counts)) {
      r <- tryCatch(run_adtnorm(raw_counts[, mk, drop = FALSE]), error = function(e) NULL)
      if (is.null(r)) { failed <- c(failed, mk); next }
      cols[[mk]] <- r[, 1]
    }
    if (!length(cols)) {
      if (verbose)
        cat(" ADTnorm could not model any marker; keeping ADT data uncorrected.\n")
      return(seurat_object)
    }
    norm_mat <- do.call(cbind, cols)
    rownames(norm_mat) <- rownames(raw_counts)
  }

  corrected <- t(norm_mat)

  rownames(corrected) <- setNames(orig_markers, safe_markers)[rownames(corrected)]

  if (length(failed)) {
    failed_orig <- setNames(orig_markers, safe_markers)[failed]
    keep_orig <- intersect(rownames(seurat_object[["ADT"]]), rownames(corrected))
    if (length(keep_orig) < 2L) {
      if (verbose)
        cat(" Too few markers survived ADTnorm; keeping ADT data uncorrected.\n")
      return(seurat_object)
    }
    seurat_object[["ADT"]] <- subset(seurat_object[["ADT"]], features = keep_orig)
    if (verbose)
      cat(sprintf(" Removed %d marker(s) ADTnorm could not model: %s\n",
                  length(failed_orig), paste(failed_orig, collapse = ", ")))
  }

  seurat_object <- SetAssayData(
    seurat_object,
    assay = "ADT",
    layer = "data",
    new.data = corrected[rownames(seurat_object[["ADT"]]), colnames(seurat_object), drop = FALSE])

  if (verbose)
    cat(sprintf(" ADTnorm complete - %d marker(s) corrected, %d removed (%d pre-filtered, %d unmodellable).\n",
                nrow(corrected), length(dropped) + length(failed),
                length(dropped), length(failed)))

  return(seurat_object)
}

process_full <- function(
    seurat_object,
    has_adt,
    normalization_method,
    n_variable_features,
    pca_dims_rna,
    pca_dims_adt,
    clustering_resolution,
    wnn_resolution,
    k_param,
    clustering_algorithm = 1,
    verbose
) {
  assay_name <- ifelse(normalization_method == "SCT", "SCT", "RNA")
  seurat_object <- .compute_rna_pca_full(
    seurat_object, assay_name, n_variable_features, verbose)
  seurat_object <- .compute_embeddings(
    seurat_object, has_adt, "pca", assay_name,
    pca_dims_rna, pca_dims_adt, clustering_resolution,
    wnn_resolution, k_param, clustering_algorithm, verbose)
  return(seurat_object)
}

process_with_sketch <- function(
    seurat_object,
    has_adt,
    normalization_method,
    sketch_ncells,
    sketch_method = "LeverageScore",
    n_variable_features,
    pca_dims_rna,
    pca_dims_adt,
    clustering_resolution,
    wnn_resolution,
    k_param,
    clustering_algorithm = 1,
    verbose
) {
  assay_name <- ifelse(normalization_method == "SCT", "SCT", "RNA")
  seurat_object <- .compute_rna_pca_sketch(
    seurat_object, assay_name, n_variable_features,
    sketch_ncells, sketch_method, verbose)
  seurat_object <- .compute_embeddings(
    seurat_object, has_adt, "full.pca", assay_name,
    pca_dims_rna, pca_dims_adt, clustering_resolution,
    wnn_resolution, k_param, clustering_algorithm, verbose)
  return(seurat_object)
}
