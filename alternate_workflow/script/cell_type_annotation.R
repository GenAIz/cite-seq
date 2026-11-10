library(Seurat)
library(SingleR)
library(celldex)
library(BiocParallel)
library(SingleCellExperiment)
library(reticulate)

.sf_fetch_layer <- function(object, assay, layer) {
  m <- tryCatch(LayerData(object, assay = assay, layer = layer),
                error = function(e) NULL)
  if (is.null(m) || length(dim(m)) != 2 || any(dim(m) == 0)) NULL else m
}

#' Cell Type Annotation for Seurat Objects
#'
#' @description
#' Performs cell type annotation using SingleR, STCAT, or both methods.
#' This function provides flexible cell type assignment with automatic
#' reference dataset handling and optional refinement. Includes automatic
#' variable gene detection and subsetting for memory-efficient STCAT processing.
#'
#' ANNOTATION SEQUENCE (SingleR-gated STCAT): SingleR (main + fine) runs
#' first, and only the cells whose SingleR labels reference T cells are sent
#' to STCAT. Non-T cells receive the fallback label (SingleR_fine by default)
#' in the Prediction column and are never exported to Python.
#'
#' @param seurat_object A Seurat object with normalized data
#' @param method Character string specifying annotation method.
#' Options: "SingleR", "STCAT", "both" (default: "both")
#' @param reference Reference dataset for SingleR annotation.
#' Options: "monaco" (default), "hpca", "blueprint_encode", "dice", "immgen",
#' or a custom SingleCellExperiment reference object
#' @param cluster_column Character string specifying the metadata column
#' containing cluster assignments (default: "RNA_clusters")
#' @param assay Character string specifying the assay to use (default: "RNA")
#' @param singler_label_type Character string or vector specifying label types.
#' Options: "main", "fine", or c("main", "fine") (default: c("main", "fine"))
#' @param singler_workers Integer specifying number of parallel workers for
#' SingleR (default: 4)
#' @param singler_fine_chunk_size Integer: number of cells per chunk for the
#' per-cell fine annotation pass. Memory scales with the chunk, not the
#' dataset (default: 20000)
#' @param stcat_python_script Path to STCAT Python script. If NULL, will look
#' for "run_stcat.py" in working directory (default: NULL)
#' @param stcat_batch_column Metadata column identifying batches. When it has
#' more than one level among the gated T cells, the T-subset embedding
#' passed to STCAT is Harmony-corrected on that column so batch structure
#' does not re-enter STCAT's neighbour graph (default: NULL = no correction)
#' @param prefix Character string to prefix annotation column names.
#' Generally not needed unless running multiple annotation workflows
#' (default: NULL for no prefix)
#' @param fallback_column When using STCAT, which SingleR annotation to use
#' for cells STCAT does not annotate: cells outside the T cell gate, gated
#' cells STCAT does not return, and "None T" predictions (default: "SingleR_fine")
#' @param verbose Logical indicating whether to print progress messages
#' (default: TRUE)
#' @param cleanup_temp Logical indicating whether to clean up temporary files
#' after STCAT processing (default: TRUE)
#'
#' @return A Seurat object with added metadata columns:
#' For SingleR:
#' - SingleR_main: Main cell type labels, cluster-level (if requested)
#' - SingleR_fine: Fine cell type labels, PER CELL (if requested)
#' - SingleR_main_score / SingleR_fine_score: confidence scores
#' - vocab_rescue: audit flag naming cells relabelled by the vocabulary
#' rescue (Platelets/Erythrocytes); NA for all other cells
#'
#' For STCAT:
#' - Prediction: STCAT predicted cell types (fallback label for non-T cells);
#' Monaco fallback names that duplicate STCAT concepts are harmonized to
#' the STCAT vocabulary (see .prediction_harmonization_map)
#' - Uncertainty_score: STCAT uncertainty scores. Exact provenance flag:
#' non-NA = label produced by STCAT, NA = label from SingleR (fallback)
#'
#' Working columns (SingleR_fine_cluster*, SingleR_*_delta, Cluster,
#' vocab_rescue_score, seurat_clusters) are removed before return;
#' see .drop_working_columns.
#'
#' @details
#' This function performs cell type annotation using one or both of the
#' following methods:
#'
#' 1. SingleR: Uses reference datasets to assign cell types at cluster level,
#' then transfers annotations to individual cells with confidence scores.
#' Fine labels are subsequently refined PER CELL in bounded chunks
#' after the cluster-level pass; the cluster-level fine labels are kept in
#' SingleR_fine_cluster columns.
#'
#' 2. STCAT: Uses Python-based STCAT for spatial transcriptomics cell type
#' annotation. Requires proper Python environment setup via reticulate.
#' Automatically handles variable gene selection:
#' - If VariableFeatures exist, uses those genes
#' - If >6000 total genes and no VariableFeatures, computes top 5000
#' - If <=6000 genes, uses all genes
#' This reduces memory usage significantly for large datasets.
#'
#' When method = "both", SingleR is run first, then STCAT refines the T cell
#' predictions. Only cells whose SingleR_main OR SingleR_fine label references
#' T cells (see .is_tcell_label) are exported to STCAT; all other cells get
#' the fallback label (SingleR_fine) in Prediction. STCAT's own "None T"
#' answers on gated cells are also replaced with the fallback.
#'
#' If method = "STCAT" and the SingleR columns STCAT depends on (the fallback
#' column and the SingleR_main/SingleR_fine gate columns) don't exist, SingleR
#' will be automatically run first to create them.
#'
#' @note
#' For STCAT functionality:
#' - Requires a properly configured Python environment with STCAT installed
#' - Requires a Python script (run_stcat.py) to execute STCAT analysis
#' - Creates temporary h5ad files for data exchange with Python
#' - Automatically subsets to variable genes to reduce memory usage
#'
#' @examples
#' \dontrun{
#' # SingleR only with Monaco reference
#' seurat_obj <- cell_type_annotation(
#' seurat_object = seurat_obj,
#' method = "SingleR",
#' reference = "monaco"
#' )
#'
#' # STCAT only (will auto-run SingleR if fallback not found)
#' seurat_obj <- cell_type_annotation(
#' seurat_object = seurat_obj,
#' method = "STCAT",
#' )
#'
#' # Both methods with custom parameters
#' seurat_obj <- cell_type_annotation(
#' seurat_object = seurat_obj,
#' method = "both",
#' reference = "hpca",
#' cluster_column = "wsnn_clusters",
#' singler_workers = 8
#' )
#' }
#'
#' @import Seurat
#' @import SingleR
#' @import celldex
#' @import SingleCellExperiment
#' @import BiocParallel
#' @importFrom reticulate py_config import py_run_string conda_python
#' @importFrom methods as
#' @export
cell_type_annotation <- function(
    seurat_object,
    method = c("both", "SingleR", "STCAT"),
    reference = "monaco",
    cluster_column = "RNA_clusters",
    assay = "RNA",
    singler_label_type = c("main", "fine"),
    singler_workers = 4,
    singler_fine_chunk_size = 20000L,
    stcat_python_script = NULL,
    stcat_batch_column = NULL,
    prefix = NULL,
    fallback_column = "SingleR_fine",
    verbose = TRUE,
    cleanup_temp = TRUE
) {
  
  if (verbose) message("=== Cell Type Annotation Pipeline ===\n")
  
  method <- match.arg(method)
  
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }
  
  if (!assay %in% names(seurat_object@assays)) {
    stop(sprintf("Assay '%s' not found in Seurat object", assay))
  }
  
  available_cols <- colnames(seurat_object@meta.data)
  if (!cluster_column %in% available_cols) {

    cluster_candidates <- grep("cluster", available_cols, ignore.case = TRUE, value = TRUE)
    error_msg <- sprintf("Cluster column '%s' not found in metadata.\nAvailable columns: %s",
                         cluster_column, paste(available_cols, collapse = ", "))
    if (length(cluster_candidates) > 0) {
      error_msg <- paste0(error_msg, sprintf("\nCluster-like columns detected: %s",
                                             paste(cluster_candidates, collapse = ", ")))
    }
    stop(error_msg)
  }
  
  if (!all(singler_label_type %in% c("main", "fine"))) {
    stop("singler_label_type must be 'main', 'fine', or both")
  }
  
  options(future.globals.maxSize = 8 * 1024^3)
  
  actual_fallback <- if (!is.null(prefix)) {
    paste0(prefix, "_", fallback_column)
  } else {
    fallback_column
  }
  
  gate_cols_wanted <- if (!is.null(prefix)) {
    paste0(prefix, "_", c("SingleR_main", "SingleR_fine"))
  } else {
    c("SingleR_main", "SingleR_fine")
  }
  has_gate <- any(gate_cols_wanted %in% colnames(seurat_object@meta.data))
  has_fallback <- actual_fallback %in% colnames(seurat_object@meta.data)
  if (method == "STCAT" && (!has_fallback || !has_gate)) {
    if (verbose) {
      if (!has_fallback)
        message(sprintf("\nNote: Fallback column '%s' not found in metadata.", actual_fallback))
      if (!has_gate)
        message("Note: No SingleR_main/SingleR_fine columns found for the T cell gate.")
      message("Automatically running SingleR annotation first...\n")
    }

    method <- "both"
  }
  
  ref_data <- NULL
  if (method %in% c("SingleR", "both")) {
    if (verbose) message("Loading reference dataset...")
    
    ref_data <- .load_reference_dataset(reference)
    
    if (verbose) {
      message(sprintf("Reference dataset loaded: %d cells, %d genes\n",
                      ncol(ref_data), nrow(ref_data)))
    }
  }
  
  if (method %in% c("SingleR", "both")) {
    if (verbose) message("\n--- Running SingleR Annotation ---")
    
    seurat_object <- .run_singler_annotation(
      seurat_object = seurat_object,
      ref_data = ref_data,
      cluster_column = cluster_column,
      assay = assay,
      label_types = singler_label_type,
      workers = singler_workers,
      prefix = prefix,
      verbose = verbose
    )
    
    if (verbose) message("SingleR annotation complete\n")

    if ("fine" %in% singler_label_type) {
      if (verbose) message("--- Refining fine labels per cell (chunked SingleR) ---")
      seurat_object <- .run_singler_fine_percell(
        seurat_object = seurat_object,
        ref_data = ref_data,
        assay = assay,
        prefix = prefix,
        chunk_size = singler_fine_chunk_size,
        verbose = verbose
      )
      if (verbose) message("Per-cell fine annotation complete\n")
    }

    seurat_object <- .rescue_missing_vocab(
      seurat_object = seurat_object,
      assay = assay,
      prefix = prefix,
      verbose = verbose
    )
  }

  if (method %in% c("STCAT", "both")) {
    if (verbose) message("\n--- Running STCAT Annotation ---")
    
    seurat_object <- .run_stcat_annotation(
      seurat_object = seurat_object,
      assay = assay,
      python_script = stcat_python_script,
      prefix = prefix,
      fallback_column = actual_fallback,
      batch_column = stcat_batch_column,
      cleanup = cleanup_temp,
      verbose = verbose
    )
    
    if (verbose) message("STCAT annotation complete\n")
  }
  
  if (method %in% c("STCAT", "both")) {
    seurat_object <- .harmonize_prediction_labels(seurat_object, prefix = prefix,
                                                  verbose = verbose)
  }

  seurat_object <- .drop_working_columns(seurat_object, prefix = prefix,
                                         cluster_column = cluster_column,
                                         verbose = verbose)

  if (verbose) {
    message("\n=== Annotation Summary ===")
    .print_annotation_summary(seurat_object, method, prefix)
  }
  
  return(seurat_object)
}

.load_reference_dataset <- function(reference) {

  if (inherits(reference, "SingleCellExperiment")) {
    return(reference)
  }

  loader <- switch(
    tolower(reference),
    "monaco" = celldex::MonacoImmuneData,
    "hpca" = celldex::HumanPrimaryCellAtlasData,
    "blueprint_encode" = celldex::BlueprintEncodeData,
    "dice" = celldex::DatabaseImmuneCellExpressionData,
    "immgen" = celldex::ImmGenData,
    stop(sprintf("Unknown reference: '%s'. Use 'monaco', 'hpca', 'blueprint_encode', 'dice', 'immgen', or provide custom SCE object", reference))
  )

  tryCatch(loader(), error = function(e) {
    stop(sprintf("Failed to load reference '%s' (celldex): %s",
                 reference, conditionMessage(e)), call. = FALSE)
  })
}

.run_singler_annotation <- function(
    seurat_object,
    ref_data,
    cluster_column,
    assay,
    label_types,
    workers,
    prefix,
    verbose
) {

  DefaultAssay(seurat_object) <- assay

  all_clusters <- as.character(seurat_object@meta.data[[cluster_column]])
  if (any(is.na(all_clusters))) {
    stop(sprintf(
      "Cluster column '%s' contains NA values. Ensure all cells have cluster assignments.",
      cluster_column
    ))
  }
  cluster_levels <- sort(unique(all_clusters))
  n_clusters <- length(cluster_levels)
  n_cells <- length(all_clusters)

  if (verbose) message(sprintf(
    " Building pseudo-bulk profiles: %d clusters from %d cells...",
    n_clusters, n_cells
  ))

  expr_mat <- .sf_fetch_layer(seurat_object, assay, "counts")
  expr_layer <- "counts"
  if (is.null(expr_mat)) {
    expr_mat <- .sf_fetch_layer(seurat_object, assay, "data")
    if (is.null(expr_mat)) stop(sprintf(
      "Could not retrieve 'counts' or 'data' layer from assay '%s'.", assay))
    expr_layer <- "data"
  }

  assign_mat <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = match(all_clusters, cluster_levels),
    x = 1,
    dims = c(n_cells, n_clusters),
    dimnames = list(colnames(expr_mat), cluster_levels)
  )

  gene_names <- rownames(expr_mat)

  pseudo_bulk <- tryCatch({
    pb <- expr_mat %*% assign_mat
    as.matrix(pb)
  }, error = function(e) {
    stop(sprintf("Pseudo-bulk aggregation failed: %s", e$message))
  })
  rm(expr_mat, assign_mat); gc()
  dimnames(pseudo_bulk) <- list(gene_names, cluster_levels)

  if (expr_layer == "counts") {

    col_tot <- colSums(pseudo_bulk)
    col_tot[col_tot == 0] <- 1
    pseudo_bulk <- log1p(sweep(pseudo_bulk, 2, col_tot, "/") * 1e6)
  } else {

    cluster_sizes <- table(all_clusters)[cluster_levels]
    pseudo_bulk <- sweep(pseudo_bulk, 2, as.numeric(cluster_sizes), "/")
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = pseudo_bulk)
  )
  rm(pseudo_bulk); gc()

  common_genes <- intersect(rownames(sce), rownames(ref_data))
  if (length(common_genes) == 0) {
    stop(sprintf(paste0(
      "No common genes between test data (%d genes) and reference (%d genes).\n",
      " Test gene examples: %s\n",
      " Reference gene examples: %s\n",
      " This usually means the Seurat object uses Ensembl IDs while the\n",
      " reference uses gene symbols (or vice versa). Check load_data gene ID conversion."
    ),
    nrow(sce), nrow(ref_data),
    paste(head(rownames(sce), 5), collapse = ", "),
    paste(head(rownames(ref_data), 5), collapse = ", ")
    ))
  }
  if (verbose) message(sprintf(" %d common genes between test and reference", length(common_genes)))

  bpparam <- SerialParam()

  for (label_type in label_types) {

    if (verbose) message(sprintf("Annotating with %s labels...", label_type))

    label_col <- switch(
      label_type,
      "main" = "label.main",
      "fine" = "label.fine",
      stop(sprintf("Invalid label type: '%s'. Use 'main' or 'fine'.", label_type))
    )

    if (!label_col %in% colnames(colData(ref_data))) {
      warning(sprintf(
        "Label column '%s' not found in reference data. Skipping %s annotation.",
        label_col, label_type
      ))
      next
    }

    pred <- SingleR(
      test = sce,
      ref = ref_data,
      labels = colData(ref_data)[[label_col]],
      BPPARAM = bpparam

    )

    col_label <- if (!is.null(prefix)) {
      sprintf("%s_SingleR_%s", prefix, label_type)
    } else {
      sprintf("SingleR_%s", label_type)
    }

    cluster_label_map <- setNames(as.character(pred$labels), rownames(pred))
    cell_labels <- unname(cluster_label_map[all_clusters])
    seurat_object@meta.data[[col_label]] <- cell_labels

    scores_mat <- pred$scores
    top_score <- apply(scores_mat, 1, max, na.rm = TRUE)
    delta <- apply(scores_mat, 1, function(r) {
      r <- sort(r[is.finite(r)], decreasing = TRUE)
      if (length(r) >= 2) r[1] - r[2] else NA_real_
    })
    names(top_score) <- rownames(pred)
    names(delta) <- rownames(pred)

    score_col <- if (!is.null(prefix)) sprintf("%s_SingleR_%s_score", prefix, label_type)
                 else sprintf("SingleR_%s_score", label_type)
    delta_col <- if (!is.null(prefix)) sprintf("%s_SingleR_%s_delta", prefix, label_type)
                 else sprintf("SingleR_%s_delta", label_type)
    seurat_object@meta.data[[score_col]] <- unname(top_score[all_clusters])
    seurat_object@meta.data[[delta_col]] <- unname(delta[all_clusters])

    if (verbose) message(sprintf(
      " Identified %d unique %s cell types across %d clusters (median score %.3f)",
      length(unique(cell_labels)), label_type, n_clusters,
      stats::median(top_score, na.rm = TRUE)
    ))
  }

  return(seurat_object)
}

.run_singler_fine_percell <- function(
    seurat_object,
    ref_data,
    assay,
    prefix,
    chunk_size = 20000L,
    verbose = TRUE
) {

  if (!"label.fine" %in% colnames(colData(ref_data))) {
    warning("Reference has no 'label.fine' column. Skipping per-cell fine refinement.")
    return(seurat_object)
  }

  col_fine <- if (!is.null(prefix)) sprintf("%s_SingleR_fine", prefix) else "SingleR_fine"
  col_clu <- paste0(col_fine, "_cluster")

  for (suf in c("", "_score", "_delta")) {
    src <- paste0(col_fine, suf)
    if (src %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[[paste0(col_clu, suf)]] <- seurat_object@meta.data[[src]]
    }
  }

  expr <- .sf_fetch_layer(seurat_object, assay, "data")
  from_counts <- is.null(expr)
  if (from_counts) {
    expr <- .sf_fetch_layer(seurat_object, assay, "counts")
    if (is.null(expr)) stop(
      "Neither 'data' nor 'counts' layer available for per-cell fine annotation.")
    if (verbose) message(" 'data' layer unavailable - computing log1p(CP10K) from counts per chunk")
  }

  n_cells <- ncol(expr)
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq(1L, n_cells, by = chunk_size)
  n_chunks <- length(starts)

  if (verbose) message(sprintf(
    " Annotating %d cells per cell in %d chunk(s) of up to %d cells...",
    n_cells, n_chunks, chunk_size))

  labels <- character(n_cells)
  scores <- rep(NA_real_, n_cells)
  deltas <- rep(NA_real_, n_cells)

  bpparam <- SerialParam()
  ref_labels <- colData(ref_data)[["label.fine"]]

  for (i in seq_along(starts)) {
    idx <- starts[i]:min(starts[i] + chunk_size - 1L, n_cells)
    chunk <- expr[, idx, drop = FALSE]

    chunk <- as(chunk, "dgCMatrix")
    if (from_counts) {
      libs <- Matrix::colSums(chunk)
      libs[libs == 0] <- 1
      chunk <- log1p(Matrix::t(Matrix::t(chunk) / libs) * 1e4)
    }

    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = chunk))
    pred <- SingleR(test = sce, ref = ref_data, labels = ref_labels, BPPARAM = bpparam)

    labels[idx] <- as.character(pred$labels)
    sm <- pred$scores
    scores[idx] <- unname(apply(sm, 1, max, na.rm = TRUE))
    deltas[idx] <- unname(apply(sm, 1, function(r) {
      r <- sort(r[is.finite(r)], decreasing = TRUE)
      if (length(r) >= 2) r[1] - r[2] else NA_real_
    }))

    if (verbose && n_chunks > 1) {
      message(sprintf(" chunk %d/%d done (%d cells)", i, n_chunks, length(idx)))
    }
    rm(chunk, sce, pred); gc(FALSE)
  }

  seurat_object@meta.data[[col_fine]] <- labels
  seurat_object@meta.data[[paste0(col_fine, "_score")]] <- scores
  seurat_object@meta.data[[paste0(col_fine, "_delta")]] <- deltas

  if (verbose) {
    n_changed <- if (col_clu %in% colnames(seurat_object@meta.data)) {
      sum(labels != as.character(seurat_object@meta.data[[col_clu]]), na.rm = TRUE)
    } else NA_integer_
    message(sprintf(
      " Per-cell fine annotation: %d unique labels (median score %.3f)%s",
      length(unique(labels)), stats::median(scores, na.rm = TRUE),
      if (!is.na(n_changed)) sprintf(" - %d cells (%.1f%%) differ from their cluster-level label",
                                     n_changed, 100 * n_changed / n_cells) else ""))
  }

  return(seurat_object)
}

.is_tcell_label <- function(labels) {
  lab <- as.character(labels)
  res <- rep(FALSE, length(lab))
  ok <- which(!is.na(lab))
  if (length(ok) == 0) return(res)

  pattern <- paste(
    "\\bT[ _-]?cells?\\b",
    "\\bT[ _-]?lymphocytes?\\b",
    "\\bMAIT\\b",
    "\\bNKT\\b",
    "\\bT[ _-]?regs?\\b",
    "T regulatory",
    "\\bTh[0-9]{1,2}\\b",
    "\\bTfh\\b",
    "Follicular helper",
    "gamma[ -]?delta",
    sep = "|"
  )
  res[ok] <- grepl(pattern, lab[ok], ignore.case = TRUE, perl = TRUE)
  res
}

.prediction_harmonization_map <- function() {
  c(
    "MAIT cells" = "MAIT",
    "Naive CD4 T cells" = "CD4 Tn",
    "Naive CD8 T cells" = "CD8 Tn",
    "Central memory CD8 T cells" = "CD8 Tcm",
    "Effector memory CD8 T cells" = "CD8 Tem",
    "Terminal effector CD8 T cells" = "CD8 Temra",
    "Terminal effector CD4 T cells" = "CD4 Temra",
    "T regulatory cells" = "CD4 Treg",
    "Th1 cells" = "CD4 Th1",
    "Th17 cells" = "CD4 Th17",
    "Vd2 gd T cells" = "Tgd",
    "Non-Vd2 gd T cells" = "Tgd"
  )
}

.harmonize_prediction_labels <- function(seurat_object, prefix = NULL, verbose = TRUE) {
  col_prediction <- if (!is.null(prefix)) paste0(prefix, "_Prediction") else "Prediction"
  if (!col_prediction %in% colnames(seurat_object@meta.data)) return(seurat_object)

  map <- .prediction_harmonization_map()

  pred <- as.character(seurat_object@meta.data[[col_prediction]])
  hit <- !is.na(pred) & pred %in% names(map)
  if (any(hit)) {
    pred[hit] <- unname(map[pred[hit]])
    seurat_object@meta.data[[col_prediction]] <- pred
    if (verbose) {
      counts <- table(pred[hit])
      message(sprintf(
        "Harmonized %d Prediction labels to the STCAT vocabulary: %s",
        sum(hit),
        paste(sprintf("%s->%s: %d", names(counts), unname(map[names(counts)]), counts),
              collapse = ", ")))
    }
  } else if (verbose) {
    message("Prediction label harmonization: nothing to rename")
  }
  return(seurat_object)
}

.drop_working_columns <- function(seurat_object, prefix = NULL,
                                  cluster_column = NULL, verbose = TRUE) {
  drop_base <- c("SingleR_fine_cluster", "SingleR_fine_cluster_score",
                 "SingleR_fine_cluster_delta", "Cluster", "vocab_rescue_score",
                 "SingleR_main_delta", "SingleR_fine_delta")
  drop_cols <- c(
    if (!is.null(prefix)) paste0(prefix, "_", drop_base) else drop_base,
    "seurat_clusters"
  )

  drop_cols <- setdiff(drop_cols, cluster_column)
  drop_cols <- intersect(drop_cols, colnames(seurat_object@meta.data))
  if (length(drop_cols) > 0) {
    for (dc in drop_cols) seurat_object@meta.data[[dc]] <- NULL
    if (verbose) message(sprintf(
      "Removed working metadata columns: %s", paste(drop_cols, collapse = ", ")))
  }
  return(seurat_object)
}

.vocab_rescue_signatures <- function() {
  list(
    "Platelets" = c("PPBP", "PF4", "NRGN", "GNG11", "TUBB1", "CAVIN2",
                       "GP9", "CLU", "RGS18", "HIST1H2AC"),
    "Erythrocytes" = c("HBB", "HBA1", "HBA2", "HBM", "HBD", "ALAS2",
                       "AHSP", "CA1", "SLC4A1", "SELENBP1")
  )
}

.rescue_missing_vocab <- function(
    seurat_object,
    assay,
    prefix = NULL,
    signatures = .vocab_rescue_signatures(),
    min_z = 2,
    min_margin = 0.5,
    verbose = TRUE
) {

  col_main <- if (!is.null(prefix)) paste0(prefix, "_SingleR_main") else "SingleR_main"
  col_fine <- if (!is.null(prefix)) paste0(prefix, "_SingleR_fine") else "SingleR_fine"
  target_cols <- intersect(c(col_main, col_fine), colnames(seurat_object@meta.data))
  if (length(target_cols) == 0) {
    if (verbose) message(" No SingleR columns present - vocabulary rescue skipped")
    return(seurat_object)
  }

  expr <- .sf_fetch_layer(seurat_object, assay, "data")
  from_counts <- is.null(expr)
  if (from_counts) {
    expr <- .sf_fetch_layer(seurat_object, assay, "counts")
    if (is.null(expr)) {
      warning("Vocabulary rescue skipped: no 'data' or 'counts' layer available.")
      return(seurat_object)
    }
  }

  genes_all <- intersect(unique(unlist(signatures)), rownames(expr))
  usable <- Filter(function(gs) length(intersect(gs, genes_all)) >= 3, signatures)
  if (length(usable) == 0) {
    if (verbose) message(" Vocabulary rescue: no signature has >=3 genes in this assay - skipped")
    return(seurat_object)
  }

  libs <- NULL
  if (from_counts) {
    libs <- colSums(expr)
    libs[libs == 0] <- 1
  }

  sub <- as(expr[genes_all, , drop = FALSE], "dgCMatrix")
  if (from_counts) {
    sub <- log1p(Matrix::t(Matrix::t(sub) / libs) * 1e4)
  }

  mu <- Matrix::rowMeans(sub)
  sd_ <- sqrt(pmax(Matrix::rowMeans(sub^2) - mu^2, 0)); sd_[sd_ == 0] <- 1

  score <- sapply(names(usable), function(ct) {
    gs <- intersect(usable[[ct]], genes_all)
    as.numeric(Matrix::colSums(sub[gs, , drop = FALSE] / sd_[gs]) / length(gs) -
               sum(mu[gs] / sd_[gs]) / length(gs))
  })
  if (is.null(dim(score))) score <- matrix(score, ncol = length(usable),
                                           dimnames = list(NULL, names(usable)))

  top_idx <- max.col(score, ties.method = "first")
  top_score <- score[cbind(seq_len(nrow(score)), top_idx)]
  second <- if (ncol(score) >= 2) {
    apply(score, 1, function(r) sort(r, decreasing = TRUE)[2])
  } else rep(-Inf, nrow(score))

  hit <- top_score >= min_z & (top_score - second) >= min_margin
  rescue_label <- ifelse(hit, colnames(score)[top_idx], NA_character_)

  rc <- if (!is.null(prefix)) paste0(prefix, "_vocab_rescue") else "vocab_rescue"
  seurat_object@meta.data[[rc]] <- rescue_label
  seurat_object@meta.data[[paste0(rc, "_score")]] <- ifelse(hit, round(top_score, 4), NA_real_)

  if (any(hit)) {
    for (cc in target_cols) {
      seurat_object@meta.data[[cc]][hit] <- rescue_label[hit]
    }
    if (verbose) {
      counts <- table(rescue_label[hit])
      message(sprintf(
        " Vocabulary rescue: %d cells relabelled (%s) - recorded in %s",
        sum(hit),
        paste(sprintf("%s: %d", names(counts), counts), collapse = ", "),
        rc))
    }
  } else if (verbose) {
    message(" Vocabulary rescue: no cells met the evidence threshold")
  }

  return(seurat_object)
}

.compute_tcell_subset_embedding <- function(
    seurat_object,
    tcell_barcodes,
    assay,
    batch_column = NULL,
    n_hvg = 1500L,
    dims = 30L,
    harmony_theta = 2,
    harmony_niter = 10,
    verbose = TRUE
) {
  tryCatch({
    if (verbose) message(sprintf(
      " Re-embedding the %d-cell T subset (HVG -> PCA%s)...",
      length(tcell_barcodes),
      if (!is.null(batch_column)) " -> Harmony if multi-batch" else ""))

    sub <- subset(seurat_object, cells = tcell_barcodes)
    DefaultAssay(sub) <- assay

    has_data <- !is.null(.sf_fetch_layer(sub, assay, "data"))
    if (!has_data) sub <- NormalizeData(sub, verbose = FALSE)

    sub <- FindVariableFeatures(sub, nfeatures = n_hvg, verbose = FALSE)
    hvg <- VariableFeatures(sub)
    if (length(hvg) < 10) stop("fewer than 10 variable features on the T subset")

    sub <- ScaleData(sub, features = hvg, verbose = FALSE)
    npcs <- max(2L, min(as.integer(dims), length(tcell_barcodes) - 1L, length(hvg) - 1L))
    sub <- RunPCA(sub, features = hvg, npcs = npcs, verbose = FALSE)
    red <- "pca"

    if (!is.null(batch_column) && batch_column %in% colnames(sub@meta.data)) {
      n_batches <- length(unique(as.character(sub@meta.data[[batch_column]])))
      if (n_batches >= 2) {
        if (requireNamespace("harmony", quietly = TRUE)) {

          rh <- tryCatch(utils::getS3method("RunHarmony", "Seurat"),
                           error = function(e) harmony::RunHarmony)
          fmls <- names(formals(rh))
          red_name <- if ("reduction.use" %in% fmls) "reduction.use" else "reduction"
          iter_name <- if (utils::packageVersion("harmony") >= "1.0.0") "max_iter"
                       else "max.iter.harmony"
          rh_args <- list(
            object = sub,
            group.by.vars = batch_column,
            reduction.save = "harmony.tsub",
            theta = harmony_theta,
            verbose = FALSE)
          rh_args[[red_name]] <- "pca"
          rh_args[[iter_name]] <- harmony_niter
          if ("project.dim" %in% fmls) rh_args$project.dim <- FALSE
          sub <- do.call(harmony::RunHarmony, rh_args)
          red <- "harmony.tsub"
          if (verbose) message(sprintf(
            " T-subset embedding Harmony-corrected on '%s' (%d batches)",
            batch_column, n_batches))
        } else {
          warning(sprintf(paste0(
            "Batch column '%s' has %d levels among the gated T cells but the ",
            "'harmony' package is not installed - the T-subset embedding sent ",
            "to STCAT is NOT batch-corrected."), batch_column, n_batches))
        }
      } else if (verbose) {
        message(sprintf(" Single batch ('%s') among gated cells - Harmony not needed",
                        batch_column))
      }
    }

    Embeddings(sub, red)
  }, error = function(e) {
    warning(sprintf(
      "T-subset re-embedding failed (%s); falling back to the full-data embedding.",
      conditionMessage(e)))
    NULL
  })
}

.ensure_python <- function(verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(NULL))

  rp <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(rp) && file.exists(rp)) {
    if (verbose) message(" Python interpreter: ", rp, " (RETICULATE_PYTHON)")
    return(invisible(rp))
  }

  env_name <- Sys.getenv("PIPELINE_PY_ENV", "")
  if (!nzchar(env_name)) {
    if (verbose)
      message(" RETICULATE_PYTHON / PIPELINE_PY_ENV unset - using reticulate's default interpreter.")
    return(invisible(NULL))
  }

  conda_path <- Sys.which("conda")
  if (!nzchar(conda_path))
    conda_path <- Find(file.exists, path.expand(c(
      "/opt/conda/bin/conda",
      "~/miniforge3/bin/conda", "~/mambaforge/bin/conda",
      "~/miniconda3/bin/conda", "~/anaconda3/bin/conda")))

  py <- NULL
  if (!is.null(conda_path) && nzchar(conda_path)) {
    py <- tryCatch(reticulate::conda_python(env_name, conda = conda_path),
                   error = function(e) NULL)
    if (is.null(py) || !file.exists(py)) {
      cand <- file.path(dirname(dirname(conda_path)), "envs", env_name, "bin", "python")
      py <- if (file.exists(cand)) cand else NULL
    }
  }

  if (!is.null(py) && file.exists(py)) {
    Sys.setenv(RETICULATE_PYTHON = py)
    if (verbose) message(" Python interpreter: ", py, " (conda env '", env_name, "')")
    return(invisible(py))
  }

  if (verbose)
    message(" Could not resolve conda env '", env_name,
            "' - using reticulate's default interpreter.")
  invisible(NULL)
}
.run_stcat_annotation <- function(
    seurat_object,
    assay,
    python_script,
    prefix,
    fallback_column,
    batch_column = NULL,
    cleanup = TRUE,
    verbose = TRUE
) {

  meta_cols <- colnames(seurat_object@meta.data)
  gate_cols <- if (!is.null(prefix)) {
    paste0(prefix, "_", c("SingleR_main", "SingleR_fine"))
  } else {
    c("SingleR_main", "SingleR_fine")
  }
  gate_cols <- intersect(gate_cols, meta_cols)
  if (length(gate_cols) == 0) {
    stop(paste0(
      "T cell gating requires SingleR annotations before STCAT, but neither ",
      "SingleR_main nor SingleR_fine is present in the metadata. ",
      "Run with method = 'both' (SingleR first, then STCAT)."))
  }

  t_mask <- Reduce(`|`, lapply(gate_cols, function(cc) {
    .is_tcell_label(seurat_object@meta.data[[cc]])
  }))
  tcell_barcodes <- colnames(seurat_object)[t_mask]

  if (verbose) {
    message(sprintf(
      "T cell gate: %d / %d cells (%.1f%%) carry T cell-referencing labels in [%s] - only these are sent to STCAT",
      length(tcell_barcodes), ncol(seurat_object),
      100 * length(tcell_barcodes) / ncol(seurat_object),
      paste(gate_cols, collapse = ", ")))
  }

  col_prediction <- if (!is.null(prefix)) paste0(prefix, "_Prediction") else "Prediction"
  col_cluster <- if (!is.null(prefix)) paste0(prefix, "_Cluster") else "Cluster"
  col_uncertainty <- if (!is.null(prefix)) paste0(prefix, "_Uncertainty_score") else "Uncertainty_score"

  if (length(tcell_barcodes) == 0) {
    warning("No cells passed the T cell gate; skipping STCAT entirely.")
    seurat_object@meta.data[[col_prediction]] <- NA_character_
    seurat_object@meta.data[[col_cluster]] <- NA_character_
    seurat_object@meta.data[[col_uncertainty]] <- NA_real_
    if (!is.null(fallback_column) && fallback_column %in% meta_cols) {

      fb_labels <- as.character(seurat_object@meta.data[[fallback_column]])
      seurat_object@meta.data[[col_prediction]] <- fb_labels
      seurat_object@meta.data[[col_cluster]] <- fb_labels
      if (verbose) message(sprintf(
        " All cells received the fallback label (%s) in %s",
        fallback_column, col_prediction))
    } else if (!is.null(fallback_column)) {
      warning(sprintf("Fallback column '%s' not found; %s left NA.",
                      fallback_column, col_prediction))
    }
    return(seurat_object)
  }

  temp_dir <- tempfile(pattern = "stcat_")
  dir.create(temp_dir, recursive = TRUE)
  
  if (verbose) message(sprintf("Created temporary directory: %s", temp_dir))
  
  if (cleanup) {
    on.exit({
      if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
        if (verbose) message("Cleaned up temporary directory")
      }
    }, add = TRUE)
  }
  
  if (verbose) message("Configuring Python environment...")
  
  tryCatch({

    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("Package 'reticulate' is required for STCAT. Please install it.")
    }
    
    .ensure_python(verbose)

    reticulate::py_config()

    py_exe <- reticulate::py_config()$python
    probe <- suppressWarnings(system2(
      py_exe, c("-c", shQuote("import scanpy, anndata, numpy, STCAT")),
      stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(probe, "status")) && attr(probe, "status") != 0) {
      stop(sprintf(paste0(
        "The Python interpreter cannot import the modules STCAT needs.\n",
        " interpreter: %s\n error: %s\n",
        "The image sets RETICULATE_PYTHON to the bundled environment; if this is ",
        "running outside the image, set RETICULATE_PYTHON or PIPELINE_PY_ENV."),
        py_exe, paste(probe, collapse = " | ")))
    }

    if (verbose) message("Python environment configured successfully (", py_exe, ")")
    
  }, error = function(e) {
    stop(sprintf("Failed to configure Python environment: %s", e$message))
  })
  
  if (verbose) message("Converting Seurat object to AnnData...")
  
  DefaultAssay(seurat_object) <- assay
  
  counts_raw <- tryCatch(
    GetAssayData(seurat_object, layer = "counts"),
    error = function(e) tryCatch(
      GetAssayData(seurat_object, slot = "counts"),
      error = function(e2) stop(sprintf("Failed to extract counts: %s", e2$message))
    )
  )

  counts_raw <- counts_raw[, tcell_barcodes, drop = FALSE]

  genes_to_use <- rownames(counts_raw)
  if (verbose) message(sprintf(
    " Sending all %d genes to STCAT (full gene set - needed for marker coverage).",
    length(genes_to_use)))

  input_h5ad <- file.path(temp_dir, "input_data.h5ad")
  output_h5ad <- file.path(temp_dir, "stcat_results.h5ad")

  if (verbose) message(sprintf(
    " Writing %d genes x %d cells to h5ad via BPCells (disk-to-disk, no RAM spike)...",
    length(genes_to_use), ncol(counts_raw)
  ))

  if (inherits(counts_raw, "IterableMatrix")) {
    BPCells::write_matrix_anndata_hdf5(counts_raw, input_h5ad, group = "X")
  } else {
    if (verbose) message(" Counts matrix is not BPCells IterableMatrix - writing via anndata")

    if (!inherits(counts_raw, "dgCMatrix"))
      counts_raw <- as(counts_raw, "dgCMatrix")
    counts_t <- Matrix::t(counts_raw)

    ad <- reticulate::import("anndata", convert = FALSE)
    scipy <- reticulate::import("scipy.sparse", convert = FALSE)

    mat_py <- scipy$csc_matrix(
      reticulate::tuple(
        reticulate::np_array(counts_t@x, dtype = "float32"),
        reticulate::np_array(counts_t@i, dtype = "int32"),
        reticulate::np_array(counts_t@p, dtype = "int32")
      ),
      shape = reticulate::tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
    )
    mat_csr <- mat_py$tocsr()

    adata <- ad$AnnData(X = mat_csr)
    adata$var_names <- reticulate::r_to_py(colnames(counts_t))
    adata$obs_names <- reticulate::r_to_py(rownames(counts_t))
    adata$write_h5ad(input_h5ad)
    rm(adata, mat_py, mat_csr, counts_t); gc()
  }
  rm(counts_raw); gc()

  if (verbose) {
    message(sprintf(" Dimensions: %d gated T cells x %d genes (of %d cells total)",
                    length(tcell_barcodes), length(genes_to_use), ncol(seurat_object)))
  }

  emb_csv <- NULL
  emb_max_dims <- 30L
  emb <- .compute_tcell_subset_embedding(
    seurat_object = seurat_object,
    tcell_barcodes = tcell_barcodes,
    assay = assay,
    batch_column = batch_column,
    dims = emb_max_dims,
    verbose = verbose
  )
  if (!is.null(emb)) {
    if (ncol(emb) > emb_max_dims) emb <- emb[, seq_len(emb_max_dims), drop = FALSE]
    emb_csv <- file.path(temp_dir, "embedding.csv")
    write.csv(emb, emb_csv, row.names = TRUE)
    if (verbose) message(sprintf(
      " Saved T-subset embedding (%d dims, %d cells) for STCAT neighbour graph",
      ncol(emb), nrow(emb)))
  } else {
    for (.red in c("harmony.full", "harmony", "integrated.rpca.full", "integrated.dr", "pca", "rna.pca")) {
      if (.red %in% names(seurat_object@reductions)) {
        emb <- Embeddings(seurat_object, .red)

        emb_cells <- intersect(rownames(emb), tcell_barcodes)
        emb <- emb[emb_cells, , drop = FALSE]
        if (ncol(emb) > emb_max_dims) emb <- emb[, seq_len(emb_max_dims), drop = FALSE]
        emb_csv <- file.path(temp_dir, "embedding.csv")
        write.csv(emb, emb_csv, row.names = TRUE)
        if (verbose) message(sprintf(
          " Saved Seurat '%s' embedding (%d dims) for STCAT neighbour graph (fallback)",
          .red, ncol(emb)
        ))
        break
      }
    }
  }
  if (is.null(emb_csv) && verbose) {
    message(" No embedding available - STCAT will compute its own neighbour graph")
  }

  if (verbose) message(sprintf("h5ad written to: %s", input_h5ad))
  
  if (verbose) message("Running STCAT analysis...")
  
  if (is.null(python_script)) {
    python_script <- "python/run_stcat.py"
  }
  
  if (!file.exists(python_script)) {
    stop(sprintf("STCAT Python script not found: %s\n Please provide the path via stcat_python_script parameter.",
                 python_script))
  }
  
  python_exe <- reticulate::py_config()$python
  
  if (verbose) message(sprintf("Using Python: %s", python_exe))
  
  result <- system2(
    python_exe,
    args = c(python_script, input_h5ad, output_h5ad),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop(sprintf("STCAT Python script failed:\n%s", paste(result, collapse = "\n")))
  }
  
  if (verbose && length(result) > 0) {
    message("STCAT output:")
    message(paste(result, collapse = "\n"))
  }
  
  possible_outputs <- c(
    output_h5ad,
    file.path(temp_dir, "stcat_results.h5ad"),
    file.path(temp_dir, "output_data.h5ad"),
    file.path(dirname(input_h5ad), basename(output_h5ad))
  )
  
  output_file <- NULL
  for (possible_file in possible_outputs) {
    if (file.exists(possible_file)) {
      output_file <- possible_file
      break
    }
  }
  
  if (is.null(output_file)) {
    h5ad_files <- list.files(temp_dir, pattern = "\\.h5ad$", full.names = TRUE)

    h5ad_files <- setdiff(h5ad_files, input_h5ad)
    
    if (length(h5ad_files) > 0) {
      output_file <- h5ad_files[1]
      if (verbose) message(sprintf("Found output file: %s", output_file))
    }
  }
  
  if (is.null(output_file) || !file.exists(output_file)) {

    all_files <- list.files(temp_dir, full.names = TRUE)
    stop(sprintf("STCAT output file not created. Expected: %s\nFiles in temp directory:\n%s",
                 output_h5ad, paste(all_files, collapse = "\n")))
  }
  
  if (verbose) message("Reading STCAT results...")
  
  stcat_meta <- tryCatch({
    if (requireNamespace("anndata", quietly = TRUE)) {
      adata_result <- anndata::read_h5ad(output_file)
      as.data.frame(adata_result$obs)
    } else {

      ad <- reticulate::import("anndata")
      adata_result <- ad$read_h5ad(output_file)
      as.data.frame(adata_result$obs)
    }
  }, error = function(e) {
    stop(sprintf("Failed to read STCAT results: %s\n Make sure anndata package is installed.", e$message))
  })
  
  if (is.null(rownames(stcat_meta))) {
    stop("STCAT results missing cell barcodes (rownames)")
  }
  
  required_cols <- c("Prediction", "Cluster", "Uncertainty score")
  missing_cols <- setdiff(required_cols, colnames(stcat_meta))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in STCAT results: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  
  stcat_meta$Prediction <- as.character(stcat_meta$Prediction)
  stcat_meta$Cluster <- as.character(stcat_meta$Cluster)
  
  if (verbose) message("Merging STCAT results into Seurat object...")
  
  common_cells <- intersect(colnames(seurat_object), rownames(stcat_meta))
  
  if (length(common_cells) == 0) {
    stop("No overlapping cells between Seurat object and STCAT results")
  }
  
  if (verbose) {
    message(sprintf(" Matched %d / %d gated T cells (%d cells total in object)",
                    length(common_cells), length(tcell_barcodes), ncol(seurat_object)))
  }

  seurat_object@meta.data[[col_prediction]] <- NA
  seurat_object@meta.data[[col_cluster]] <- NA
  seurat_object@meta.data[[col_uncertainty]] <- NA
  
  stcat_meta_matched <- stcat_meta[common_cells, , drop = FALSE]
  
  seurat_object@meta.data[[col_prediction]][match(common_cells, colnames(seurat_object))] <- 
    stcat_meta_matched$Prediction
  seurat_object@meta.data[[col_cluster]][match(common_cells, colnames(seurat_object))] <- 
    stcat_meta_matched$Cluster
  seurat_object@meta.data[[col_uncertainty]][match(common_cells, colnames(seurat_object))] <- 
    stcat_meta_matched$`Uncertainty score`
  
  if (!is.null(fallback_column) && fallback_column %in% colnames(seurat_object@meta.data)) {

    fb_labels <- as.character(seurat_object@meta.data[[fallback_column]])

    non_common_cells <- setdiff(colnames(seurat_object), common_cells)
    if (length(non_common_cells) > 0) {
      non_common_idx <- match(non_common_cells, colnames(seurat_object))
      seurat_object@meta.data[[col_prediction]][non_common_idx] <- fb_labels[non_common_idx]
      seurat_object@meta.data[[col_cluster]][non_common_idx] <- fb_labels[non_common_idx]
      seurat_object@meta.data[[col_uncertainty]][non_common_idx] <- NA

      if (verbose) {
        n_gated_out <- length(setdiff(non_common_cells, tcell_barcodes))
        n_unmatched <- length(non_common_cells) - n_gated_out
        message(sprintf(
          " Applied fallback (%s) for %d cells: %d outside the T cell gate, %d gated but not returned by STCAT",
          fallback_column, length(non_common_cells), n_gated_out, n_unmatched))
      }
    }
    
    none_t_cells <- which(seurat_object@meta.data[[col_prediction]] == "None T")
    if (length(none_t_cells) > 0) {
      seurat_object@meta.data[[col_prediction]][none_t_cells] <- fb_labels[none_t_cells]
      seurat_object@meta.data[[col_cluster]][none_t_cells] <- fb_labels[none_t_cells]

      seurat_object@meta.data[[col_uncertainty]][none_t_cells] <- NA

      if (verbose) {
        message(sprintf(" Replaced %d 'None T' predictions with fallback (%s); their Uncertainty_score set to NA",
                        length(none_t_cells), fallback_column))
      }
    }
    
  } else if (!is.null(fallback_column)) {
    warning(sprintf("Fallback column '%s' not found. Cells with 'None T' and non-matching cells will retain original/NA predictions.",
                    fallback_column))
  }
  
  return(seurat_object)
}

.print_annotation_summary <- function(seurat_object, method, prefix) {
  
  meta <- seurat_object@meta.data
  
  if (method %in% c("SingleR", "both")) {
    
    col_main <- if (!is.null(prefix)) {
      paste0(prefix, "_SingleR_main")
    } else {
      "SingleR_main"
    }
    
    col_fine <- if (!is.null(prefix)) {
      paste0(prefix, "_SingleR_fine")
    } else {
      "SingleR_fine"
    }
    
    if (col_main %in% colnames(meta)) {
      n_main <- length(unique(meta[[col_main]]))
      message(sprintf("SingleR main cell types: %d", n_main))
    }
    
    if (col_fine %in% colnames(meta)) {
      n_fine <- length(unique(meta[[col_fine]]))
      message(sprintf("SingleR fine cell types: %d", n_fine))
    }
  }
  
  if (method %in% c("STCAT", "both")) {
    
    col_pred <- if (!is.null(prefix)) {
      paste0(prefix, "_Prediction")
    } else {
      "Prediction"
    }
    
    if (col_pred %in% colnames(meta)) {
      n_pred <- length(unique(meta[[col_pred]][!is.na(meta[[col_pred]])]))
      n_na <- sum(is.na(meta[[col_pred]]))
      message(sprintf("STCAT predictions: %d unique types (%d NAs)", n_pred, n_na))
    }
  }
  
  message("===========================\n")
}
