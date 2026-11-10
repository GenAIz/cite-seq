#' Aggregate Seurat Object for Machine Learning
#'
#' Aggregates single-cell RNA-seq data from a Seurat object by grouping cells
#' based on participant, timepoint, and cell type annotations. Returns both
#' count data (number of cells expressing each feature) and mean expression
#' values, suitable for downstream machine learning analyses.
#'
#' @param seurat_obj A Seurat object containing single-cell data with metadata
#'   and one or more assays (e.g., RNA, ADT, SCT).
#' @param participant_col Character string specifying the metadata column name
#'   for participant/sample IDs (e.g., "Participant ID").
#' @param timepoint_col Character string specifying the metadata column name
#'   for timepoint information (e.g., "Timepoint (months)").
#' @param celltype_col Character vector specifying one or more metadata column
#'   names for cell type annotations. If multiple columns are provided, each
#'   will be processed independently (e.g., c("Prediction", "SingleR_fine")).
#' @param columns_to_preserve Character vector of additional metadata columns
#'   to include in the output (e.g., c("Unique Sample ID", "batch")). Default
#'   is NULL.
#' @param assays Character vector specifying which assays to aggregate (e.g.,
#'   c("RNA", "ADT", "SCT")). Default is c("RNA", "ADT").
#' @param output_xlsx Character string specifying the path for Excel output file.
#'   If NULL (default), no file is written. When provided, results are saved
#'   as separate sheets in the workbook.
#'
#' @return
#' \strong{Single celltype column:} Returns a list with two elements:
#' \itemize{
#'   \item \code{counts}: data.table with cell counts expressing each feature
#'   \item \code{expression}: data.table with mean expression values
#' }
#'
#' \strong{Multiple celltype columns:} Returns a named list with results for each
#' celltype column:
#' \itemize{
#'   \item \code{counts_<celltype>}: counts for each celltype column
#'   \item \code{expression_<celltype>}: expression for each celltype column
#' }
#'
#' Both count and expression tables include:
#' \itemize{
#'   \item Grouping columns (participant, timepoint, celltype)
#'   \item Preserved metadata columns
#'   \item \code{cell_percentage}: percentage of each cell type within participant-timepoint
#'   \item \code{rna_ncells}: total cells with RNA expression > 0
#'   \item \code{adt_ncells}: total cells with ADT expression > 0 (if ADT assay present)
#'   \item Feature columns with original names (no prefixes)
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Groups cells by participant, timepoint, and cell type
#'   \item For RNA assay: uses VariableFeatures if available, otherwise all features
#'   \item Calculates counts: number of cells with expression > 0 for each feature
#'   \item Calculates mean expression: average expression across cells in each group
#'   \item Computes cell percentages within participant-timepoint combinations
#'   \item Compatible with both Seurat v4 and v5
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with single celltype column
#' result <- aggregate_seurat_for_ml_final(
#'   seurat_obj = pbmc_seurat,
#'   participant_col = "Participant ID",
#'   timepoint_col = "Timepoint (months)",
#'   celltype_col = "Prediction",
#'   assays = c("RNA", "ADT")
#' )
#' 
#' # Access results
#' counts_data <- result$counts
#' expression_data <- result$expression
#' 
#' # Example 2: Include additional metadata columns
#' result <- aggregate_seurat_for_ml_final(
#'   seurat_obj = pbmc_seurat,
#'   participant_col = "Participant ID",
#'   timepoint_col = "Timepoint (months)",
#'   celltype_col = "Prediction",
#'   columns_to_preserve = c("Unique Sample ID", "batch", "sex"),
#'   assays = c("RNA", "ADT")
#' )
#' 
#' # Example 3: Multiple celltype columns with Excel output
#' result <- aggregate_seurat_for_ml_final(
#'   seurat_obj = pbmc_seurat,
#'   participant_col = "Participant ID",
#'   timepoint_col = "Timepoint (months)",
#'   celltype_col = c("Prediction", "SingleR_fine"),
#'   columns_to_preserve = c("Unique Sample ID"),
#'   assays = c("RNA", "ADT"),
#'   output_xlsx = "aggregated_results.xlsx"
#' )
#' 
#' # Access results for each celltype annotation
#' prediction_counts <- result$counts_Prediction
#' prediction_expr <- result$expression_Prediction
#' singler_counts <- result$counts_SingleR_fine
#' singler_expr <- result$expression_SingleR_fine
#' 
#' # Example 4: Using SCT assay with variable features
#' result <- aggregate_seurat_for_ml_final(
#'   seurat_obj = sct_seurat,
#'   participant_col = "donor_id",
#'   timepoint_col = "timepoint",
#'   celltype_col = "cell_type",
#'   assays = c("SCT"),
#'   output_xlsx = "sct_aggregated.xlsx"
#' )
#' 
#' # Example 5: Validation - check specific group
#' # Validate counts for participant "40035", timepoint "12", celltype "CD4 Tn"
#' validate_gene <- "HBB"
#' manual_cells <- seurat_obj@meta.data[["Participant ID"]] == "40035" &
#'                 seurat_obj@meta.data[["Timepoint (months)"]] == "12" &
#'                 seurat_obj@meta.data[["Prediction"]] == "CD4 Tn"
#' 
#' rna_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
#' manual_count <- sum(rna_counts[validate_gene, manual_cells] > 0)
#' 
#' agg_count <- result$counts[
#'   result$counts[["Participant.ID"]] == "40035" &
#'   result$counts[["Timepoint..months."]] == "12" &
#'   result$counts[["Prediction"]] == "CD4 Tn",
#'   validate_gene
#' ]
#' 
#' print(paste("Manual count:", manual_count))
#' print(paste("Aggregated count:", agg_count))
#' }
#'
#' @importFrom Seurat GetAssayData VariableFeatures
#' @importFrom Matrix rowMeans rowSums
#' @importFrom data.table as.data.table data.table merge setnames
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'
#' @export
aggregate_seurat_for_ml_final <- function(
    seurat_obj,
    participant_col,
    timepoint_col,
    celltype_col,
    columns_to_preserve = NULL,
    assays = c("RNA", "ADT"),
    output_xlsx = NULL
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(openxlsx)
  })
  
  if (length(celltype_col) > 1) {
    message("Processing multiple celltype columns: ", paste(celltype_col, collapse = ", "))
    
    all_results <- list()
    
    for (ct_col in celltype_col) {
      message("\n=== Processing celltype column: ", ct_col, " ===")
      
      result_single <- aggregate_seurat_for_ml_final(
        seurat_obj = seurat_obj,
        participant_col = participant_col,
        timepoint_col = timepoint_col,
        celltype_col = ct_col,
        columns_to_preserve = columns_to_preserve,
        assays = assays,
        output_xlsx = NULL
      )
      
      ct_name <- make.names(ct_col)
      all_results[[paste0("counts_", ct_name)]] <- result_single$counts
      all_results[[paste0("expression_", ct_name)]] <- result_single$expression
    }
    
    if (!is.null(output_xlsx)) {
      wb <- createWorkbook()
      
      for (sheet_name in names(all_results)) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, all_results[[sheet_name]])
      }
      
      saveWorkbook(wb, output_xlsx, overwrite = TRUE)
      message("\nSaved all results to: ", output_xlsx)
      message("Sheets created: ", paste(names(all_results), collapse = ", "))
    }
    
    return(all_results)
  }
  
  meta <- seurat_obj@meta.data
  group_cols <- c(participant_col, timepoint_col, celltype_col)
  
  missing <- setdiff(group_cols, colnames(meta))
  if (length(missing) > 0) {
    stop("These grouping columns are missing from metadata: ",
         paste(missing, collapse = ", "))
  }
  
  if (!is.null(columns_to_preserve)) {
    columns_to_preserve <- setdiff(columns_to_preserve, group_cols)
    missing_pres <- setdiff(columns_to_preserve, colnames(meta))
    if (length(missing_pres) > 0) {
      stop("These columns_to_preserve are missing from metadata: ",
           paste(missing_pres, collapse = ", "))
    }
  } else {
    columns_to_preserve <- character(0)
  }
  
  meta_dt <- as.data.table(meta)
  meta_dt[, group_id := .GRP, by = group_cols]
  
  key_cols <- c("group_id", group_cols, columns_to_preserve)
  key_table <- meta_dt[, lapply(.SD, function(x) x[1]), by = group_id, .SDcols = setdiff(key_cols, "group_id")]
  
  cell_counts <- meta_dt[, .(total_cells = .N), by = group_id]
  
  if (participant_col %in% group_cols && timepoint_col %in% group_cols) {
    participant_totals <- meta_dt[, .(participant_total = .N), by = c(participant_col, timepoint_col)]
    cell_counts_with_context <- merge(
      cell_counts,
      key_table[, c("group_id", participant_col, timepoint_col), with = FALSE],
      by = "group_id"
    )
    cell_counts_with_context <- merge(
      cell_counts_with_context,
      participant_totals,
      by = c(participant_col, timepoint_col)
    )
    cell_counts_with_context[, cell_percentage := (total_cells / participant_total) * 100]
    cell_counts <- cell_counts_with_context[, .(group_id, total_cells, cell_percentage)]
  } else {

    total_cells_overall <- nrow(meta_dt)
    cell_counts[, cell_percentage := (total_cells / total_cells_overall) * 100]
  }
  
  count_results <- list()
  expression_results <- list()
  ncells_summary <- data.table(group_id = unique(meta_dt$group_id))
  
  for (assay_name in assays) {
    if (!assay_name %in% names(seurat_obj@assays)) {
      warning(paste("Skipping assay (not found):", assay_name))
      next
    }
    
    assay <- seurat_obj[[assay_name]]
    
    if (assay_name == "RNA" && length(VariableFeatures(seurat_obj)) > 0) {
      features_to_use <- VariableFeatures(seurat_obj)
      message("Using ", length(features_to_use), " variable features for RNA assay")
    } else {
      features_to_use <- rownames(assay)
      message("Using all ", length(features_to_use), " features for ", assay_name, " assay")
    }
    
    missing_feats   <- setdiff(features_to_use, rownames(assay))
    features_to_use <- intersect(features_to_use, rownames(assay))
    if (length(missing_feats) > 0)
      message("  Dropping ", length(missing_feats), " feature(s) not present in assay '",
              assay_name, "': ", paste(utils::head(missing_feats, 5), collapse = ", "),
              if (length(missing_feats) > 5) ", ..." else "")
    if (length(features_to_use) == 0)
      stop(sprintf("No usable features for assay '%s'.", assay_name))
    
    m_expr <- tryCatch({
      GetAssayData(assay, layer = "data")
    }, error = function(e) {
      GetAssayData(assay, slot = "data")
    })

    if (inherits(m_expr, "IterableMatrix")) {
      ridx   <- match(features_to_use, rownames(m_expr))
      ridx   <- ridx[!is.na(ridx)]
      m_expr <- as(m_expr[ridx, ], "dgCMatrix")
    } else {
      m_expr <- m_expr[features_to_use, , drop = FALSE]
    }
    
    m_counts <- tryCatch({
      GetAssayData(assay, layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(assay, slot = "counts")
      }, error = function(e2) {
        m_expr
      })
    })
    if (inherits(m_counts, "IterableMatrix")) {
      ridx     <- match(features_to_use, rownames(m_counts))
      ridx     <- ridx[!is.na(ridx)]
      m_counts <- as(m_counts[ridx, ], "dgCMatrix")
    } else {
      m_counts <- m_counts[features_to_use, , drop = FALSE]
    }
    
    gvec <- meta_dt$group_id
    
    dt_map <- data.table(
      cell = colnames(seurat_obj),
      group_id = gvec,
      idx = seq_len(ncol(seurat_obj))
    )
    idx_list <- split(dt_map$idx, dt_map$group_id)
    
    features_to_use <- intersect(rownames(m_counts), rownames(m_expr))
    if (length(features_to_use) == 0)
      stop(sprintf("No features left for assay '%s' after subsetting.", assay_name))
    m_counts <- m_counts[features_to_use, , drop = FALSE]
    m_expr   <- m_expr[features_to_use, , drop = FALSE]
    
    .aggregate_groups <- function(mat, f_one, f_many) {
      out <- matrix(0,
                    nrow = length(idx_list), ncol = nrow(mat),
                    dimnames = list(names(idx_list), rownames(mat)))
      for (i in seq_along(idx_list)) {
        idx <- idx_list[[i]]
        out[i, ] <- if (length(idx) == 1L) f_one(mat[, idx]) else f_many(mat[, idx, drop = FALSE])
      }
      out
    }
    
    counts_aggregated <- .aggregate_groups(
      m_counts,
      f_one  = function(v) as.numeric(v > 0),
      f_many = function(m) Matrix::rowSums(m > 0))
    
    df_counts <- as.data.table(counts_aggregated)
    df_counts[, group_id := as.integer(names(idx_list))]
    
    ncells_col <- paste0(tolower(assay_name), "_ncells")
    ncells_per_group <- df_counts[, .(ncells = rowSums(.SD, na.rm = TRUE)), .SDcols = features_to_use]
    ncells_per_group[, group_id := df_counts$group_id]
    setnames(ncells_per_group, "ncells", ncells_col)
    
    ncells_summary <- merge(ncells_summary, ncells_per_group, by = "group_id", all.x = TRUE)
    
    count_results[[assay_name]] <- df_counts
    
    expr_aggregated <- .aggregate_groups(
      m_expr,
      f_one  = function(v) as.numeric(v),
      f_many = function(m) Matrix::rowMeans(m))
    
    df_expr <- as.data.table(expr_aggregated)
    df_expr[, group_id := as.integer(names(idx_list))]
    
    expression_results[[assay_name]] <- df_expr
  }
  
  counts_dt <- key_table
  counts_dt <- merge(counts_dt, cell_counts, by = "group_id", all.x = TRUE)
  counts_dt <- merge(counts_dt, ncells_summary, by = "group_id", all.x = TRUE)
  
  for (assay_name in names(count_results)) {
    counts_dt <- merge(counts_dt, count_results[[assay_name]], by = "group_id", all.x = TRUE)
  }
  
  counts_dt[, group_id := NULL]
  counts_dt[, total_cells := NULL]
  
  expression_dt <- key_table
  expression_dt <- merge(expression_dt, cell_counts, by = "group_id", all.x = TRUE)
  expression_dt <- merge(expression_dt, ncells_summary, by = "group_id", all.x = TRUE)
  
  for (assay_name in names(expression_results)) {
    expression_dt <- merge(expression_dt, expression_results[[assay_name]], by = "group_id", all.x = TRUE)
  }
  
  expression_dt[, group_id := NULL]
  expression_dt[, total_cells := NULL]
  
  base_cols <- c(participant_col, timepoint_col, celltype_col, columns_to_preserve)
  ncells_cols <- grep("_ncells$", colnames(counts_dt), value = TRUE)
  meta_cols <- c(base_cols, "cell_percentage", ncells_cols)
  
  feature_cols_counts <- setdiff(colnames(counts_dt), meta_cols)
  counts_dt <- counts_dt[, c(meta_cols, feature_cols_counts), with = FALSE]
  
  feature_cols_expr <- setdiff(colnames(expression_dt), meta_cols)
  expression_dt <- expression_dt[, c(meta_cols, feature_cols_expr), with = FALSE]
  
  setnames(counts_dt, make.names(colnames(counts_dt), unique = TRUE))
  setnames(expression_dt, make.names(colnames(expression_dt), unique = TRUE))
  
  if (!is.null(output_xlsx)) {
    wb <- createWorkbook()
    
    addWorksheet(wb, "counts")
    writeData(wb, "counts", counts_dt)
    
    addWorksheet(wb, "expression")
    writeData(wb, "expression", expression_dt)
    
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    message("Saved output to: ", output_xlsx)
  }
  
  return(list(
    counts = counts_dt,
    expression = expression_dt
  ))
}
