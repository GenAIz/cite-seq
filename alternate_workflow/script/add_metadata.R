library(Seurat)
library(dplyr)
library(tibble)

#' Add External Metadata to Seurat Object
#'
#' @description
#' Integrates external metadata from various file formats (Excel, CSV, TSV) into a
#' Seurat object by matching on a specified key column. The function performs a left join
#' to preserve all cells in the Seurat object and adds matching metadata where available.
#' Compatible with both Seurat v4 and v5.
#'
#' @details
#' This function reads metadata from external files and joins it with the existing
#' Seurat object metadata based on a common identifier. The workflow includes:
#' 
#' **Supported File Formats**:
#' - Excel files (.xlsx, .xls): Requires readxl package
#' - CSV files (.csv): Uses read.csv
#' - TSV/Tab-delimited files (.tsv, .txt): Uses read.delim
#' 
#' **Join Behavior**:
#' - Performs left join to retain all cells in Seurat object
#' - Cells without matching metadata receive NA values
#' - Duplicate entries in metadata file trigger warnings
#' - Existing columns with same names are overwritten with warning
#' 
#' **Key Validation**:
#' - Checks that join keys exist in both datasets
#' - Validates join success and reports matching statistics
#' - Handles missing values in join keys
#' - Optionally filters metadata to only relevant entries
#' 
#' **Safety Features**:
#' - Preserves original row names (cell barcodes)
#' - Creates backup of original metadata (optional)
#' - Validates file existence before reading
#' - Comprehensive error handling and informative messages
#'
#' @param seurat_object Seurat object to add metadata to. Must contain existing metadata
#' with at least one column matching the object_key parameter.
#' @param metadata_file Character. Path to the metadata file. Supports:
#' - Excel: .xlsx, .xls
#' - CSV: .csv
#' - TSV/TXT: .tsv, .txt
#' @param object_key Character. Column name in the object metadata to use for joining
#' (e.g., "MULTI_ID", "orig.ident", "sample_id"). Must exist in seurat_object@meta.data.
#' @param metadata_key Character. Column name in the external metadata file to match with
#' object_key (e.g., "Specimen Barcode", "Sample_ID", "Barcode"). Must exist in metadata file.
#' @param sheet Character or numeric. For Excel files, specify the sheet name or number
#' (default: 1, first sheet). Ignored for non-Excel files.
#' @param col_types Character. Column type specification for reading the file:
#' - For Excel (readxl): "text", "numeric", "date", "logical", or "guess"
#' - For CSV/TSV: Standard R column types
#' - Default: "text" (reads all columns as character to preserve IDs)
#' @param subset_to_existing Logical. If TRUE (default), filters the metadata file to
#' only include rows matching values present in the Seurat object. Reduces memory
#' usage for large metadata files. If FALSE, reads entire metadata file.
#' @param overwrite_existing Logical. If TRUE, silently overwrites existing columns with
#' same names. If FALSE (default), warns before overwriting but proceeds with join.
#' @param add_match_status Logical. If TRUE (default), adds a column "metadata_matched"
#' indicating whether each cell successfully matched metadata (TRUE/FALSE). Useful for
#' QC and identifying cells without metadata.
#' @param backup_metadata Logical. If TRUE, stores original metadata in
#' seurat_object@misc$original_metadata before adding new data (default: FALSE).
#' Useful for reverting changes if needed.
#' @param na_values Character vector. Values to treat as NA when reading the file
#' (default: c("", "NA", "N/A", "na", "n/a")). File-format specific.
#' @param verbose Logical. If TRUE (default), prints detailed progress messages including
#' number of matches, columns added, and warnings. Set to FALSE for silent execution.
#'
#' @return Modified Seurat object with additional metadata columns. Original cell order
#' and row names are preserved. Returns:
#' \describe{
#' \item{@meta.data}{Updated with joined metadata columns}
#' \item{@misc$original_metadata}{Original metadata (if backup_metadata = TRUE)}
#' \item{@misc$metadata_join_stats}{Join statistics including match rate}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with Excel file
#' # Add metadata from Excel sheet matching on MULTI_ID
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "sample_metadata.xlsx",
#' object_key = "MULTI_ID",
#' metadata_key = "Specimen Barcode",
#' sheet = "unique ID"
#' )
#' 
#' # Example 2: CSV file with custom column types
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "patient_data.csv",
#' object_key = "sample_id",
#' metadata_key = "Sample_ID",
#' col_types = "text"
#' )
#' 
#' # Example 3: TSV file with all features enabled
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "clinical_data.tsv",
#' object_key = "patient_id",
#' metadata_key = "Patient_ID",
#' subset_to_existing = TRUE,
#' add_match_status = TRUE,
#' backup_metadata = TRUE,
#' verbose = TRUE
#' )
#' 
#' # Example 4: CSV file without subsetting
#' # Load entire metadata file (useful for small files)
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "metadata.csv",
#' object_key = "orig.ident",
#' metadata_key = "sample",
#' subset_to_existing = FALSE
#' )
#' 
#' # Example 5: Excel file with numeric sheet index
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "multi_sheet_metadata.xlsx",
#' object_key = "MULTI_ID",
#' metadata_key = "Barcode",
#' sheet = 2, # Second sheet
#' col_types = "guess"
#' )
#' 
#' # Example 6: Silent execution for batch processing
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "batch_metadata.csv",
#' object_key = "sample",
#' metadata_key = "Sample_Name",
#' verbose = FALSE
#' )
#' 
#' # Example 7: With custom NA values
#' seurat_obj <- add_metadata(
#' seurat_object = seurat_obj,
#' metadata_file = "data_with_missing.xlsx",
#' object_key = "MULTI_ID",
#' metadata_key = "ID",
#' sheet = "Sheet1",
#' na_values = c("", "NA", "missing", "unknown", ".")
#' )
#' 
#' # Check results:
#' # View new columns
#' head(seurat_obj@meta.data)
#' 
#' # Check matching statistics
#' table(seurat_obj$metadata_matched)
#' 
#' # View join statistics
#' seurat_obj@misc$metadata_join_stats
#' 
#' # Restore original metadata if needed (when backup_metadata = TRUE)
#' seurat_obj@meta.data <- seurat_obj@misc$original_metadata
#' }
#'
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom readxl read_excel
#' @importFrom utils read.csv read.delim
#' 
#' @export
add_metadata <- function(
    seurat_object,
    metadata_file,
    object_key,
    metadata_key,
    sheet = 1,
    col_types = "text",
    subset_to_existing = TRUE,
    overwrite_existing = FALSE,
    add_match_status = TRUE,
    backup_metadata = FALSE,
    na_values = c("", "NA", "N/A", "na", "n/a"),
    verbose = TRUE
) {
  
  if (verbose) cat("=== Adding Metadata to Seurat Object ===\n\n")
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  if (!object_key %in% colnames(seurat_object@meta.data)) {
    stop(sprintf("Key '%s' not found in the object metadata. Available columns: %s",
                 object_key,
                 paste(colnames(seurat_object@meta.data), collapse = ", ")))
  }
  
  if (backup_metadata) {
    if (verbose) cat("Creating backup of original metadata...\n")
    seurat_object@misc$original_metadata <- seurat_object@meta.data
  }
  
  if (verbose) cat(sprintf("Reading metadata from: %s\n", basename(metadata_file)))
  
  file_ext <- tolower(tools::file_ext(metadata_file))
  
  metadata <- tryCatch({
    if (file_ext %in% c("xlsx", "xls")) {

      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required for Excel files. Install with: install.packages('readxl')")
      }
      if (verbose) cat(sprintf("Reading Excel sheet: %s\n", sheet))
      readxl::read_excel(
        path = metadata_file,
        sheet = sheet,
        col_types = col_types,
        na = na_values
      )
    } else if (file_ext == "csv") {

      read.csv(
        file = metadata_file,
        stringsAsFactors = FALSE,
        na.strings = na_values
      )
    } else if (file_ext %in% c("tsv", "txt")) {

      read.delim(
        file = metadata_file,
        stringsAsFactors = FALSE,
        na.strings = na_values
      )
    } else {
      stop(sprintf("Unsupported file format: .%s. Supported formats: xlsx, xls, csv, tsv, txt", 
                   file_ext))
    }
  }, error = function(e) {
    stop("Error reading metadata file: ", e$message)
  })
  
  metadata <- as.data.frame(metadata)
  
  if (verbose) {
    cat(sprintf("Metadata loaded: %d rows, %d columns\n", 
                nrow(metadata), ncol(metadata)))
  }
  
  if (!metadata_key %in% colnames(metadata)) {
    stop(sprintf("Key '%s' not found in metadata file. Available columns: %s",
                 metadata_key,
                 paste(colnames(metadata), collapse = ", ")))
  }
  
  duplicates <- duplicated(metadata[[metadata_key]])
  if (any(duplicates)) {
    n_duplicates <- sum(duplicates)
    warning(sprintf("Found %d duplicate entries in metadata key '%s'. Only the first occurrence of each is used.",
                    n_duplicates, metadata_key))
    if (verbose) {
      cat(sprintf("WARNING: %d duplicate entries found in metadata - dropping them\n", n_duplicates))
    }
    metadata <- metadata[!duplicates, , drop = FALSE]
  }
  
  if (subset_to_existing) {
    if (verbose) cat("Subsetting metadata to existing values in Seurat object...\n")
    
    unique_ids <- unique(seurat_object@meta.data[[object_key]])
    unique_ids <- unique_ids[!is.na(unique_ids)]
    
    n_before <- nrow(metadata)
    metadata <- metadata[metadata[[metadata_key]] %in% unique_ids, ]
    n_after <- nrow(metadata)
    
    if (verbose) {
      cat(sprintf("Metadata subset: %d -> %d rows (%.1f%% retained)\n",
                  n_before, n_after, 100 * n_after / n_before))
    }
    
    if (n_after == 0) {
      stop("No matching entries found between Seurat object and metadata file. Check key columns.")
    }
  }
  
  existing_cols <- colnames(seurat_object@meta.data)
  new_cols <- setdiff(colnames(metadata), metadata_key)
  overlapping_cols <- intersect(existing_cols, new_cols)
  
  if (length(overlapping_cols) > 0) {
    if (overwrite_existing) {
      warning(sprintf("Overwriting %d existing metadata column(s) with values from the file: %s",
                      length(overlapping_cols), paste(overlapping_cols, collapse = ", ")))
      if (verbose)
        cat(sprintf("Overwriting %d existing column(s): %s\n",
                    length(overlapping_cols), paste(overlapping_cols, collapse = ", ")))
      seurat_object@meta.data <- seurat_object@meta.data[
        , setdiff(colnames(seurat_object@meta.data), overlapping_cols), drop = FALSE]
    } else {
      warning(sprintf(paste0("%d column(s) already exist on the object and were NOT imported ",
                             "(set overwrite_existing = TRUE to replace them): %s"),
                      length(overlapping_cols), paste(overlapping_cols, collapse = ", ")))
      if (verbose)
        cat(sprintf("Keeping existing values; skipping %d incoming column(s): %s\n",
                    length(overlapping_cols), paste(overlapping_cols, collapse = ", ")))
      metadata <- metadata[, setdiff(colnames(metadata), overlapping_cols), drop = FALSE]
    }
  }
  
  if (verbose) cat(sprintf("\nJoining metadata on '%s' = '%s'...\n", 
                           object_key, metadata_key))
  
  cell_meta <- seurat_object@meta.data %>%
    tibble::rownames_to_column("cell")
  
  n_cols_before <- ncol(cell_meta) - 1
  
  if (add_match_status) {

    metadata$temp_match_indicator <- TRUE
  }
  
  n_cells_before_join <- nrow(cell_meta)
  cell_meta <- cell_meta %>%
    dplyr::left_join(metadata, by = setNames(metadata_key, object_key))
  
  if (nrow(cell_meta) != n_cells_before_join) {
    stop(sprintf(paste0("Metadata join changed the cell count (%d -> %d). The key '%s' is not ",
                        "unique in the metadata file, so cells were duplicated. Deduplicate the ",
                        "metadata before joining."),
                 n_cells_before_join, nrow(cell_meta), metadata_key))
  }
  
  n_cols_after <- ncol(cell_meta) - 1
  n_cols_added <- n_cols_after - n_cols_before
  
  if (add_match_status) {

    cell_meta$metadata_matched <- !is.na(cell_meta$temp_match_indicator)

    cell_meta$temp_match_indicator <- NULL
    
    if (verbose) {
      n_matched <- sum(cell_meta$metadata_matched)
      n_total <- nrow(cell_meta)
      cat(sprintf("Match status: %d/%d cells matched (%.1f%%)\n",
                  n_matched, n_total, 100 * n_matched / n_total))
    }
  }
  
  cell_meta <- cell_meta %>%
    tibble::column_to_rownames("cell")
  
  seurat_object@meta.data <- cell_meta
  
  join_stats <- list(
    metadata_file = metadata_file,
    object_key = object_key,
    metadata_key = metadata_key,
    n_cells_total = nrow(cell_meta),
    n_cells_matched = if(add_match_status) sum(cell_meta$metadata_matched) else NA,
    match_rate = if(add_match_status) mean(cell_meta$metadata_matched) else NA,
    n_columns_added = n_cols_added,
    columns_added = setdiff(colnames(cell_meta), c(existing_cols, "metadata_matched")),
    timestamp = Sys.time()
  )
  
  seurat_object@misc$metadata_join_stats <- join_stats
  
  if (verbose) {
    cat("\n=== Metadata Integration Complete ===\n")
    cat(sprintf("Columns added: %d\n", n_cols_added))
    cat(sprintf("New columns: %s\n", 
                paste(join_stats$columns_added, collapse = ", ")))
    if (add_match_status) {
      cat(sprintf("Cells matched: %d/%d (%.1f%%)\n",
                  join_stats$n_cells_matched,
                  join_stats$n_cells_total,
                  100 * join_stats$match_rate))
    }
    cat(sprintf("Total metadata columns: %d\n", ncol(seurat_object@meta.data)))
    if (backup_metadata) {
      cat("Original metadata backed up to @misc$original_metadata\n")
    }
    cat("Join statistics stored in @misc$metadata_join_stats\n")
  }
  
  return(seurat_object)
}

#' View Metadata Join Statistics
#'
#' @description
#' Displays statistics from the most recent metadata join operation performed
#' using add_metadata().
#'
#' @param seurat_object Seurat object with metadata join statistics
#'
#' @return Invisibly returns the join statistics list
#' @export
#'
#' @examples
#' \dontrun{
#' view_metadata_join_stats(seurat_obj)
#' }
view_metadata_join_stats <- function(seurat_object) {
  if (is.null(seurat_object@misc$metadata_join_stats)) {
    message("No metadata join statistics found. Run add_metadata() first.")
    return(invisible(NULL))
  }
  
  stats <- seurat_object@misc$metadata_join_stats
  
  cat("=== Metadata Join Statistics ===\n")
  cat(sprintf("File: %s\n", basename(stats$metadata_file)))
  cat(sprintf("Join keys: %s = %s\n", stats$object_key, stats$metadata_key))
  cat(sprintf("Timestamp: %s\n", stats$timestamp))
  cat(sprintf("Total cells: %s\n", format(stats$n_cells_total, big.mark = ",")))
  if (!is.na(stats$n_cells_matched)) {
    cat(sprintf("Matched cells: %s (%.1f%%)\n", 
                format(stats$n_cells_matched, big.mark = ","),
                100 * stats$match_rate))
  }
  cat(sprintf("Columns added: %d\n", stats$n_columns_added))
  cat(sprintf("New columns: %s\n", paste(stats$columns_added, collapse = ", ")))
  
  invisible(stats)
}

#' Restore Original Metadata
#'
#' @description
#' Restores the original metadata from backup created when backup_metadata = TRUE
#' was used in add_metadata().
#'
#' @param seurat_object Seurat object with backed up metadata
#'
#' @return Seurat object with restored original metadata
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- restore_original_metadata(seurat_obj)
#' }
restore_original_metadata <- function(seurat_object) {
  if (is.null(seurat_object@misc$original_metadata)) {
    stop("No backup metadata found. Set backup_metadata = TRUE when adding metadata.")
  }
  
  cat("Restoring original metadata...\n")
  n_cols_before <- ncol(seurat_object@meta.data)
  
  seurat_object@meta.data <- seurat_object@misc$original_metadata
  
  n_cols_after <- ncol(seurat_object@meta.data)
  cat(sprintf("Metadata restored: %d -> %d columns\n", n_cols_before, n_cols_after))
  
  return(seurat_object)
}
