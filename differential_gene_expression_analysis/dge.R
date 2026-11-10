#' Candidate spellings of a group level after Seurat aggregation
#'
#' `Seurat::AggregateExpression(group.by = c(a, b))` joins group levels with
#' an underscore, so a column is named `<a>_<b>`. Recovering `(a, b)` by
#' splitting on `_` is impossible whenever a level itself contains an
#' underscore ("Severe_DISEASE", "Non_classical", "Day_0"): earlier code split
#' and took the first token, silently assigning every sample to an invented
#' group. Columns are instead matched against candidate names generated from
#' the levels that actually exist in the metadata; this helper enumerates the
#' spellings Seurat and `make.names()` can produce for one level.
#'
#' @param v A group level (any scalar; coerced to character).
#' @return Character vector of unique candidate spellings of `v`.
.aggregated_level_variants <- function(v) {
  v <- as.character(v)
  unique(c(
    v,
    make.names(v),
    gsub("_", "-", v),
    gsub("[^A-Za-z0-9]+", "-", v),
    gsub("[^A-Za-z0-9]+", ".", v),
    paste0("g", v),
    paste0("X", v),
    paste0("X", gsub("[^A-Za-z0-9._]", ".", v))
  ))
}

#' Map aggregated pseudobulk columns back to their group levels
#'
#' Matches each aggregated column name against whole candidate strings built
#' from the group levels present in `meta` (via
#' [.aggregated_level_variants()]), so embedded underscores in levels are
#' irrelevant. Anything unmatched or ambiguous is a hard error: guessing here
#' is what caused the original silently-wrong-biology defect.
#'
#' @param col_names Column names produced by `Seurat::AggregateExpression`.
#' @param meta Cell metadata holding the grouping columns.
#' @param group_cols The grouping columns used for aggregation, in order.
#' @param verbose Print a one-line summary of the mapping.
#' @return A data frame with one row per column of `col_names` and one
#'   column per entry of `group_cols`, holding the recovered levels.
.map_aggregated_columns <- function(col_names, meta, group_cols, verbose = FALSE) {
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols) > 0)
    stop("Grouping column(s) not found in metadata: ", paste(missing_cols, collapse = ", "))

  combos <- unique(as.data.frame(lapply(meta[, group_cols, drop = FALSE], as.character),
                                 stringsAsFactors = FALSE))
  colnames(combos) <- group_cols

  cand <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(combos))) {
    per  <- lapply(group_cols, function(g) .aggregated_level_variants(combos[i, g]))
    grid <- expand.grid(per, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    nms  <- unique(do.call(paste, c(as.list(grid), sep = "_")))
    for (n in nms) assign(n, c(get0(n, envir = cand, ifnotfound = integer(0)), i), envir = cand)
  }

  idx <- vapply(col_names, function(cn) {
    hit <- get0(cn, envir = cand, ifnotfound = integer(0))
    hit <- unique(hit)
    if (length(hit) == 0) return(NA_integer_)
    if (length(hit) > 1)  return(-1L)
    hit
  }, integer(1))

  bad <- col_names[is.na(idx)]
  amb <- col_names[!is.na(idx) & idx < 0]
  if (length(bad) > 0 || length(amb) > 0) {
    stop(sprintf(paste0(
      "Could not map %d aggregated column(s) back to their %s levels.\n",
      "  unmatched: %s\n  ambiguous: %s\n",
      "This usually means a group level contains a character Seurat rewrites in a ",
      "way not covered here. Fix the level names rather than letting samples be ",
      "assigned to the wrong group."),
      length(bad) + length(amb), paste(group_cols, collapse = " x "),
      if (length(bad)) paste(utils::head(bad, 8), collapse = ", ") else "none",
      if (length(amb)) paste(utils::head(amb, 8), collapse = ", ") else "none"))
  }

  out <- combos[idx, , drop = FALSE]
  rownames(out) <- col_names
  if (verbose)
    message(sprintf("  Mapped %d aggregated column(s) to %s levels",
                    length(col_names), paste(group_cols, collapse = " x ")))
  out
}

#' Perform Differential Gene Expression Analysis
#'
#' Supports Seurat-native tests (FindAllMarkers/FindMarkers) and bulk methods
#' (DESeq2/limma/edgeR) via pseudobulking with proper design matrix and contrast
#' specifications.
#'
#' For pseudobulk methods, two design strategies are available:
#' \itemize{
#'   \item \strong{stratified}: Subset data by timepoint and fit separate models
#'     with condition as the factor. Simpler approach, analyzes each
#'     timepoint independently.
#'   \item \strong{interaction}: Combine timepoint and condition into a
#'     single interaction term (Timepoint_Condition). More powerful,
#'     accounts for repeated measures via \code{duplicateCorrelation()}, all
#'     contrasts in one model.
#' }
#'
#' @param seurat_obj Seurat object or named list of Seurat objects.
#' @param test.use Character. One of "wilcox", "bimod", "roc", "t", "negbinom",
#'   "poisson", "LR", "MAST", "DESeq2", "limma", "edger", "edgeR".
#'   Default: NULL (uses MAST).
#' @param pseudobulking Logical. Required TRUE for DESeq2/limma/edger.
#' @param pseudobulk_mode One of "cluster_sample", "sample_only",
#'   "celltype_condition".
#' @param design_strategy Character. One of "stratified" or "interaction".
#'   Only used for pseudobulk methods. Default: "stratified".
#'   \itemize{
#'     \item "stratified": Fit separate models per timepoint using condition.
#'     \item "interaction": Single model with Timepoint x Condition interaction,
#'       accounts for repeated measures.
#'   }
#' @param cluster_col Column for cluster/cell type (default "SingleR_main").
#' @param sample_col Column for sample/replicate (default "Participant ID").
#' @param condition_col Column for condition (default "Terminal Group").
#' @param timepoint_col Column for timepoint (default "Timepoint (months)"). Used for
#'   stratified and interaction designs.
#' @param participant_col Column for participant ID for repeated measures blocking
#'   (default same as sample_col). Used in interaction design for
#'   \code{duplicateCorrelation()}.
#' @param pseudobulk_assay Assay to aggregate counts from (default "RNA").
#' @param contrasts Named list of contrast specifications. If NULL, defaults are
#'   generated based on available groups. For custom contrasts, provide a list
#'   where names are contrast labels and values are character formulas
#'   (e.g., \code{list("GroupA_vs_GroupB" = "GroupA - GroupB")}).
#' @param timepoints_to_analyze Character or numeric vector of timepoints to
#'   analyze. Default NULL analyzes all available timepoints.
#' @param conditions_to_analyze Character vector of conditions to include.
#'   Default NULL includes all available conditions.
#' @param reference_condition Character. Reference condition for contrasts
#'   (denominator). If NULL, all pairwise comparisons are made.
#' @param ident.1,ident.2 For simple two-group comparisons (optional). When
#'   both are set, only this single contrast is computed; for cell-level
#'   methods, ident.1 alone compares that group against all other cells.
#' @param latent.vars Vector of latent covariates for Seurat-native methods
#'   (e.g., \code{c("nCount_RNA", "percent.mt")}). Only used for MAST/LR.
#' @param block_repeated_measures Logical. Whether to use
#'   \code{duplicateCorrelation()} for repeated measures in interaction design.
#'   Default TRUE.
#' @param min_samples_per_group Minimum samples per group for a contrast to be
#'   computed (default 2).
#' @param min_cells_per_pseudobulk Minimum cells required to form a pseudobulk
#'   sample (default 10).
#' @param lfc_threshold Log fold-change threshold for testing (default 0).
#' @param p_adjust_method Method for p-value adjustment (default "BH").
#' @param return_fit Logical. If TRUE, return fitted model objects in addition
#'   to results tables (default FALSE). Note: when TRUE, the function returns
#'   a list instead of a data frame.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param ... Additional arguments passed to FindMarkers/FindAllMarkers.
#'
#' @return A data frame containing all DE results with the following columns:
#' \itemize{
#'   \item \code{gene}: Gene name/symbol
#'   \item \code{cluster}: Cell cluster/type
#'   \item \code{timepoint}: Timepoint (for stratified design)
#'   \item \code{contrast}: Contrast name (e.g., "24mR_vs_24mO")
#'   \item Method-specific statistics (limma: logFC, AveExpr, t, P.Value, adj.P.Val, B;
#'     edgeR: logFC, logCPM, F, PValue, FDR; DESeq2: baseMean, log2FoldChange, pvalue, padj)
#' }
#' For multiple input objects, an additional \code{object} column is included.
#'
#' @section Input Validation:
#' The function performs comprehensive validation of input data:
#' \itemize{
#'   \item Checks that all specified metadata columns exist
#'   \item Reports unique values in condition, timepoint, and cluster columns
#'   \item Warns about NA values in key columns
#'   \item Validates minimum sample sizes per group
#'   \item Reports data structure summary before analysis
#' }
#'
#' @section Design Matrix Options:
#' \subsection{Option 1 - Stratified by Timepoint}{
#' Uses condition directly as the factor, then subsets data by Timepoint:
#' \preformatted{
#' design <- model.matrix(~ 0 + condition, data = metadata_subset)
#' }
#' }
#' \subsection{Option 2 - Interaction Model}{
#' Combines Timepoint and Condition into a single factor:
#' \preformatted{
#' metadata$Group <- interaction(Timepoint, Condition, sep = "_")
#' design <- model.matrix(~ 0 + Group, data = metadata)
#' }
#' }
#'
#' @section Contrast Specification:
#' Contrasts can be specified as a named list of formulas:
#' \preformatted{
#' contrasts <- list(
#'   "GroupA_vs_GroupB_at_T1" = "T1_GroupA - T1_GroupB",
#'   "GroupA_vs_GroupB_at_T2" = "T2_GroupA - T2_GroupB"
#' )
#' }
#'
#' @importFrom Seurat AggregateExpression FindAllMarkers FindMarkers
#' @importFrom limma makeContrasts voom lmFit contrasts.fit eBayes topTable
#'   duplicateCorrelation
#' @importFrom edgeR DGEList calcNormFactors filterByExpr estimateDisp glmQLFit
#'   glmQLFTest topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink counts
#' @importFrom dplyr bind_rows mutate
#' @importFrom stats model.matrix as.formula setNames
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Stratified design - compare conditions at each timepoint
#' results <- perform_dge_analysis(
#'   seurat_obj = my_seurat,
#'   test.use = "limma",
#'   pseudobulking = TRUE,
#'   design_strategy = "stratified",
#'   cluster_col = "SingleR_main",
#'   sample_col = "Participant.ID",
#'   condition_col = "Terminal.Group",
#'   timepoint_col = "Timepoint.(months)",
#'   timepoints_to_analyze = c(12, 24),
#'   conditions_to_analyze = c("24mR", "24mO"),
#'   ident.1 = "24mR",
#'   ident.2 = "24mO",
#'   verbose = TRUE
#' )
#' 
#' # Example 2: Compare timepoints within a condition group
#' # First subset to the condition of interest
#' obj_24mR <- subset(my_seurat, subset = Terminal.Group == "24mR")
#' 
#' results_timepoint <- perform_dge_analysis(
#'   seurat_obj = obj_24mR,
#'   test.use = "limma",
#'   pseudobulking = TRUE,
#'   design_strategy = "stratified",
#'   cluster_col = "SingleR_main",
#'   sample_col = "Participant.ID",
#'   condition_col = "Timepoint.(months)",   # Treat timepoint as condition
#'   timepoint_col = "Timepoint.(months)",
#'   timepoints_to_analyze = NULL,           # Don't stratify - triggers direct mode
#'   conditions_to_analyze = c("3", "12", "24"),
#'   reference_condition = "3",              # Baseline timepoint
#'   verbose = TRUE
#' )
#'
#' # Example 3: Interaction design with repeated measures correction
#' results_interaction <- perform_dge_analysis(
#'   seurat_obj = my_seurat,
#'   test.use = "limma",
#'   pseudobulking = TRUE,
#'   design_strategy = "interaction",
#'   cluster_col = "SingleR_main",
#'   sample_col = "Participant.ID",
#'   condition_col = "Terminal.Group",
#'   timepoint_col = "Timepoint.(months)",
#'   participant_col = "Participant.ID",
#'   block_repeated_measures = TRUE,
#'   verbose = TRUE
#' )
#' 
#' # Export results to CSV
#' write.csv(results, "DEG_results.csv", row.names = FALSE)
#' }
perform_dge_analysis <- function(seurat_obj,
                                 test.use = NULL,
                                 pseudobulking = FALSE,
                                 pseudobulk_mode = c("cluster_sample", "sample_only", "celltype_condition"),
                                 design_strategy = c("stratified", "interaction"),
                                 cluster_col = "SingleR_main",
                                 sample_col = "Participant ID",
                                 condition_col = "Terminal Group",
                                 timepoint_col = "Timepoint (months)",
                                 participant_col = NULL,
                                 pseudobulk_assay = "RNA",
                                 contrasts = NULL,
                                 timepoints_to_analyze = NULL,
                                 conditions_to_analyze = NULL,
                                 reference_condition = NULL,
                                 ident.1 = NULL,
                                 ident.2 = NULL,
                                 latent.vars = NULL,
                                 block_repeated_measures = TRUE,
                                 min_samples_per_group = 2,
                                 min_cells_per_pseudobulk = 10,
                                 lfc_threshold = 0,
                                 p_adjust_method = "BH",
                                 return_fit = FALSE,
                                 verbose = TRUE,
                                 ...) {
  
  pseudobulk_mode <- match.arg(pseudobulk_mode)
  design_strategy <- match.arg(design_strategy)
  
  if (is.null(participant_col)) {
    participant_col <- sample_col
  }
  
  if (is.null(test.use)) test.use <- "MAST"
  
  test.use_input <- test.use
  test.use <- tolower(test.use)
  valid_tests <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson",
                   "lr", "mast", "deseq2", "limma", "edger")
  if (!test.use %in% valid_tests) {
    stop("Invalid test.use: '", test.use_input, "'. Must be one of: ",
         paste(valid_tests, collapse = ", "))
  }
  
  bulk_methods <- c("deseq2", "limma", "edger")
  is_bulk_method <- test.use %in% bulk_methods
  
  if (is_bulk_method && !pseudobulking) {
    stop("Bulk methods (DESeq2/limma/edgeR) require pseudobulking = TRUE.")
  }
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  multi_input <- is.list(seurat_obj) && !inherits(seurat_obj, "Seurat")
  objs <- if (multi_input) seurat_obj else list(seurat_obj)
  obj_names <- names(objs) %||% paste0("Object_", seq_along(objs))
  
  #' Validate that required columns exist in metadata
  validate_columns <- function(obj, required_cols, obj_name = "Object") {
    md_cols <- colnames(obj@meta.data)
    missing <- setdiff(required_cols, md_cols)
    
    if (length(missing) > 0) {
      stop(
        "Missing required metadata columns in ", obj_name, ":\n",
        "  Missing: ", paste(missing, collapse = ", "), "\n",
        "  Available columns: ", paste(md_cols, collapse = ", ")
      )
    }
    invisible(TRUE)
  }
  
  #' Summarize metadata structure
  summarize_metadata <- function(obj, cols, obj_name = "Object") {
    if (!verbose) return(invisible(NULL))
    
    message("\n", strrep("=", 60))
    message("DATA SUMMARY: ", obj_name)
    message(strrep("=", 60))
    message("Total cells: ", format(ncol(obj), big.mark = ","))
    
    for (col in cols) {
      if (!col %in% colnames(obj@meta.data)) next
      
      values <- obj@meta.data[[col]]
      n_na <- sum(is.na(values))
      unique_vals <- unique(values[!is.na(values)])
      
      message("\nColumn: '", col, "'")
      message("  Unique values (n=", length(unique_vals), "): ",
              paste(head(sort(as.character(unique_vals)), 10), collapse = ", "),
              if (length(unique_vals) > 10) " ..." else "")
      
      if (n_na > 0) {
        message("  WARNING: ", n_na, " NA values (",
                round(100 * n_na / length(values), 1), "%)")
      }
      
      if (length(unique_vals) <= 20) {
        tbl <- sort(table(values), decreasing = TRUE)
        message("  Distribution:")
        for (i in seq_along(tbl)) {
          message("    ", names(tbl)[i], ": ", format(tbl[i], big.mark = ","),
                  " (", round(100 * tbl[i] / sum(tbl), 1), "%)")
        }
      }
    }
    message(strrep("=", 60), "\n")
    invisible(NULL)
  }
  
  #' Validate condition values
  validate_conditions <- function(obj, condition_col, conditions_to_analyze = NULL,
                                  ident.1 = NULL, ident.2 = NULL) {
    available_conditions <- unique(obj@meta.data[[condition_col]])
    available_conditions <- available_conditions[!is.na(available_conditions)]
    
    if (!is.null(ident.1) && !ident.1 %in% available_conditions) {
      stop("ident.1 = '", ident.1, "' not found in column '", condition_col, "'.\n",
           "  Available values: ", paste(available_conditions, collapse = ", "))
    }
    if (!is.null(ident.2) && !ident.2 %in% available_conditions) {
      stop("ident.2 = '", ident.2, "' not found in column '", condition_col, "'.\n",
           "  Available values: ", paste(available_conditions, collapse = ", "))
    }
    
    if (!is.null(conditions_to_analyze)) {
      invalid <- setdiff(as.character(conditions_to_analyze), as.character(available_conditions))
      if (length(invalid) > 0) {
        stop("Invalid conditions_to_analyze: ", paste(invalid, collapse = ", "), "\n",
             "  Available values: ", paste(available_conditions, collapse = ", "))
      }
    }
    
    if (!is.null(reference_condition) && !reference_condition %in% as.character(available_conditions)) {
      stop("reference_condition = '", reference_condition,
           "' not found in column '", condition_col, "'.\n",
           "  Available values: ", paste(available_conditions, collapse = ", "))
    }
    
    available_conditions
  }
  
  #' Validate timepoint values
  validate_timepoints <- function(obj, timepoint_col, timepoints_to_analyze = NULL) {
    available_timepoints <- unique(obj@meta.data[[timepoint_col]])
    available_timepoints <- available_timepoints[!is.na(available_timepoints)]
    
    if (!is.null(timepoints_to_analyze)) {

      timepoints_chr <- as.character(timepoints_to_analyze)
      available_chr <- as.character(available_timepoints)
      
      invalid <- timepoints_chr[!timepoints_chr %in% available_chr]
      if (length(invalid) > 0) {
        stop("Invalid timepoints_to_analyze: ", paste(invalid, collapse = ", "), "\n",
             "  Available values: ", paste(available_timepoints, collapse = ", "))
      }
    }
    
    available_timepoints
  }
  
  #' Validate cluster values
  validate_clusters <- function(obj, cluster_col) {
    clusters <- unique(obj@meta.data[[cluster_col]])
    clusters <- clusters[!is.na(clusters)]
    
    if (length(clusters) == 0) {
      stop("No valid cluster values found in column '", cluster_col, "'")
    }
    
    clusters
  }
  
  #' Report pseudobulk sample structure
  report_pseudobulk_structure <- function(obj, cluster_col, sample_col,
                                          condition_col, timepoint_col = NULL) {
    if (!verbose) return(invisible(NULL))
    
    message("\n", strrep("-", 40))
    message("PSEUDOBULK SAMPLE STRUCTURE")
    message(strrep("-", 40))
    
    md <- obj@meta.data
    
    sample_condition <- unique(md[, c(sample_col, condition_col), drop = FALSE])
    sample_condition <- sample_condition[!is.na(sample_condition[[condition_col]]), ]
    cond_tbl <- table(sample_condition[[condition_col]])
    message("\nSamples per condition:")
    for (i in seq_along(cond_tbl)) {
      message("  ", names(cond_tbl)[i], ": n=", cond_tbl[i])
    }
    
    if (!is.null(timepoint_col) && timepoint_col %in% colnames(md)) {
      sample_tp_cond <- unique(md[, c(sample_col, timepoint_col, condition_col), drop = FALSE])
      sample_tp_cond <- sample_tp_cond[!is.na(sample_tp_cond[[condition_col]]) &
                                         !is.na(sample_tp_cond[[timepoint_col]]), ]
      
      message("\nSamples per Timepoint x Condition:")
      tp_cond_tbl <- table(sample_tp_cond[[timepoint_col]], sample_tp_cond[[condition_col]])
      print(tp_cond_tbl)
      
      low_n <- which(tp_cond_tbl > 0 & tp_cond_tbl < min_samples_per_group, arr.ind = TRUE)
      if (nrow(low_n) > 0) {
        message("\nWARNING: The following groups have < ", min_samples_per_group, " samples:")
        for (i in seq_len(nrow(low_n))) {
          tp <- rownames(tp_cond_tbl)[low_n[i, 1]]
          cond <- colnames(tp_cond_tbl)[low_n[i, 2]]
          n <- tp_cond_tbl[low_n[i, 1], low_n[i, 2]]
          message("  Timepoint ", tp, ", ", cond, ": n=", n)
        }
      }
    }
    
    cluster_tbl <- table(md[[cluster_col]])
    message("\nCells per cluster:")
    for (cl in names(sort(cluster_tbl, decreasing = TRUE))) {
      message("  ", cl, ": ", format(cluster_tbl[cl], big.mark = ","))
    }
    
    message(strrep("-", 40), "\n")
    invisible(NULL)
  }
  
  #' Sum counts into pseudobulk columns
  #'
  #' @param obj Seurat object (already subset to the cells in scope).
  #' @param group_by_cols Metadata columns whose combinations define the
  #'   pseudobulk units (passed to `Seurat::AggregateExpression`).
  #' @return Aggregated counts matrix (genes x pseudobulk units), or `NULL`
  #'   when nothing could be aggregated.
  aggregate_pseudobulk <- function(obj, group_by_cols) {
    if (verbose) {
      message("  Aggregating by: ", paste(group_by_cols, collapse = " x "))
    }
    
    missing <- setdiff(group_by_cols, colnames(obj@meta.data))
    if (length(missing) > 0) {
      stop("Cannot aggregate: missing columns ", paste(missing, collapse = ", "))
    }
    
    ag <- Seurat::AggregateExpression(
      object = obj,
      assays = pseudobulk_assay,
      group.by = group_by_cols,
      slot = "counts",
      return.seurat = FALSE,
      verbose = FALSE
    )[[pseudobulk_assay]]
    
    if (is.null(ag) || ncol(ag) == 0) {
      warning("Aggregation produced no samples.")
      return(NULL)
    }
    
    as.matrix(ag)
  }
  
  #' Map AggregateExpression column names back to their group levels via .map_aggregated_columns().
  build_pseudobulk_meta <- function(obj, counts_mat, group_by_cols) {
    cn <- colnames(counts_mat)

    meta_df <- .map_aggregated_columns(
      col_names  = cn,
      meta       = obj@meta.data,
      group_cols = group_by_cols,
      verbose    = verbose
    )
    rownames(meta_df) <- cn
    
    cols_to_get <- unique(c(sample_col, condition_col, timepoint_col, participant_col))
    cols_to_get <- intersect(cols_to_get, colnames(obj@meta.data))
    
    if (length(cols_to_get) > 0) {

      orig_meta <- unique(obj@meta.data[, cols_to_get, drop = FALSE])
      
      if (sample_col %in% colnames(orig_meta) && sample_col %in% colnames(meta_df)) {

        orig_meta[[sample_col]] <- as.character(orig_meta[[sample_col]])
        
        orig_meta <- orig_meta[!duplicated(orig_meta[[sample_col]]), , drop = FALSE]
        rownames(orig_meta) <- orig_meta[[sample_col]]
        
        sample_ids <- as.character(meta_df[[sample_col]])
        
        found_ids <- sample_ids %in% rownames(orig_meta)
        
        if (verbose && sum(found_ids) == 0) {
          message("      WARNING: No sample IDs matched in metadata lookup")
          message("        Parsed IDs: ", paste(head(sample_ids, 5), collapse = ", "))
          message("        Original IDs: ", paste(head(rownames(orig_meta), 5), collapse = ", "))
        }
        
        for (col in setdiff(cols_to_get, colnames(meta_df))) {
          meta_df[[col]] <- orig_meta[sample_ids, col]
        }
      }
    }
    
    meta_df
  }
  
  #' Generate pairwise contrasts dynamically from available groups
  generate_contrasts <- function(groups, design_levels, strategy,
                                 ref_condition = NULL, timepoints = NULL) {
    contrasts_list <- list()
    
    if (strategy == "stratified") {

      if (!is.null(ref_condition)) {

        other_groups <- setdiff(groups, ref_condition)
        for (g in other_groups) {
          name <- paste0(g, "_vs_", ref_condition)
          formula <- paste0("`", g, "` - `", ref_condition, "`")
          contrasts_list[[name]] <- formula
        }
      } else {

        if (length(groups) >= 2) {
          pairs <- utils::combn(groups, 2, simplify = FALSE)
          for (p in pairs) {
            name <- paste0(p[1], "_vs_", p[2])
            formula <- paste0("`", p[1], "` - `", p[2], "`")
            contrasts_list[[name]] <- formula
          }
        }
      }
      
    } else {

      if (!is.null(timepoints)) {
        conditions <- unique(gsub("^[^_]+_", "", groups))
        
        for (tp in timepoints) {
          tp_groups <- groups[grepl(paste0("^", tp, "_"), groups)]
          tp_conditions <- gsub(paste0("^", tp, "_"), "", tp_groups)
          
          if (!is.null(ref_condition) && ref_condition %in% tp_conditions) {

            other_conds <- setdiff(tp_conditions, ref_condition)
            for (cond in other_conds) {
              g1 <- paste0(tp, "_", cond)
              g2 <- paste0(tp, "_", ref_condition)
              if (g1 %in% groups && g2 %in% groups) {
                name <- paste0(cond, "_vs_", ref_condition, "_at_", tp)
                formula <- paste0("`", g1, "` - `", g2, "`")
                contrasts_list[[name]] <- formula
              }
            }
          } else {

            if (length(tp_conditions) >= 2) {
              pairs <- utils::combn(tp_conditions, 2, simplify = FALSE)
              for (p in pairs) {
                g1 <- paste0(tp, "_", p[1])
                g2 <- paste0(tp, "_", p[2])
                if (g1 %in% groups && g2 %in% groups) {
                  name <- paste0(p[1], "_vs_", p[2], "_at_", tp)
                  formula <- paste0("`", g1, "` - `", g2, "`")
                  contrasts_list[[name]] <- formula
                }
              }
            }
          }
        }
      }
    }
    
    contrasts_list
  }
  
  #' Filter contrasts to only those with valid groups in design
  filter_valid_contrasts <- function(contrasts_list, design_levels) {
    Filter(function(formula) {

      cleaned <- gsub("[`]", "", formula)
      cleaned <- gsub("\\s*-\\s*", " ", cleaned)
      cleaned <- gsub("\\s*\\+\\s*", " ", cleaned)
      groups_in_formula <- unlist(strsplit(trimws(cleaned), "\\s+"))
      groups_in_formula <- groups_in_formula[groups_in_formula != ""]
      all(groups_in_formula %in% design_levels)
    }, contrasts_list)
  }
  
  #' Fit limma-voom on a pseudobulk counts matrix and test the contrasts
  #'
  #' @param counts_mat Pseudobulk counts (genes x units).
  #' @param meta_df Per-unit metadata aligned to `counts_mat` columns.
  #' @param design Model matrix.
  #' @param contrast_matrix Contrast matrix from `limma::makeContrasts`.
  #' @param block_col Optional metadata column for `duplicateCorrelation()`
  #'   blocking (repeated measures).
  #' @param cluster_name Cluster label used in progress messages and results.
  #' @return Data frame of per-gene results for every contrast, or `NULL`.
  run_limma_voom <- function(counts_mat, meta_df, design, contrast_matrix,
                             block_col = NULL, cluster_name = NULL) {
    if (verbose && !is.null(cluster_name)) {
      message("  Running limma-voom for: ", cluster_name)
    }
    
    dge <- edgeR::DGEList(counts = counts_mat)
    
    group_var <- if ("Group" %in% colnames(meta_df)) {
      meta_df$Group
    } else if ("condition_valid" %in% colnames(meta_df)) {
      meta_df$condition_valid
    } else if (condition_col %in% colnames(meta_df)) {
      meta_df[[condition_col]]
    } else {
      NULL
    }
    
    keep <- edgeR::filterByExpr(dge, group = group_var)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (nrow(dge) == 0) {
      if (verbose) message("    No genes passed filtering")
      return(NULL)
    }
    
    if (verbose) message("    Genes after filtering: ", nrow(dge))
    
    dge <- edgeR::calcNormFactors(dge)
    
    v <- limma::voom(dge, design, plot = FALSE)
    
    if (!is.null(block_col) && block_col %in% colnames(meta_df) && block_repeated_measures) {
      block_values <- meta_df[[block_col]]
      n_unique_blocks <- length(unique(block_values))
      
      if (n_unique_blocks > 1 && n_unique_blocks < nrow(meta_df)) {
        if (verbose) message("    Applying duplicateCorrelation (", n_unique_blocks, " blocks)")
        corfit <- limma::duplicateCorrelation(v, design, block = block_values)
        
        if (!is.na(corfit$consensus) && corfit$consensus > -1 && corfit$consensus < 1) {
          fit <- limma::lmFit(v, design, block = block_values,
                              correlation = corfit$consensus)
          if (verbose) message("    Consensus correlation: ", round(corfit$consensus, 3))
        } else {
          if (verbose) message("    Invalid correlation estimate, fitting without blocking")
          fit <- limma::lmFit(v, design)
        }
      } else {
        fit <- limma::lmFit(v, design)
      }
    } else {
      fit <- limma::lmFit(v, design)
    }
    
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    
    results_list <- lapply(colnames(contrast_matrix), function(coef) {
      tt <- limma::topTable(fit2, coef = coef, number = Inf,
                            adjust.method = p_adjust_method)
      tt$gene <- rownames(tt)
      tt$contrast <- coef
      tt
    })
    names(results_list) <- colnames(contrast_matrix)
    
    if (verbose) {
      for (coef in names(results_list)) {
        n_sig <- sum(results_list[[coef]]$adj.P.Val < 0.05, na.rm = TRUE)
        message("    ", coef, ": ", n_sig, " DEGs (adj.P < 0.05)")
      }
    }
    
    if (return_fit) {
      return(list(results = results_list, fit = fit2, voom = v))
    }
    results_list
  }
  
  #' Fit edgeR quasi-likelihood models on a pseudobulk counts matrix
  #'
  #' @param counts_mat Pseudobulk counts (genes x units).
  #' @param meta_df Per-unit metadata aligned to `counts_mat` columns.
  #' @param design Model matrix.
  #' @param contrast_matrix Contrast matrix from `limma::makeContrasts`.
  #' @param cluster_name Cluster label used in progress messages and results.
  #' @return Data frame of per-gene results for every contrast, or `NULL`.
  run_edger <- function(counts_mat, meta_df, design, contrast_matrix,
                        cluster_name = NULL) {
    if (verbose && !is.null(cluster_name)) {
      message("  Running edgeR for: ", cluster_name)
    }
    
    dge <- edgeR::DGEList(counts = counts_mat)
    
    group_var <- if ("Group" %in% colnames(meta_df)) {
      meta_df$Group
    } else if ("condition_valid" %in% colnames(meta_df)) {
      meta_df$condition_valid
    } else if (condition_col %in% colnames(meta_df)) {
      meta_df[[condition_col]]
    } else {
      NULL
    }
    
    keep <- edgeR::filterByExpr(dge, group = group_var)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (nrow(dge) == 0) {
      if (verbose) message("    No genes passed filtering")
      return(NULL)
    }
    
    if (verbose) message("    Genes after filtering: ", nrow(dge))
    
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design)
    
    fit <- edgeR::glmQLFit(dge, design)
    
    results_list <- lapply(colnames(contrast_matrix), function(coef_name) {
      contrast_vec <- contrast_matrix[, coef_name]
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
      tt <- edgeR::topTags(qlf, n = Inf, adjust.method = p_adjust_method)$table
      tt$gene <- rownames(tt)
      tt$contrast <- coef_name
      tt
    })
    names(results_list) <- colnames(contrast_matrix)
    
    if (verbose) {
      for (coef in names(results_list)) {
        n_sig <- sum(results_list[[coef]]$FDR < 0.05, na.rm = TRUE)
        message("    ", coef, ": ", n_sig, " DEGs (FDR < 0.05)")
      }
    }
    
    if (return_fit) {
      return(list(results = results_list, fit = fit, dge = dge))
    }
    results_list
  }
  
  #' Fit DESeq2 on a pseudobulk counts matrix and test the contrasts
  #'
  #' @param counts_mat Pseudobulk counts (genes x units).
  #' @param meta_df Per-unit metadata aligned to `counts_mat` columns.
  #' @param design_col Metadata column holding the condition factor.
  #' @param contrasts_list Named list of contrast formulas to evaluate.
  #' @param cluster_name Cluster label used in progress messages and results.
  #' @return Data frame of per-gene results for every contrast, or `NULL`.
  run_deseq2 <- function(counts_mat, meta_df, design_col, contrasts_list,
                         cluster_name = NULL) {
    if (verbose && !is.null(cluster_name)) {
      message("  Running DESeq2 for: ", cluster_name)
    }
    
    counts_mat <- round(counts_mat)
    mode(counts_mat) <- "integer"
    
    design_formula <- stats::as.formula(paste("~", design_col))
    
    meta_df[[design_col]] <- factor(meta_df[[design_col]])
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_mat,
      colData = meta_df,
      design = design_formula
    )
    
    keep <- rowSums(DESeq2::counts(dds) >= 10) >= min_samples_per_group
    dds <- dds[keep, ]
    
    if (nrow(dds) == 0) {
      if (verbose) message("    No genes passed filtering")
      return(NULL)
    }
    
    if (verbose) message("    Genes after filtering: ", nrow(dds))
    
    dds <- DESeq2::DESeq(dds, quiet = !verbose)
    
    results_list <- lapply(names(contrasts_list), function(contrast_name) {
      contrast_spec <- contrasts_list[[contrast_name]]
      
      res <- tryCatch({
        DESeq2::results(dds, contrast = contrast_spec,
                        pAdjustMethod = p_adjust_method)
      }, error = function(e) {
        if (verbose) message("    Error extracting results for ", contrast_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(res)) return(NULL)
      
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df$contrast <- contrast_name
      res_df
    })
    
    results_list <- Filter(Negate(is.null), results_list)
    names(results_list) <- names(contrasts_list)[seq_along(results_list)]
    
    if (verbose) {
      for (coef in names(results_list)) {
        n_sig <- sum(results_list[[coef]]$padj < 0.05, na.rm = TRUE)
        message("    ", coef, ": ", n_sig, " DEGs (padj < 0.05)")
      }
    }
    
    if (return_fit) {
      return(list(results = results_list, dds = dds))
    }
    results_list
  }
  
  #' Two-group pseudobulk comparison for one cluster (`ident.1` vs `ident.2`)
  #'
  #' @param obj_cluster Seurat object subset to one cluster's cells.
  #' @param cluster_name The cluster's label, carried into the results.
  #' @return Data frame of per-gene results, or `NULL` when the comparison
  #'   cannot be formed (group sizes, missing levels).
  run_direct_comparison <- function(obj_cluster, cluster_name) {
    if (verbose) message("  Aggregating by: ", condition_col, " x ", sample_col)
    
    counts_mat <- tryCatch({
      ag <- Seurat::AggregateExpression(
        object = obj_cluster,
        assays = pseudobulk_assay,
        group.by = c(condition_col, sample_col),
        slot = "counts",
        return.seurat = FALSE,
        verbose = FALSE
      )[[pseudobulk_assay]]
      as.matrix(ag)
    }, error = function(e) {
      if (verbose) message("  Error in aggregation: ", e$message)
      NULL
    })
    
    if (is.null(counts_mat) || ncol(counts_mat) == 0) {
      if (verbose) message("  Skipping: no pseudobulk samples generated")
      return(NULL)
    }
    
    if (verbose) message("  Samples aggregated: ", ncol(counts_mat))
    
    col_names <- colnames(counts_mat)
    
    meta_df <- data.frame(
      sample_id = col_names,
      row.names = col_names,
      stringsAsFactors = FALSE
    )
    
    mapped <- .map_aggregated_columns(
      col_names  = col_names,
      meta       = obj_cluster@meta.data,
      group_cols = c(condition_col, sample_col),
      verbose    = verbose
    )
    meta_df[[condition_col]] <- mapped[[condition_col]]
    meta_df[[sample_col]]    <- mapped[[sample_col]]
    
    if (participant_col != sample_col) {

      orig_meta <- unique(obj_cluster@meta.data[, c(sample_col, participant_col), drop = FALSE])
      orig_meta[[sample_col]] <- as.character(orig_meta[[sample_col]])
      orig_meta <- orig_meta[!duplicated(orig_meta[[sample_col]]), ]
      rownames(orig_meta) <- orig_meta[[sample_col]]
      meta_df[[participant_col]] <- orig_meta[meta_df[[sample_col]], participant_col]
    } else {
      meta_df[[participant_col]] <- meta_df[[sample_col]]
    }
    
    if (verbose) {
      message("  Conditions found: ", paste(unique(meta_df[[condition_col]]), collapse = ", "))
      message("  Unique participants: ", length(unique(meta_df[[participant_col]])))
      
      cond_counts <- table(meta_df[[condition_col]])
      message("  Samples per condition: ", paste(names(cond_counts), "=", cond_counts, collapse = ", "))
    }
    
    if (!is.null(conditions_to_analyze)) {
      keep <- meta_df[[condition_col]] %in% as.character(conditions_to_analyze)
      meta_df <- meta_df[keep, , drop = FALSE]
      counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
    }
    
    if (nrow(meta_df) == 0) {
      if (verbose) message("  Skipping: no samples match conditions_to_analyze")
      return(NULL)
    }
    
    group_counts <- table(meta_df[[condition_col]])
    valid_groups <- names(group_counts)[group_counts >= min_samples_per_group]
    
    if (verbose) {
      message("  Group counts: ", paste(names(group_counts), "=", group_counts, collapse = ", "))
      message("  Valid groups (n >= ", min_samples_per_group, "): ", paste(valid_groups, collapse = ", "))
    }
    
    if (length(valid_groups) < 2) {
      if (verbose) message("  Skipping: < 2 groups with min samples")
      return(NULL)
    }
    
    keep <- meta_df[[condition_col]] %in% valid_groups
    meta_df <- meta_df[keep, , drop = FALSE]
    counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
    
    original_groups <- valid_groups
    valid_names <- make.names(valid_groups)
    names_map <- setNames(original_groups, valid_names)
    reverse_names_map <- setNames(valid_names, original_groups)
    
    meta_df$condition_valid <- reverse_names_map[as.character(meta_df[[condition_col]])]
    meta_df$condition_valid <- factor(meta_df$condition_valid, levels = valid_names)
    
    if (verbose) {
      message("  Groups: ", paste(original_groups, "=", table(meta_df[[condition_col]]), collapse = ", "))
      if (any(valid_names != original_groups)) {
        message("  (Using R-safe names: ", paste(original_groups, "->", valid_names, collapse = ", "), ")")
      }
    }
    
    design <- stats::model.matrix(~ 0 + condition_valid, data = meta_df)
    colnames(design) <- gsub("^condition_valid", "", colnames(design))
    
    if (!is.null(ident.1) && !is.null(ident.2)) {
      ident.1_valid <- reverse_names_map[as.character(ident.1)]
      ident.2_valid <- reverse_names_map[as.character(ident.2)]
      if (!is.null(ident.1_valid) && !is.null(ident.2_valid)) {
        contrast_specs <- list()
        contrast_specs[[paste0(ident.1, "_vs_", ident.2)]] <-
          paste0(ident.1_valid, " - ", ident.2_valid)
      } else {
        if (verbose) message("  Skipping: ident.1 or ident.2 not in valid groups")
        return(NULL)
      }
    } else if (!is.null(contrasts)) {
      contrast_specs <- lapply(contrasts, function(formula) {
        .substitute_group_names(formula, original_groups, reverse_names_map)
      })
    } else {
      contrast_specs <- list()
      if (!is.null(reference_condition) && as.character(reference_condition) %in% original_groups) {
        ref_valid <- reverse_names_map[as.character(reference_condition)]
        other_groups <- setdiff(original_groups, as.character(reference_condition))
        for (g in other_groups) {
          g_valid <- reverse_names_map[g]
          name <- paste0(g, "_vs_", reference_condition)
          contrast_specs[[name]] <- paste0(g_valid, " - ", ref_valid)
        }
      } else {
        pairs <- utils::combn(original_groups, 2, simplify = FALSE)
        for (p in pairs) {
          name <- paste0(p[1], "_vs_", p[2])
          p1_valid <- reverse_names_map[p[1]]
          p2_valid <- reverse_names_map[p[2]]
          contrast_specs[[name]] <- paste0(p1_valid, " - ", p2_valid)
        }
      }
    }
    
    if (length(contrast_specs) == 0) {
      if (verbose) message("  Skipping: no valid contrasts")
      return(NULL)
    }
    
    if (verbose) {
      message("  Contrasts to test: ", paste(names(contrast_specs), collapse = ", "))
    }
    
    contrast_formulas <- unlist(contrast_specs)
    
    contrast_matrix <- tryCatch({
      limma::makeContrasts(contrasts = contrast_formulas, levels = design)
    }, error = function(e) {
      if (verbose) message("  Error building contrasts: ", e$message)
      NULL
    })
    
    if (is.null(contrast_matrix)) return(NULL)
    
    participant_counts <- table(meta_df[[participant_col]])
    has_repeated_measures <- any(participant_counts > 1)
    
    if (verbose && has_repeated_measures) {
      message("  Detected repeated measures (", sum(participant_counts > 1), 
              " participants with multiple conditions)")
    }
    
    if (has_repeated_measures && !identical(test.use, "limma")) {
      warning(sprintf(paste0(
        "test.use = '%s' does not model repeated measures. %d of %d participants ",
        "contribute more than one sample to this comparison, so specimens are ",
        "being treated as independent replicates and p-values will be ",
        "anti-conservative. Use test.use = 'limma' (which blocks via ",
        "duplicateCorrelation), or stratify so that each participant appears once ",
        "per comparison."),
        test.use, sum(participant_counts > 1), length(participant_counts)),
        call. = FALSE)
      if (verbose)
        message("  WARNING: ", test.use, " does not block on ", participant_col,
                " - see warning above.")
    }
    
    res <- tryCatch({
      if (test.use == "limma") {
        block_col_to_use <- if (has_repeated_measures && block_repeated_measures) {
          participant_col
        } else {
          NULL
        }
        run_limma_voom(counts_mat, meta_df, design, contrast_matrix,
                       block_col = block_col_to_use, cluster_name = cluster_name)
      } else if (test.use == "edger") {
        run_edger(counts_mat, meta_df, design, contrast_matrix,
                  cluster_name = cluster_name)
      } else if (test.use == "deseq2") {
        deseq_contrasts <- lapply(names(contrast_specs), function(cn) {
          formula <- contrast_specs[[cn]]
          groups <- trimws(unlist(strsplit(formula, " - ")))
          c("condition_valid", groups[1], groups[2])
        })
        names(deseq_contrasts) <- names(contrast_specs)
        run_deseq2(counts_mat, meta_df, "condition_valid", deseq_contrasts,
                   cluster_name = cluster_name)
      }
    }, error = function(e) {
      if (verbose) message("  Error in DE analysis: ", e$message)
      NULL
    })
    
    if (!is.null(res) && is.list(res) && !is.data.frame(res)) {
      all_dfs <- lapply(names(res), function(contrast_name) {
        df <- res[[contrast_name]]
        if (is.data.frame(df)) {
          df$cluster <- cluster_name
        }
        df
      })
      res <- dplyr::bind_rows(all_dfs)
    } else if (is.data.frame(res)) {
      res$cluster <- cluster_name
    }
    
    res
  }
  
  #' Stratified pseudobulk design for one cluster
  #'
  #' Fits a separate condition model per timepoint (the `"stratified"` `design_strategy`), aggregating per sample within the cluster.
  #'
  #' @param obj Seurat object (all clusters; subset internally).
  #' @param cluster_name The cluster to analyze.
  #' @return Data frame of per-gene results across timepoints and contrasts,
  #'   or `NULL` when every stratum was skipped.
  run_stratified_pseudobulk <- function(obj, cluster_name) {
    if (verbose) message("\n--- Cluster: ", cluster_name, " (stratified design) ---")
    
    cluster_cells <- obj@meta.data[[cluster_col]] == cluster_name
    n_cells <- sum(cluster_cells, na.rm = TRUE)
    
    if (n_cells < min_cells_per_pseudobulk) {
      if (verbose) message("  Skipping: insufficient cells (", n_cells, ")")
      return(NULL)
    }
    
    obj_cluster <- subset(obj, cells = colnames(obj)[cluster_cells])
    
    if (!is.null(conditions_to_analyze)) {
      cond_cells <- obj_cluster@meta.data[[condition_col]] %in% as.character(conditions_to_analyze)
      cond_cells[is.na(cond_cells)] <- FALSE
      if (sum(cond_cells) < min_cells_per_pseudobulk) {
        if (verbose) message("  Skipping: insufficient cells after condition filter (", sum(cond_cells), ")")
        return(NULL)
      }
      obj_cluster <- subset(obj_cluster, cells = colnames(obj_cluster)[cond_cells])
      if (verbose) {
        message("  Filtered to conditions: ", paste(conditions_to_analyze, collapse = ", "),
                " (", ncol(obj_cluster), " cells)")

        cond_table <- table(obj_cluster@meta.data[[condition_col]])
        message("  Condition distribution after filter: ", 
                paste(names(cond_table), "=", cond_table, collapse = ", "))
        n_participants <- length(unique(obj_cluster@meta.data[[sample_col]]))
        message("  Unique participants after filter: ", n_participants)
      }
    }
    
    direct_mode <- is.null(timepoints_to_analyze) || (condition_col == timepoint_col)
    
    if (direct_mode) {
      if (verbose) message("  Running in DIRECT COMPARISON mode (no timepoint stratification)")
      return(run_direct_comparison(obj_cluster, cluster_name))
    }
    
    available_tp <- unique(obj_cluster@meta.data[[timepoint_col]])
    available_tp <- available_tp[!is.na(available_tp)]
    
    if (!is.null(timepoints_to_analyze)) {
      available_tp <- intersect(as.character(available_tp),
                                as.character(timepoints_to_analyze))
    }
    
    if (length(available_tp) == 0) {
      if (verbose) message("  Skipping: no matching timepoints")
      return(NULL)
    }
    
    if (verbose) message("  Timepoints to analyze: ", paste(available_tp, collapse = ", "))
    
    results_by_timepoint <- list()
    
    for (tp in available_tp) {
      if (verbose) message("\n  Timepoint: ", tp)
      
      tp_cells <- as.character(obj_cluster@meta.data[[timepoint_col]]) == as.character(tp)
      tp_cells[is.na(tp_cells)] <- FALSE
      
      if (sum(tp_cells) < min_cells_per_pseudobulk) {
        if (verbose) message("    Skipping: insufficient cells (", sum(tp_cells), ")")
        next
      }
      
      obj_tp <- subset(obj_cluster, cells = colnames(obj_cluster)[tp_cells])
      
      if (verbose) {
        cond_at_tp <- table(obj_tp@meta.data[[condition_col]])
        message("    Cells: ", ncol(obj_tp))
        message("    Conditions at this timepoint: ", 
                paste(names(cond_at_tp), "=", cond_at_tp, collapse = ", "))
        samples_at_tp <- unique(obj_tp@meta.data[[sample_col]])
        message("    Unique samples: ", length(samples_at_tp), 
                " (", paste(head(samples_at_tp, 6), collapse = ", "), 
                if(length(samples_at_tp) > 6) "..." else "", ")")
      }
      
      counts_mat <- aggregate_pseudobulk(obj_tp, c(sample_col))
      if (is.null(counts_mat) || ncol(counts_mat) == 0) {
        if (verbose) message("    Skipping: no pseudobulk samples generated")
        next
      }
      
      meta_df <- build_pseudobulk_meta(obj_tp, counts_mat, c(sample_col))
      meta_df <- meta_df[colnames(counts_mat), , drop = FALSE]
      
      if (verbose) {
        message("    Samples aggregated: ", nrow(meta_df))
        message("    Sample IDs: ", paste(meta_df[[sample_col]], collapse = ", "))
        if (condition_col %in% colnames(meta_df)) {
          message("    Conditions found: ", paste(unique(meta_df[[condition_col]]), collapse = ", "))
          message("    Condition values: ", paste(meta_df[[condition_col]], collapse = ", "))
        } else {
          message("    WARNING: condition_col '", condition_col, "' not in meta_df columns: ",
                  paste(colnames(meta_df), collapse = ", "))
        }
      }
      
      if (!condition_col %in% colnames(meta_df) || all(is.na(meta_df[[condition_col]]))) {
        if (verbose) message("    Attempting direct condition lookup from original metadata...")
        
        orig_sample_cond <- unique(obj_tp@meta.data[, c(sample_col, condition_col), drop = FALSE])
        rownames(orig_sample_cond) <- as.character(orig_sample_cond[[sample_col]])
        
        sample_ids <- as.character(meta_df[[sample_col]])
        meta_df[[condition_col]] <- orig_sample_cond[sample_ids, condition_col]
        
        if (verbose) {
          message("    After direct lookup - Conditions: ", 
                  paste(meta_df[[condition_col]], collapse = ", "))
        }
      }
      
      valid_cond <- meta_df[[condition_col]]
      valid_cond <- valid_cond[!is.na(valid_cond)]
      valid_cond <- unique(valid_cond)
      
      if (!is.null(conditions_to_analyze)) {
        valid_cond <- intersect(as.character(valid_cond), as.character(conditions_to_analyze))
      }
      
      if (verbose) {
        message("    Valid conditions after filtering: ", paste(valid_cond, collapse = ", "))
      }
      
      keep_samples <- !is.na(meta_df[[condition_col]]) &
        meta_df[[condition_col]] %in% valid_cond
      meta_df <- meta_df[keep_samples, , drop = FALSE]
      counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
      
      if (nrow(meta_df) == 0) {
        if (verbose) message("    Skipping: no samples match conditions_to_analyze")
        next
      }
      
      group_counts <- table(meta_df[[condition_col]])
      valid_groups <- names(group_counts)[group_counts >= min_samples_per_group]
      
      if (verbose) {
        message("    Group counts: ", paste(names(group_counts), "=", group_counts, collapse = ", "))
        message("    Valid groups (n >= ", min_samples_per_group, "): ", paste(valid_groups, collapse = ", "))
      }
      
      if (length(valid_groups) < 2) {
        if (verbose) message("    Skipping: < 2 groups with min samples")
        next
      }
      
      keep_samples <- meta_df[[condition_col]] %in% valid_groups
      meta_df <- meta_df[keep_samples, , drop = FALSE]
      counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
      
      original_groups <- valid_groups
      valid_names <- make.names(valid_groups)
      names_map <- setNames(original_groups, valid_names)
      reverse_names_map <- setNames(valid_names, original_groups)
      
      meta_df$condition_valid <- reverse_names_map[as.character(meta_df[[condition_col]])]
      meta_df$condition_valid <- factor(meta_df$condition_valid, levels = valid_names)
      
      if (verbose) {
        message("    Groups: ", paste(original_groups, "=", table(meta_df[[condition_col]]), collapse = ", "))
        if (any(valid_names != original_groups)) {
          message("    (Using R-safe names: ", paste(original_groups, "->", valid_names, collapse = ", "), ")")
        }
      }
      
      design <- stats::model.matrix(~ 0 + condition_valid, data = meta_df)
      colnames(design) <- gsub("^condition_valid", "", colnames(design))
      
      if (!is.null(ident.1) && !is.null(ident.2)) {

        ident.1_valid <- reverse_names_map[as.character(ident.1)]
        ident.2_valid <- reverse_names_map[as.character(ident.2)]
        if (!is.null(ident.1_valid) && !is.null(ident.2_valid)) {
          contrast_specs <- list()

          contrast_specs[[paste0(ident.1, "_vs_", ident.2)]] <-
            paste0(ident.1_valid, " - ", ident.2_valid)
        } else {
          if (verbose) message("    Skipping: ident.1 or ident.2 not in valid groups")
          next
        }
      } else if (!is.null(contrasts)) {

        contrast_specs <- lapply(contrasts, function(formula) {
          .substitute_group_names(formula, original_groups, reverse_names_map)
        })
      } else {

        contrast_specs <- list()
        if (!is.null(reference_condition)) {
          ref_valid <- reverse_names_map[as.character(reference_condition)]
          other_groups <- setdiff(original_groups, as.character(reference_condition))
          for (g in other_groups) {
            g_valid <- reverse_names_map[g]
            name <- paste0(g, "_vs_", reference_condition)
            contrast_specs[[name]] <- paste0(g_valid, " - ", ref_valid)
          }
        } else {
          pairs <- utils::combn(original_groups, 2, simplify = FALSE)
          for (p in pairs) {
            name <- paste0(p[1], "_vs_", p[2])
            p1_valid <- reverse_names_map[p[1]]
            p2_valid <- reverse_names_map[p[2]]
            contrast_specs[[name]] <- paste0(p1_valid, " - ", p2_valid)
          }
        }
      }
      
      if (length(contrast_specs) == 0) {
        if (verbose) message("    Skipping: no valid contrasts")
        next
      }
      
      if (verbose) {
        message("    Contrasts to test: ", paste(names(contrast_specs), collapse = ", "))
      }
      
      contrast_formulas <- unlist(contrast_specs)
      
      contrast_matrix <- tryCatch({
        limma::makeContrasts(contrasts = contrast_formulas, levels = design)
      }, error = function(e) {
        if (verbose) message("    Error building contrasts: ", e$message)
        NULL
      })
      
      if (is.null(contrast_matrix)) next
      
      res <- tryCatch({
        if (test.use == "limma") {
          run_limma_voom(counts_mat, meta_df, design, contrast_matrix,
                         cluster_name = paste(cluster_name, tp, sep = "_tp"))
        } else if (test.use == "edger") {
          run_edger(counts_mat, meta_df, design, contrast_matrix,
                    cluster_name = paste(cluster_name, tp, sep = "_tp"))
        } else if (test.use == "deseq2") {

          deseq_contrasts <- lapply(names(contrast_specs), function(cn) {

            formula <- contrast_specs[[cn]]
            groups <- trimws(unlist(strsplit(formula, " - ")))
            c("condition_valid", groups[1], groups[2])
          })
          names(deseq_contrasts) <- names(contrast_specs)
          run_deseq2(counts_mat, meta_df, "condition_valid", deseq_contrasts,
                     cluster_name = paste(cluster_name, tp, sep = "_tp"))
        }
      }, error = function(e) {
        if (verbose) message("    Error in DE analysis: ", e$message)
        NULL
      })
      
      if (!is.null(res)) {

        if (is.list(res) && !is.data.frame(res)) {
          res <- lapply(res, function(df) {
            if (is.data.frame(df)) {
              df$timepoint <- tp
              df$cluster <- cluster_name
            }
            df
          })
        }
        results_by_timepoint[[as.character(tp)]] <- res
      }
    }
    
    if (length(results_by_timepoint) == 0) return(NULL)
    
    all_dfs <- list()
    for (tp_name in names(results_by_timepoint)) {
      tp_res <- results_by_timepoint[[tp_name]]
      if (is.list(tp_res) && !is.data.frame(tp_res)) {
        for (contrast_name in names(tp_res)) {
          df <- tp_res[[contrast_name]]
          if (is.data.frame(df)) {
            all_dfs[[length(all_dfs) + 1]] <- df
          }
        }
      } else if (is.data.frame(tp_res)) {
        all_dfs[[length(all_dfs) + 1]] <- tp_res
      }
    }
    
    if (length(all_dfs) == 0) return(NULL)
    dplyr::bind_rows(all_dfs)
  }
  
  #' Interaction pseudobulk design for one cluster
  #'
  #' Single model with a timepoint-by-condition interaction term (the `"interaction"` `design_strategy`), with repeated measures handled via `duplicateCorrelation()` when `block_repeated_measures` is `TRUE`.
  #'
  #' @param obj Seurat object (all clusters; subset internally).
  #' @param cluster_name The cluster to analyze.
  #' @return Data frame of per-gene results, or `NULL`.
  run_interaction_pseudobulk <- function(obj, cluster_name) {
    if (verbose) message("\n--- Cluster: ", cluster_name, " (interaction design) ---")
    
    cluster_cells <- obj@meta.data[[cluster_col]] == cluster_name
    n_cells <- sum(cluster_cells, na.rm = TRUE)
    
    if (n_cells < min_cells_per_pseudobulk) {
      if (verbose) message("  Skipping: insufficient cells (", n_cells, ")")
      return(NULL)
    }
    
    obj_cluster <- subset(obj, cells = colnames(obj)[cluster_cells])
    
    if (!is.null(timepoints_to_analyze)) {
      tp_cells <- as.character(obj_cluster@meta.data[[timepoint_col]]) %in%
        as.character(timepoints_to_analyze)
      tp_cells[is.na(tp_cells)] <- FALSE
      obj_cluster <- subset(obj_cluster, cells = colnames(obj_cluster)[tp_cells])
    }
    
    if (!is.null(conditions_to_analyze)) {
      cond_cells <- obj_cluster@meta.data[[condition_col]] %in% conditions_to_analyze
      cond_cells[is.na(cond_cells)] <- FALSE
      obj_cluster <- subset(obj_cluster, cells = colnames(obj_cluster)[cond_cells])
    }
    
    if (ncol(obj_cluster) < min_cells_per_pseudobulk) {
      if (verbose) message("  Skipping: insufficient cells after filtering")
      return(NULL)
    }
    
    counts_mat <- aggregate_pseudobulk(obj_cluster, c(timepoint_col, sample_col))
    if (is.null(counts_mat) || ncol(counts_mat) == 0) {
      if (verbose) message("  Skipping: no pseudobulk samples")
      return(NULL)
    }
    
    meta_df <- build_pseudobulk_meta(obj_cluster, counts_mat, c(timepoint_col, sample_col))
    meta_df <- meta_df[colnames(counts_mat), , drop = FALSE]
    
    complete_rows <- !is.na(meta_df[[timepoint_col]]) & !is.na(meta_df[[condition_col]])
    meta_df <- meta_df[complete_rows, , drop = FALSE]
    counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
    
    if (nrow(meta_df) == 0) {
      if (verbose) message("  Skipping: no complete samples")
      return(NULL)
    }
    
    meta_df$Group <- paste(meta_df[[timepoint_col]], meta_df[[condition_col]], sep = "_")
    meta_df$Group <- factor(meta_df$Group)
    
    group_counts <- table(meta_df$Group)
    valid_groups <- names(group_counts)[group_counts >= min_samples_per_group]
    
    if (length(valid_groups) < 2) {
      if (verbose) message("  Skipping: < 2 groups with min samples")
      return(NULL)
    }
    
    keep_samples <- meta_df$Group %in% valid_groups
    meta_df <- meta_df[keep_samples, , drop = FALSE]
    counts_mat <- counts_mat[, rownames(meta_df), drop = FALSE]
    
    original_groups <- valid_groups
    valid_names <- make.names(valid_groups)
    names_map <- setNames(original_groups, valid_names)
    reverse_names_map <- setNames(valid_names, original_groups)
    
    meta_df$Group_valid <- reverse_names_map[as.character(meta_df$Group)]
    meta_df$Group_valid <- factor(meta_df$Group_valid, levels = valid_names)
    
    if (verbose) {
      message("  Groups: ", paste(original_groups, "=", table(meta_df$Group), collapse = ", "))
      if (any(valid_names != original_groups)) {
        message("  (Using R-safe names internally)")
      }
    }
    
    design <- stats::model.matrix(~ 0 + Group_valid, data = meta_df)
    colnames(design) <- gsub("^Group_valid", "", colnames(design))
    
    valid_timepoints <- unique(gsub("_.*", "", original_groups))
    
    if (!is.null(ident.1) && !is.null(ident.2)) {

      contrast_specs <- list()
      for (tp in valid_timepoints) {
        g1_orig <- paste0(tp, "_", ident.1)
        g2_orig <- paste0(tp, "_", ident.2)
        if (g1_orig %in% original_groups && g2_orig %in% original_groups) {
          g1_valid <- reverse_names_map[g1_orig]
          g2_valid <- reverse_names_map[g2_orig]
          name <- paste0(ident.1, "_vs_", ident.2, "_at_", tp)
          contrast_specs[[name]] <- paste0(g1_valid, " - ", g2_valid)
        }
      }
    } else if (!is.null(contrasts)) {

      contrast_specs <- lapply(contrasts, function(formula) {
        .substitute_group_names(formula, original_groups, reverse_names_map)
      })
    } else {

      contrast_specs <- list()
      for (tp in valid_timepoints) {
        tp_groups_orig <- original_groups[grepl(paste0("^", tp, "_"), original_groups)]
        tp_conditions <- gsub(paste0("^", tp, "_"), "", tp_groups_orig)
        
        if (!is.null(reference_condition) && reference_condition %in% tp_conditions) {
          ref_orig <- paste0(tp, "_", reference_condition)
          ref_valid <- reverse_names_map[ref_orig]
          other_conds <- setdiff(tp_conditions, reference_condition)
          for (cond in other_conds) {
            g_orig <- paste0(tp, "_", cond)
            g_valid <- reverse_names_map[g_orig]
            name <- paste0(cond, "_vs_", reference_condition, "_at_", tp)
            contrast_specs[[name]] <- paste0(g_valid, " - ", ref_valid)
          }
        } else if (length(tp_conditions) >= 2) {
          pairs <- utils::combn(tp_conditions, 2, simplify = FALSE)
          for (p in pairs) {
            g1_orig <- paste0(tp, "_", p[1])
            g2_orig <- paste0(tp, "_", p[2])
            g1_valid <- reverse_names_map[g1_orig]
            g2_valid <- reverse_names_map[g2_orig]
            name <- paste0(p[1], "_vs_", p[2], "_at_", tp)
            contrast_specs[[name]] <- paste0(g1_valid, " - ", g2_valid)
          }
        }
      }
    }
    
    if (length(contrast_specs) == 0) {
      if (verbose) message("  Skipping: no valid contrasts for available groups")
      return(NULL)
    }
    
    if (verbose) {
      message("  Contrasts to test: ", paste(names(contrast_specs), collapse = ", "))
    }
    
    contrast_formulas <- unlist(contrast_specs)
    
    contrast_matrix <- tryCatch({
      limma::makeContrasts(contrasts = contrast_formulas, levels = design)
    }, error = function(e) {
      if (verbose) message("  Error building contrasts: ", e$message)
      NULL
    })
    
    if (is.null(contrast_matrix)) return(NULL)
    
    block_col <- if (block_repeated_measures) participant_col else NULL
    
    res <- tryCatch({
      if (test.use == "limma") {
        run_limma_voom(counts_mat, meta_df, design, contrast_matrix,
                       block_col = block_col, cluster_name = cluster_name)
      } else if (test.use == "edger") {
        run_edger(counts_mat, meta_df, design, contrast_matrix,
                  cluster_name = cluster_name)
      } else if (test.use == "deseq2") {

        deseq_contrasts <- lapply(names(contrast_specs), function(cn) {
          formula <- contrast_specs[[cn]]
          groups <- trimws(unlist(strsplit(formula, " - ")))
          c("Group_valid", groups[1], groups[2])
        })
        names(deseq_contrasts) <- names(contrast_specs)
        run_deseq2(counts_mat, meta_df, "Group_valid", deseq_contrasts,
                   cluster_name = cluster_name)
      }
    }, error = function(e) {
      if (verbose) message("  Error in DE analysis: ", e$message)
      NULL
    })
    
    if (!is.null(res) && is.list(res) && !is.data.frame(res)) {
      all_dfs <- lapply(names(res), function(contrast_name) {
        df <- res[[contrast_name]]
        if (is.data.frame(df)) {
          df$cluster <- cluster_name
        }
        df
      })
      res <- dplyr::bind_rows(all_dfs)
    }
    
    res
  }
  
  #' Seurat-native single-cell testing across clusters
  #'
  #' Runs `FindAllMarkers` (all-vs-rest per cluster) or, when `ident.1` is set, `FindMarkers` comparing it with `ident.2`, or with all other cells when `ident.2` is unset, using the selected single-cell test (`MAST`, `wilcox`, `bimod`, `roc`, `t`, `negbinom`, `poisson`, `LR`).
  #'
  #' @param obj Seurat object with the cluster column as active identity.
  #' @return Data frame of per-gene results.
  run_seurat_native <- function(obj) {

    if (is.null(ident.1) && is.null(ident.2) &&
        !is.null(cluster_col) && cluster_col %in% colnames(obj@meta.data)) {
      Seurat::Idents(obj) <- obj@meta.data[[cluster_col]]
    }

    if (!is.null(ident.1)) {
      wanted <- as.character(c(ident.1, ident.2))
      if (!all(wanted %in% levels(Seurat::Idents(obj)))) {
        for (id_col in c(condition_col, cluster_col)) {
          if (!is.null(id_col) && id_col %in% colnames(obj@meta.data) &&
              all(wanted %in% obj@meta.data[[id_col]])) {
            if (verbose) {
              message("Setting Idents from '", id_col,
                      "' (contains ident.1/ident.2 levels)")
            }
            Seurat::Idents(obj) <- obj@meta.data[[id_col]]
            break
          }
        }
      }
    }
    if (!is.null(latent.vars) && length(latent.vars) > 0L &&
        test.use %in% c("mast", "lr")) {
      if (verbose) {
        message("Using latent.vars: ", paste(latent.vars, collapse = ", "))
      }
    }
    
    test_map <- c("wilcox" = "wilcox", "bimod" = "bimod", "roc" = "roc",
                  "t" = "t", "negbinom" = "negbinom", "poisson" = "poisson",
                  "lr" = "LR", "mast" = "MAST")
    seurat_test <- test_map[test.use]
    
    if (is.null(ident.1) && is.null(ident.2)) {
      if (verbose) message("Running FindAllMarkers (test.use = ", seurat_test, ")")
      markers <- Seurat::FindAllMarkers(obj, test.use = seurat_test,
                                        latent.vars = latent.vars, ...)
    } else {
      if (verbose) {
        message("Running FindMarkers (test.use = ", seurat_test,
                ") comparing ", ident.1, " vs ",
                if (is.null(ident.2)) "rest" else ident.2)
      }
      markers <- Seurat::FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
                                     test.use = seurat_test,
                                     latent.vars = latent.vars, ...)
    }
    
    if (!("gene" %in% colnames(markers))) markers$gene <- rownames(markers)
    markers
  }
  
  all_results <- list()
  
  for (i in seq_along(objs)) {
    obj <- objs[[i]]
    obj_name <- obj_names[i]
    
    if (verbose) {
      message("\n", strrep("=", 60))
      message("PROCESSING: ", obj_name)
      message(strrep("=", 60))
    }
    
    if (is_bulk_method) {
      required_cols <- c(cluster_col, sample_col, condition_col)
      if (design_strategy == "interaction" || !is.null(timepoints_to_analyze)) {
        required_cols <- c(required_cols, timepoint_col)
      }
      if (design_strategy == "interaction" && block_repeated_measures) {
        required_cols <- c(required_cols, participant_col)
      }
    } else {
      required_cols <- c(cluster_col)
    }
    
    validate_columns(obj, unique(required_cols), obj_name)
    
    summarize_metadata(obj, unique(required_cols), obj_name)
    
    if (is_bulk_method) {
      validate_conditions(obj, condition_col, conditions_to_analyze, ident.1, ident.2)
      
      if (design_strategy == "interaction" || !is.null(timepoints_to_analyze)) {
        validate_timepoints(obj, timepoint_col, timepoints_to_analyze)
      }
      
      validate_clusters(obj, cluster_col)
      
      report_pseudobulk_structure(obj, cluster_col, sample_col, condition_col,
                                  if (design_strategy == "interaction") timepoint_col else NULL)
    }
    
    if (is_bulk_method) {
      clusters <- unique(obj@meta.data[[cluster_col]])
      clusters <- clusters[!is.na(clusters)]
      cluster_results <- list()
      
      for (cl in clusters) {
        res <- if (design_strategy == "stratified") {
          run_stratified_pseudobulk(obj, cl)
        } else {
          run_interaction_pseudobulk(obj, cl)
        }
        
        if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
          cluster_results[[as.character(cl)]] <- res
        }
      }
      
      if (length(cluster_results) > 0) {
        all_results[[obj_name]] <- dplyr::bind_rows(cluster_results)
      } else {
        all_results[[obj_name]] <- NULL
      }
      
    } else {
      all_results[[obj_name]] <- run_seurat_native(obj)
    }
  }
  
  final_dfs <- list()
  for (obj_name in names(all_results)) {
    res <- all_results[[obj_name]]
    if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
      if (length(objs) > 1) {
        res$object <- obj_name
      }
      final_dfs[[obj_name]] <- res
    }
  }
  
  if (length(final_dfs) > 0) {
    final_result <- dplyr::bind_rows(final_dfs)
    
    priority_cols <- c("gene", "cluster", "timepoint", "contrast", "logFC", "AveExpr",
                       "t", "P.Value", "adj.P.Val", "B",
                       "logCPM", "F", "PValue", "FDR",
                       "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                       "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"
    )
    existing_priority <- intersect(priority_cols, colnames(final_result))
    other_cols <- setdiff(colnames(final_result), priority_cols)
    final_result <- final_result[, c(existing_priority, other_cols), drop = FALSE]
    
  } else {
    final_result <- data.frame()
    if (verbose) message("No results produced.")
  }
  
  if (verbose && nrow(final_result) > 0) {
    message("\n", strrep("=", 60))
    message("ANALYSIS COMPLETE")
    message(strrep("=", 60))
    message("Total rows: ", format(nrow(final_result), big.mark = ","))
    message("Clusters analyzed: ", paste(unique(final_result$cluster), collapse = ", "))
    if ("timepoint" %in% colnames(final_result)) {
      message("Timepoints: ", paste(unique(final_result$timepoint), collapse = ", "))
    }
    if ("contrast" %in% colnames(final_result)) {
      message("Contrasts: ", paste(unique(final_result$contrast), collapse = ", "))
    }
    n_sig <- if ("adj.P.Val" %in% colnames(final_result)) {
      sum(final_result$adj.P.Val < 0.05, na.rm = TRUE)
    } else if ("FDR" %in% colnames(final_result)) {
      sum(final_result$FDR < 0.05, na.rm = TRUE)
    } else if ("padj" %in% colnames(final_result)) {
      sum(final_result$padj < 0.05, na.rm = TRUE)
    } else if ("p_val_adj" %in% colnames(final_result)) {
      sum(final_result$p_val_adj < 0.05, na.rm = TRUE)
    } else {
      NA
    }
    if (!is.na(n_sig)) {
      message("Significant genes (adj.P < 0.05): ", format(n_sig, big.mark = ","))
    }
  }
  
  return(final_result)
}

#' Escape regex metacharacters in a group label
#'
#' @param x Character vector.
#' @return `x` with every regex metacharacter backslash-escaped.
.escape_regex <- function(x) gsub("([][{}()*+?.\\^$|])", "\\\\\\1", x)

#' Substitute original group labels for their design-matrix names
#'
#' Contrast formulas are written with original group labels
#' while the design matrix uses `make.names()`
#' forms. Labels here are full of regex metacharacters ("CD4+",
#' "Tcm (CD27+)"), so each label is escaped before matching; lookarounds are
#' used instead of `\\b` (which fails for labels that begin or end with
#' punctuation); and substitution happens in one pass, longest labels first,
#' so "CD4+ Tcm" wins over "CD4+" and a replacement can never be re-matched.
#'
#' @param formula A contrast formula as a character scalar.
#' @param original_groups The original group labels appearing in formulas.
#' @param reverse_names_map Named map from original label to design name.
#' @return The formula with every label replaced by its design name.
.substitute_group_names <- function(formula, original_groups, reverse_names_map) {
  formula <- as.character(formula)
  if (length(original_groups) == 0 || !nzchar(formula)) return(formula)

  ord <- original_groups[order(nchar(original_groups), decreasing = TRUE)]
  alt <- paste(vapply(ord, .escape_regex, character(1)), collapse = "|")

  pat <- sprintf("`(?:%s)`|(?<![A-Za-z0-9._])(?:%s)(?![A-Za-z0-9._])", alt, alt)

  m <- gregexpr(pat, formula, perl = TRUE)
  hits <- regmatches(formula, m)[[1]]
  if (length(hits) == 0) return(formula)

  repl <- vapply(hits, function(h) {
    key <- sub("^`", "", sub("`$", "", h))
    val <- reverse_names_map[[key]]
    if (is.null(val) || is.na(val)) h else as.character(val)
  }, character(1), USE.NAMES = FALSE)

  regmatches(formula, m) <- list(repl)
  formula
}
