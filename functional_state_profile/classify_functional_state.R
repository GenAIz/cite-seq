#' Match ADT marker names against assay feature names
#'
#' TWIN of `.adt_match_features()` in
#' `alternate_workflow/script/process_seurat_pipeline.R`; `make drift` compares
#' them, so edit both together. ADT panels are not universal - vendor
#' decoration differs per dataset - so protein marker files are written in
#' canonical form and resolved here. Markers are matched as whole tokens
#' after punctuation is normalised ("Ab.CD3", "Hu-CD3", "CD56_TotalSeqC" all
#' resolve to token CD3), while feature names that merely contain the marker
#' as a substring (ABCD3, NCAM1-AS1) do not match. Multi-token markers
#' (HLA-DR) match as a consecutive token run, with the glued spelling
#' (HLADR) accepted too.
#'
#' @param features Character vector of assay feature names.
#' @param wanted Character vector of canonical marker names to look for.
#' @return Unique entries of `features` that matched any wanted marker.
.adt_match_features <- function(features, wanted) {

  #' Split names into uppercase alphanumeric tokens (punctuation as separators).
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

#' Assign each cell its dominant functional state
#'
#' Pure function over the (already name-aligned) score matrix, factored out
#' so the `min_score` / margin logic is unit-testable without a Seurat
#' object. The dominant state is `which.max` per row; with `min_score` set,
#' cells whose winning score falls below it are labelled `"unassigned"`
#' instead of being given a confident-looking label. The raw winning score
#' and the margin (top score minus runner-up) are reported either way, so
#' low-confidence and near-tie calls stay auditable.
#'
#' @param raw_scores Numeric matrix of UCell scores, cells x states, with
#' state names as column names.
#' @param min_score Optional number in `[0, 1]`; winners scoring below it are
#' labelled `"unassigned"`. `NULL` disables thresholding.
#' @return A list with `state` (character), `score` (winning score) and
#' `margin` (top score minus runner-up) per cell.
.assign_dominant_states <- function(raw_scores, min_score = NULL) {
  if (!is.matrix(raw_scores) || ncol(raw_scores) < 1 || nrow(raw_scores) < 1) {
    stop("raw_scores must be a non-empty matrix")
  }
  if (!is.null(min_score) && "unassigned" %in% colnames(raw_scores)) {
    stop("min_score cannot be used with a state named 'unassigned': ",
         "the label would be ambiguous")
  }

  top_idx <- max.col(raw_scores, ties.method = "first")
  top_score <- raw_scores[cbind(seq_len(nrow(raw_scores)), top_idx)]

  if (ncol(raw_scores) >= 2L) {
    masked <- raw_scores
    masked[cbind(seq_len(nrow(masked)), top_idx)] <- -Inf
    second <- apply(masked, 1, max)
  } else {
    second <- rep(0, nrow(raw_scores))
  }

  state <- colnames(raw_scores)[top_idx]
  if (!is.null(min_score)) {
    state[top_score < min_score] <- "unassigned"
  }

  list(state = unname(state),
       score = unname(top_score),
       margin = unname(top_score - second))
}

#' Classify T Cell Functional States Using UCell Scoring
#'
#' Assigns per-cell scores for the T cell functional states defined in the
#' marker JSON using the UCell algorithm (Mann-Whitney U statistic).
#' Each cell receives a raw UCell score (0-1) per state; scores are independent
#' measurements and are NOT normalized to sum to 100% per cell. This rank-based
#' approach is robust to dataset composition and captures cells in transitional
#' or mixed functional states.
#'
#' @param seurat_obj A Seurat object with normalized expression data
#' @param marker_json_path Path to JSON file containing marker gene definitions.
#' The JSON should have state names as keys, with "core" and "additional"
#' arrays of gene symbols for each state.
#' @param use_core_only Logical. If TRUE, uses only core markers (high
#' confidence). If FALSE, uses both core and additional markers.
#' Default: FALSE
#' @param use_hvg Logical or character. Controls whether marker genes are
#' filtered to Highly Variable Genes (HVGs) before scoring:
#' \itemize{
#' \item TRUE: Restrict markers to HVGs only. HVGs are detected
#' automatically from the Seurat object (\code{VariableFeatures()}).
#' If no HVGs are found, falls back to all genes with a warning.
#' \item FALSE (default): Use all detected marker genes regardless of
#' variability status.
#' \item "auto": Use HVGs if they exist in the object, otherwise silently
#' fall back to all genes.
#' }
#' @param assay Character. Name of assay to use for expression data.
#' Default: "RNA"
#' @param maxRank Integer. Maximum rank for UCell scoring. This caps the rank
#' values to account for data sparsity. Default: 1500. Adjust based on median
#' number of expressed genes per cell (typically set to that order of magnitude).
#' @param min_score Numeric in [0, 1] or NULL. Minimum UCell score the winning
#' state must reach for a cell to receive its label. Cells whose top score
#' falls below the threshold are labelled \code{"unassigned"} instead of
#' being given a confident-looking label from an unconditional which.max.
#' \code{NULL} (default) preserves the historical behaviour: every cell is
#' labelled with its top-scoring state regardless of how weak that score is.
#' @param ncores Integer. Number of cores for parallel processing.
#' Default: 1 (single-threaded)
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return A Seurat object with the following columns added to
#' \code{@meta.data}:
#' \itemize{
#' \item One column per state (e.g., \code{exhaustion}, \code{naive-like})
#' containing the raw UCell score (0-1 range)
#' \item \code{dominant_state}: The state with the highest UCell score, or
#' \code{"unassigned"} when \code{min_score} is set and no state reaches it
#' \item \code{dominant_state_score}: The raw UCell score of the top state
#' (reported even for unassigned cells, so the near-misses are auditable)
#' \item \code{dominant_state_margin}: Top score minus runner-up score.
#' Small margins mark cells whose label was effectively a coin toss
#' between two states
#' }
#'
#' @details
#' \strong{UCell Algorithm:}
#' UCell scores are based on the Mann-Whitney U statistic, calculated
#' independently for each cell using gene expression ranks. This approach is:
#' \itemize{
#' \item \emph{Robust to dataset composition}: Results are consistent
#' regardless of dataset size or heterogeneity
#' \item \emph{Memory efficient}: Works with sparse matrices; processes data 
#' in mini-batches to minimize RAM usage
#' \item \emph{Fast}: Rank-based computation is computationally efficient
#' }
#' 
#' Raw UCell scores (0-1 range) are kept as independent measurements without
#' row-wise normalization, preserving their statistical meaning.
#'
#' \strong{Functional States Classified:}
#' The states are defined entirely by the marker JSON; nothing in this function
#' is specific to a state list. The bundled RNA file currently defines:
#' \itemize{
#' \item stress-response: Heat shock response and cellular stress
#' \item anergy: Hyporesponsive state with impaired TCR signalling
#' \item cell-death: Apoptotic and necroptotic pathways
#' \item naive-like: Stemness and memory precursor characteristics
#' \item memory: Antigen-experienced compartment (TCM/TEM/TEMRA/TRM/TSCM)
#' \item virtual-memory: Antigen-inexperienced memory-phenotype cells
#' \item adhesion: Cell migration and tissue homing
#' \item senescence: Terminal differentiation
#' \item exhaustion: Dysfunction with inhibitory receptor expression
#' \item regulating: Regulatory T cell function
#' \item cytotoxicity: Cytotoxic effector function
#' \item proliferation: Active cell cycle
#' \item quiescence: Non-proliferative resting state
#' \item IFN-response: Type I interferon signalling
#' \item activated: Acute T cell activation
#' }
#' The companion ADT file (\code{cell_state_markers_adt.json}, run with
#' \code{assay = "ADT"}) defines the surface-marker subset of these states.
#'
#' @section maxRank Parameter:
#' The maxRank parameter caps rank values to mitigate single-cell data sparsity.
#' A good rule of thumb is to set maxRank to the median number of expressed
#' genes per cell. For example:
#' \code{VlnPlot(seurat_obj, features="nFeature_RNA", pt.size=0, log=TRUE)}
#' If your dataset has median ~1000 expressed genes, set maxRank=1000.
#' This is especially important for sparse technologies (e.g., spatial
#' transcriptomics, CITE-seq, antibody tags).
#'
#' @section Marker Gene Selection:
#' Core markers represent the highest-confidence genes consistently validated
#' across recent literature (2021-2025). Additional markers include genes with
#' strong evidence but potentially broader expression. Using
#' \code{use_core_only = TRUE} provides more specific signatures but may miss
#' subtle states. Using both core and additional markers increases sensitivity.
#'
#' @section Interpretation:
#' UCell scores range from 0 to 1, representing the strength of each signature:
#' \itemize{
#' \item Score \eqn{\geq} 0.5: Strong, clear signature for that state
#' \item Score 0.3-0.5: Moderate signature expression
#' \item Score < 0.3: Weak or absent signature
#' \item Multiple states \eqn{\geq} 0.3: Co-expression of functional programs
#' }
#' Unlike normalized percentages, raw UCell scores are directly comparable
#' between cells and datasets, and represent the statistical significance of
#' each signature in isolation (not relative to other signatures in the same cell).
#'
#' @examples
#' \dontrun{
#' # Core markers only, auto-detect HVGs
#' seurat_obj <- classify_functional_states_ucell(
#' seurat_obj = my_seurat,
#' marker_json_path = "t_cell_state_markers.json",
#' use_core_only = TRUE,
#' use_hvg = "auto"
#' )
#'
#' # Use all markers with automatic maxRank estimation
#' seurat_obj <- classify_functional_states_ucell(
#' seurat_obj = my_seurat,
#' marker_json_path = "t_cell_state_markers.json",
#' use_core_only = FALSE,
#' use_hvg = FALSE
#' )
#'
#' # Custom maxRank for sparse datasets
#' seurat_obj <- classify_functional_states_ucell(
#' seurat_obj = my_seurat,
#' marker_json_path = "t_cell_state_markers.json",
#' maxRank = 800 # Adjust based on your data
#' )
#'
#' # View results
#' head(seurat_obj@meta.data[, c("dominant_state", "exhaustion", "cytotoxicity")])
#'
#' # Find highly exhausted cells (score > 0.5)
#' exhausted_cells <- colnames(seurat_obj)[seurat_obj$exhaustion > 0.5]
#'
#' # Find transitional cells (no dominant state with clear signal)
#' transitional <- colnames(seurat_obj)[seurat_obj$dominant_state_score < 0.3]
#'
#' # Visualise
#' FeaturePlot(seurat_obj, features = "exhaustion")
#' DimPlot(seurat_obj, group.by = "dominant_state")
#' }
#'
#' @references
#' Andreatta, M., Carmona, S. J. (2021) UCell: Robust and scalable 
#' single-cell gene signature scoring. Computational and Structural 
#' Biotechnology Journal 19:3796-3805.
#' 
#' Andreatta, M., Carmona, S. J. (2026) UCell and pyUCell: single-cell 
#' gene signature scoring for R and Python. Bioinformatics 42(2):btag055.
#'
#' Marker genes compiled from recent literature (2021-2025):
#' \itemize{
#' \item Exhaustion: Guegan et al., Cell Rep Med 2024
#' (doi:10.1016/j.xcrm.2024.101831)
#' \item Naive-like: Allen Institute Immune Health Atlas 2024;
#' Szabo et al., Nat Commun 2019 (doi:10.1038/s41467-019-12464-3)
#' \item Anergy: Nature Communications 2021
#' (doi:10.1038/s41467-021-24442-5); JCI Insight 2021
#' \item IFN-response: PNAS 2024 (doi:10.1073/pnas.2311722121);
#' Nat Commun 2025 (doi:10.1038/s41467-024-55459-8)
#' \item Regulatory: van Lierop et al., Nat Commun 2021
#' (doi:10.1038/s41467-021-22975-7)
#' }
#'
#' @seealso
#' \code{\link[UCell]{AddModuleScore_UCell}} for direct UCell integration,
#' \code{\link[UCell]{ScoreSignatures_UCell}} for raw UCell scoring,
#' \code{\link[Seurat]{FindVariableFeatures}} for HVG detection,
#' \code{\link[Seurat]{VariableFeatures}} to inspect current HVGs
#'
#' @export
#' @import Seurat
#' @import jsonlite
#' @import UCell
#' @import Matrix
classify_functional_states_ucell <- function(seurat_obj,
                                             marker_json_path,
                                             use_core_only = FALSE,
                                             use_hvg = FALSE,
                                             assay = "RNA",
                                             maxRank = 1500,
                                             min_score = NULL,
                                             ncores = 1,
                                             verbose = TRUE) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!file.exists(marker_json_path)) {
    stop("marker_json_path file not found: ", marker_json_path)
  }
  
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object. ",
         "Available assays: ",
         paste(names(seurat_obj@assays), collapse = ", "))
  }
  
  if (!is.logical(use_hvg) && !identical(use_hvg, "auto")) {
    stop("use_hvg must be TRUE, FALSE, or \"auto\"")
  }
  
  if (!is.numeric(maxRank) || maxRank <= 0) {
    stop("maxRank must be a positive numeric value")
  }

  if (!is.null(min_score)) {
    if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score) ||
        min_score < 0 || min_score > 1) {
      stop("min_score must be NULL or a single number in [0, 1]")
    }
  }
  
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("UCell package is required. Install with: ",
         "BiocManager::install('UCell')")
  }
  
  if (verbose) {
    cat(rep("=", 70), "\n", sep = "")
    cat("T Cell Functional State Classifier (UCell-based)\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }
  
  markers <- jsonlite::fromJSON(marker_json_path)
  
  aliases <- markers[["_aliases"]]
  markers <- markers[!startsWith(names(markers), "_")]

  if (use_core_only) {
    marker_list <- lapply(markers, function(x) x$core)
    if (verbose) cat("Marker set : core only\n")
  } else {
    marker_list <- lapply(markers, function(x) c(x$core, x$additional))
    if (verbose) cat("Marker set : core + additional\n")
  }
  
  hvgs <- tryCatch(
    VariableFeatures(seurat_obj[[assay]]),
    error = function(e) character(0)
  )

  if (is.null(hvgs)) hvgs <- character(0)
  
  n_hvgs <- length(hvgs)
  hvgs_available <- n_hvgs > 0
  
  apply_hvg_filter <- FALSE
  
  if (isTRUE(use_hvg)) {
    if (!hvgs_available) {
      warning(
        "use_hvg = TRUE but no HVGs found in assay '", assay, "'. ",
        "Run FindVariableFeatures() first, or set use_hvg = FALSE. ",
        "Falling back to all genes.",
        call. = FALSE
      )
    } else {
      apply_hvg_filter <- TRUE
    }
  } else if (identical(use_hvg, "auto")) {
    apply_hvg_filter <- hvgs_available
  }

  if (verbose) {
    if (apply_hvg_filter) {
      cat(sprintf("HVG filter : ON (%d HVGs in assay '%s')\n",
                  n_hvgs, assay))
    } else if (hvgs_available && !apply_hvg_filter) {
      cat(sprintf("HVG filter : OFF (%d HVGs available but not used)\n",
                  n_hvgs))
    } else {
      cat(sprintf("HVG filter : OFF (no HVGs detected in assay '%s')\n",
                  assay))
    }
    
    cat(sprintf("UCell settings :\n"))
    cat(sprintf(" maxRank : %d\n", maxRank))
    cat(sprintf(" parallel cores : %d\n\n", ncores))
  }
  
  available_genes <- rownames(seurat_obj[[assay]])

  if (!identical(assay, "RNA")) {

    #' A marker set plus every configured alias of its members.
    expand <- function(mk) unique(c(mk, unlist(aliases[mk], use.names = FALSE)))
    marker_list <- lapply(marker_list, function(mk)
      .adt_match_features(available_genes, expand(mk)))
  }
  
  if (apply_hvg_filter) {
    gene_universe <- intersect(available_genes, hvgs)
    if (verbose) {
      cat(sprintf("Gene universe : %d HVGs (of %d total genes)\n\n",
                  length(gene_universe), length(available_genes)))
    }
  } else {
    gene_universe <- available_genes
    if (verbose) {
      cat(sprintf("Gene universe : %d genes (all detected)\n\n",
                  length(gene_universe)))
    }
  }
  
  if (verbose) cat("Marker gene coverage:\n")
  
  for (state in names(marker_list)) {
    present_all <- sum(marker_list[[state]] %in% available_genes)
    present_eff <- sum(marker_list[[state]] %in% gene_universe)
    total <- length(marker_list[[state]])
    
    if (verbose) {
      if (apply_hvg_filter) {
        cat(sprintf(" %-20s %3d/%3d in universe (%3d/%3d total, %.0f%% HVG overlap)\n",
                    paste0(state, ":"),
                    present_eff, total,
                    present_all, total,
                    100 * present_eff / max(present_all, 1)))
      } else {
        cat(sprintf(" %-20s %3d/%3d (%.1f%%)\n",
                    paste0(state, ":"),
                    present_eff, total,
                    100 * present_eff / total))
      }
    }
    
    if (present_eff < 3) {
      warning("Only ", present_eff, " markers in gene universe for '", state,
              "' - scores may be unreliable. Consider use_hvg = FALSE.",
              call. = FALSE)
    }
  }
  
  if (verbose) cat("\nPreparing signatures for UCell scoring...\n")
  
  signatures_filtered <- lapply(marker_list, function(genes) {
    genes[genes %in% gene_universe]
  })
  
  empty_sigs <- names(signatures_filtered)[lengths(signatures_filtered) == 0]
  if (length(empty_sigs) > 0) {

    looks_ensembl <- mean(grepl("^ENS[A-Z]*[0-9]{6,}", available_genes)) > 0.5
    hint <- if (looks_ensembl && length(empty_sigs) == length(signatures_filtered)) {
      paste0("Every state is empty and the features look like Ensembl IDs (e.g. ",
             paste(utils::head(available_genes, 2), collapse = ", "),
             "). The marker file is written in gene symbols - run the gene-ID ",
             "conversion in load_data before this step.")
    } else if (apply_hvg_filter) {
      "The HVG filter is active, so the universe is restricted to variable genes. Try use_hvg = FALSE."
    } else {
      "These marker genes are absent from the assay entirely."
    }
    stop(sprintf(paste0(
      "No marker genes survive the gene universe for state(s): %s.\n %s"),
      paste(empty_sigs, collapse = ", "), hint))
  }
  
  if (verbose) cat("\nCalculating UCell scores...\n")
  
  data_mat <- tryCatch(GetAssayData(seurat_obj, assay = assay, layer = "data"),
                       error = function(e) NULL)
  if (is.null(data_mat) || length(dim(data_mat)) != 2 || any(dim(data_mat) == 0)) {
    data_mat <- tryCatch(GetAssayData(seurat_obj, assay = assay, layer = "counts"),
                         error = function(e) NULL)
    if (is.null(data_mat) || length(dim(data_mat)) != 2 || any(dim(data_mat) == 0))
      stop(sprintf("Assay '%s' has neither a 'data' nor a 'counts' layer with values.", assay))
    if (verbose) cat(" 'data' layer empty - scoring on counts (UCell ranks are unaffected)\n")
  }
  
  if (verbose) {
    cat(sprintf(" Expression matrix: %s (%.2f GB in memory)\n",
                class(data_mat)[1],
                object.size(data_mat) / 1e9))
  }
  
  tryCatch({

    if (ncores > 1) {
      bpparam <- BiocParallel::MulticoreParam(workers = ncores)
    } else {
      bpparam <- BiocParallel::SerialParam()
    }
    
    if (verbose) cat(" Starting parallel UCell scoring (mini-batch processing)...\n")
    
    if (is.matrix(data_mat) || methods::is(data_mat, "Matrix")) {
      ucell_scores <- UCell::ScoreSignatures_UCell(
        matrix = data_mat,
        features = signatures_filtered,
        maxRank = maxRank,
        BPPARAM = bpparam
      )
    } else {
      starts <- seq(1L, ncol(data_mat), by = 20000L)
      parts <- vector("list", length(starts))
      for (i in seq_along(starts)) {
        idx <- starts[i]:min(starts[i] + 20000L - 1L, ncol(data_mat))
        parts[[i]] <- UCell::ScoreSignatures_UCell(
          matrix = methods::as(data_mat[, idx, drop = FALSE], "dgCMatrix"),
          features = signatures_filtered,
          maxRank = maxRank,
          BPPARAM = bpparam
        )
        if (verbose && length(starts) > 1)
          cat(sprintf(" chunk %d/%d scored (%d cells)\n",
                      i, length(starts), length(idx)))
      }
      ucell_scores <- do.call(rbind, parts)
    }
    
    if (verbose) {
      cat(" UCell scoring completed successfully\n")
    }
    
  }, error = function(e) {
    stop("Error during UCell scoring: ", conditionMessage(e))
  })
  
  if (verbose) cat("\nProcessing UCell scores...\n")
  
  score_cols <- colnames(ucell_scores)
  raw_scores <- as.matrix(ucell_scores)
  
  colnames(raw_scores) <- gsub("_UCell$", "", colnames(raw_scores))
  
  missing_scores <- setdiff(names(marker_list), colnames(raw_scores))
  if (length(missing_scores) > 0) {
    stop(sprintf(paste0(
      "UCell returned no scores for state(s): %s.\n",
      "Expected %d signature(s), got %d column(s): %s"),
      paste(missing_scores, collapse = ", "),
      length(marker_list), ncol(raw_scores),
      paste(colnames(raw_scores), collapse = ", ")))
  }
  raw_scores <- raw_scores[, names(marker_list), drop = FALSE]
  
  if (!is.null(rownames(raw_scores)) &&
      !identical(rownames(raw_scores), colnames(seurat_obj))) {
    common <- intersect(colnames(seurat_obj), rownames(raw_scores))
    if (length(common) != ncol(seurat_obj))
      stop(sprintf("UCell returned scores for %d of %d cells.",
                   length(common), ncol(seurat_obj)))
    raw_scores <- raw_scores[colnames(seurat_obj), , drop = FALSE]
  }
  
  assignment <- .assign_dominant_states(raw_scores, min_score = min_score)
  seurat_obj@meta.data[, names(marker_list)] <- raw_scores
  seurat_obj@meta.data$dominant_state <- assignment$state
  seurat_obj@meta.data$dominant_state_score <- assignment$score
  seurat_obj@meta.data$dominant_state_margin <- assignment$margin

  if (verbose && !is.null(min_score)) {
    n_unassigned <- sum(assignment$state == "unassigned")
    cat(sprintf("\nmin_score = %.3f: %d of %d cells (%.1f%%) below threshold -> \"unassigned\"\n",
                min_score, n_unassigned, length(assignment$state),
                100 * n_unassigned / length(assignment$state)))
  }
  
  score_ranges <- list(
    min_vals = apply(raw_scores, 2, min),
    max_vals = apply(raw_scores, 2, max),
    mean_vals = apply(raw_scores, 2, mean)
  )
  
  if (any(raw_scores < 0) || any(raw_scores > 1)) {
    warning("Some UCell scores outside expected 0-1 range - check algorithm execution.",
            call. = FALSE)
  }
  
  if (verbose) {
    cat("\nState distribution (by dominant state assignment):\n")
    state_counts <- sort(table(seurat_obj@meta.data$dominant_state),
                         decreasing = TRUE)
    state_pcts <- 100 * state_counts / sum(state_counts)
    
    for (state in names(state_counts)) {
      cat(sprintf(" %-20s %6d cells (%.1f%%)\n",
                  paste0(state, ":"), state_counts[state], state_pcts[state]))
    }
    
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Classification complete!\n")
    cat(rep("=", 70), "\n\n", sep = "")
    
    cat("Added to metadata:\n")
    cat(" - State scores :", paste(names(marker_list), collapse = ", "), "\n")
    cat(" (raw UCell scores, 0-1 range)\n")
    cat(" - dominant_state : State with highest UCell score",
        if (!is.null(min_score)) sprintf(" (\"unassigned\" below %.3f)", min_score) else "",
        "\n", sep = "")
    cat(" - dominant_state_score: Raw UCell score of dominant state\n")
    cat(" - dominant_state_margin: Top score minus runner-up score\n\n")
    
    cat("Score ranges by state:\n")
    for (i in seq_along(score_ranges$min_vals)) {
      state <- names(score_ranges$min_vals)[i]
      cat(sprintf(" %-20s Min: %.4f | Mean: %.4f | Max: %.4f\n",
                  paste0(state, ":"),
                  score_ranges$min_vals[i],
                  score_ranges$mean_vals[i],
                  score_ranges$max_vals[i]))
    }
    
    cat(sprintf("\nAlgorithm details:\n"))
    cat(sprintf(" - Method : UCell (Mann-Whitney U statistic)\n"))
    cat(sprintf(" - Scoring : Rank-based, 0-1 scale (no normalization)\n"))
    cat(sprintf(" - HVG filter used : %s\n\n", apply_hvg_filter))
  }
  
  return(seurat_obj)
}

#' Validate Functional State Classification
#'
#' Checks that UCell scores are in the expected 0-1 range and reports basic statistics about the classification results.
#'
#' @param seurat_obj A Seurat object with functional state classifications
#' produced by \code{\link{classify_functional_states_ucell}}.
#' @param states Character vector of state names to check. If \code{NULL},
#' state names are inferred from \code{dominant_state} values in the
#' metadata.
#'
#' @return Invisibly returns a logical: \code{TRUE} if all scores are in 
#' the expected 0-1 range, \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' validate_classification(seurat_obj)
#' }
#'
#' @export
validate_classification <- function(seurat_obj, states = NULL) {
  
  if (is.null(states)) {
    if (!"dominant_state" %in% colnames(seurat_obj@meta.data)) {
      stop("Cannot detect state columns. ",
           "Please specify the 'states' parameter.")
    }

    states <- setdiff(unique(seurat_obj@meta.data$dominant_state), "unassigned")
    if (length(states) == 0) {
      stop("All cells are 'unassigned'; cannot infer state columns. ",
           "Please specify the 'states' parameter.")
    }
  }
  
  missing_cols <- states[!states %in% colnames(seurat_obj@meta.data)]
  if (length(missing_cols) > 0) {
    stop("State columns not found: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("=== Classification Validation ===\n\n")
  
  score_mat <- as.matrix(seurat_obj@meta.data[, states, drop = FALSE])
  
  cat("UCell score ranges (expected: 0-1):\n")
  cat(sprintf(" Min : %.6f\n", min(score_mat)))
  cat(sprintf(" Max : %.6f\n", max(score_mat)))
  cat(sprintf(" Mean : %.6f\n", mean(score_mat)))
  cat(sprintf(" Median : %.6f\n", median(score_mat)))
  cat(sprintf(" SD : %.6f\n\n", sd(score_mat)))
  
  all_valid <- all(score_mat >= 0 & score_mat <= 1)
  
  if (all_valid) {
    cat("PASS: All scores in expected 0-1 range\n\n")
  } else {
    bad <- sum(score_mat < 0 | score_mat > 1)
    cat(sprintf("WARNING: %d/%d score values outside 0-1 range\n\n",
                bad, length(score_mat)))
  }
  
  cat("Dominant state distribution:\n")
  state_table <- sort(table(seurat_obj@meta.data$dominant_state),
                      decreasing = TRUE)
  state_pct <- 100 * state_table / sum(state_table)
  
  for (state in names(state_table)) {
    cat(sprintf(" %-20s %6d (%.1f%%)\n",
                state, state_table[state], state_pct[state]))
  }
  
  invisible(all_valid)
}
