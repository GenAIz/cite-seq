`%||%` <- function(a, b) if (!is.null(a)) a else b

make_dge_seurat <- function(n_cells = 120, n_genes = 100, n_samples = 6,
                             n_conditions = 2, n_timepoints = 2,
                             n_clusters = 3, seed = 42,
                             cluster_signal = FALSE,
                             condition_names = NULL) {
  set.seed(seed)
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
  })

  gene_names <- paste0("GENE", sprintf("%04d", seq_len(n_genes)))
  cell_names <- paste0("Cell", sprintf("%04d", seq_len(n_cells)))

  rna_counts <- matrix(rpois(n_genes * n_cells, lambda = 5),
                       nrow = n_genes, ncol = n_cells,
                       dimnames = list(gene_names, cell_names))

  conditions <- condition_names %||% rep(c("Control", "Treated"), length.out = n_conditions)
  n_conditions <- length(conditions)
  timepoints <- paste0("T", seq_len(n_timepoints))
  samples    <- paste0("S", sprintf("%02d", seq_len(n_samples)))
  clusters   <- paste0("Cluster", seq_len(n_clusters))

  condition_assign <- rep(conditions, each = n_cells / n_conditions)
  timepoint_assign <- rep(rep(timepoints, each = n_cells / (n_conditions * n_timepoints)),
                           times = n_conditions)
  sample_assign    <- rep(samples, each = n_cells / n_samples)
  cluster_assign   <- rep(clusters, length.out = n_cells)

  treated_idx <- which(condition_assign == "Treated")
  rna_counts[1:20, treated_idx] <- rna_counts[1:20, treated_idx] + rpois(20 * length(treated_idx), 10)

  if (cluster_signal) {
    c1 <- which(cluster_assign[seq_len(n_cells)] == "Cluster1")
    rna_counts[21:40, c1] <- rna_counts[21:40, c1] +
      rpois(20 * length(c1), lambda = 25)
  }

  obj <- CreateSeuratObject(counts = as(rna_counts, "sparseMatrix"),
                             project = "DGE_Test")
  obj$Condition     <- condition_assign[1:ncol(obj)]
  obj$Timepoint     <- timepoint_assign[1:ncol(obj)]
  obj$SampleID      <- sample_assign[1:ncol(obj)]
  obj$CellType      <- cluster_assign[1:ncol(obj)]
  obj$ParticipantID <- paste0("P", as.integer(factor(sample_assign[1:ncol(obj)])))

  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
  })

  obj
}

make_mini_dge_seurat <- function(seed = 42) {
  make_dge_seurat(n_cells = 60, n_genes = 50, n_samples = 6,
                   n_conditions = 2, n_timepoints = 2,
                   n_clusters = 2, seed = seed)
}

make_mini_dge_seurat_clustered <- function(seed = 42) {
  make_dge_seurat(n_cells = 60, n_genes = 50, n_samples = 6,
                   n_conditions = 2, n_timepoints = 2,
                   n_clusters = 2, seed = seed, cluster_signal = TRUE)
}
