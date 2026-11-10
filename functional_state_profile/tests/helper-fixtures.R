`%||%` <- function(a, b) if (!is.null(a)) a else b

FIXTURES <- file.path(dirname(sys.frame(1)$ofile %||% "."), "fixtures")
TEST_MARKERS_JSON <- file.path(FIXTURES, "test_markers.json")

make_fsp_seurat <- function(n_cells = 50, seed = 42) {
  set.seed(seed)
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
  })

  marker_genes <- c("PDCD1", "LAG3", "HAVCR2", "CTLA4", "TIGIT",
                     "CCR7", "LEF1", "TCF7", "SELL", "IL7R",
                     "GZMB", "PRF1", "NKG7", "GNLY", "GZMA")
  filler_genes <- paste0("FILLER", sprintf("%04d", seq_len(85)))
  gene_names   <- c(marker_genes, filler_genes)
  cell_names   <- paste0("Cell", sprintf("%04d", seq_len(n_cells)))

  rna_counts <- matrix(rpois(length(gene_names) * n_cells, lambda = 2),
                       nrow = length(gene_names), ncol = n_cells,
                       dimnames = list(gene_names, cell_names))

  rna_counts[1:3, 1:15] <- rna_counts[1:3, 1:15] + rpois(3 * 15, 8)

  rna_counts[11:13, 16:30] <- rna_counts[11:13, 16:30] + rpois(3 * 15, 8)

  rna_counts[6:8, 31:45] <- rna_counts[6:8, 31:45] + rpois(3 * 15, 8)

  obj <- CreateSeuratObject(counts = as(rna_counts, "sparseMatrix"),
                             project = "FSP_Test")
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
  })

  obj
}
