`%||%` <- function(a, b) if (!is.null(a)) a else b

make_annotatable_seurat <- function(n_cells = 50, n_genes = 200, seed = 42) {
  set.seed(seed)
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
  })

  gene_names <- c(

    "CD3E", "CD4", "CD8A", "CD8B", "CD14", "CD19", "MS4A1", "NCAM1",
    "FOXP3", "IL2RA", "PTPRC", "LYZ", "NKG7", "GZMB", "PRF1",
    "IL7R", "CCR7", "LEF1", "SELL", "TCF7",
    paste0("GENE", sprintf("%04d", seq_len(n_genes - 20)))
  )
  cell_names <- paste0("Cell", sprintf("%04d", seq_len(n_cells)))

  rna_counts <- matrix(0L, nrow = length(gene_names), ncol = n_cells,
                       dimnames = list(gene_names, cell_names))
  for (i in seq_along(gene_names)) {
    rna_counts[i, ] <- sample(0:8, n_cells, replace = TRUE)
  }

  obj <- CreateSeuratObject(counts = as(rna_counts, "sparseMatrix"),
                             project = "TestAnnotation")
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 50, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, npcs = 10, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:5, verbose = FALSE)

    obj <- FindClusters(obj, resolution = 0.5, cluster.name = "RNA_clusters",
                        verbose = FALSE)
  })

  obj
}
