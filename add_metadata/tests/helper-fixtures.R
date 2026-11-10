`%||%` <- function(a, b) if (!is.null(a)) a else b

FIXTURES <- file.path(dirname(sys.frame(1)$ofile %||% "."), "fixtures")

make_mini_seurat <- function(n = 10) {
  suppressPackageStartupMessages(library(Seurat))
  counts <- matrix(sample(0:20, 20 * n, replace = TRUE),
                   nrow = 20, ncol = n,
                   dimnames = list(paste0("G", 1:20), paste0("Cell", 1:n)))
  obj <- CreateSeuratObject(counts = as(counts, "sparseMatrix"))
  obj$orig.ident <- paste0("Sample", rep(1:5, each = 2)[1:n])
  obj
}

write_meta_csv <- function(path, samples = paste0("Sample", 1:5),
                            add_col = TRUE) {
  df <- data.frame(
    Sample    = samples,
    Treatment = rep(c("Control", "Treated"), length.out = length(samples)),
    stringsAsFactors = FALSE
  )
  if (add_col) df$Score <- seq_along(samples) * 1.5
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}
