skip_gate <- function() {
  if (!exists(".is_tcell_label", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
}

test_that(".is_tcell_label accepts every Monaco T cell-referencing label", {
  skip_gate()
  monaco_t <- c(

    "T cells", "CD4+ T cells", "CD8+ T cells",

    "Naive CD4 T cells", "Naive CD8 T cells",
    "Central memory CD8 T cells", "Effector memory CD8 T cells",
    "Terminal effector CD8 T cells", "Terminal effector CD4 T cells",
    "Follicular helper T cells", "T regulatory cells",
    "Th1 cells", "Th1/Th17 cells", "Th17 cells", "Th2 cells",
    "MAIT cells", "Vd2 gd T cells", "Non-Vd2 gd T cells"
  )
  res <- .is_tcell_label(monaco_t)
  expect_true(all(res), info = paste("missed:", paste(monaco_t[!res], collapse = ", ")))
})

test_that(".is_tcell_label accepts other references' T spellings", {
  skip_gate()
  other_t <- c(
    "T_cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "T cells, CD4+, naive",
    "NKT cells",
    "Tregs",
    "gamma-delta T cells"
  )
  res <- .is_tcell_label(other_t)
  expect_true(all(res), info = paste("missed:", paste(other_t[!res], collapse = ", ")))
})

test_that(".is_tcell_label rejects every non-T label", {
  skip_gate()
  non_t <- c(

    "NK cells", "B cells", "Monocytes", "Dendritic cells", "Progenitors",

    "Natural killer cells", "Naive B cells", "Non-switched memory B cells",
    "Switched memory B cells", "Exhausted B cells", "Plasmablasts",
    "Classical monocytes", "Non classical monocytes", "Intermediate monocytes",
    "Plasmacytoid dendritic cells", "Myeloid dendritic cells",
    "Low-density basophils", "Low-density neutrophils", "Progenitor cells",

    "Mast cells", "Tissue_stem_cells", "Stem cells",

    "Epithelial_cells", "Endothelial_cells", "Erythrocytes", "HSC",
    "Neutrophils", "Basophils", "Macrophages", "Fibroblasts"
  )
  res <- .is_tcell_label(non_t)
  expect_false(any(res), info = paste("false positives:", paste(non_t[res], collapse = ", ")))
})

test_that(".is_tcell_label handles NA, empty strings, factors and zero-length input", {
  skip_gate()
  expect_identical(.is_tcell_label(c(NA, "T cells", "")), c(FALSE, TRUE, FALSE))
  expect_identical(.is_tcell_label(character(0)), logical(0))
  expect_identical(.is_tcell_label(rep(NA_character_, 3)), rep(FALSE, 3))

  expect_identical(.is_tcell_label(factor(c("T cells", "B cells"))), c(TRUE, FALSE))
})

test_that(".is_tcell_label is case-insensitive", {
  skip_gate()
  expect_true(all(.is_tcell_label(c("t cells", "T CELLS", "mait cells"))))
})

.make_gate_so <- function(main_labels, fine_labels = NULL, prefix = NULL) {
  skip_if_not_installed("Seurat")
  n <- length(main_labels)
  counts <- matrix(rpois(5 * n, 2), nrow = 5,
                   dimnames = list(paste0("g", 1:5), paste0("c", seq_len(n))))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  main_col <- if (!is.null(prefix)) paste0(prefix, "_SingleR_main") else "SingleR_main"
  fine_col <- if (!is.null(prefix)) paste0(prefix, "_SingleR_fine") else "SingleR_fine"
  so@meta.data[[main_col]] <- main_labels
  if (!is.null(fine_labels)) so@meta.data[[fine_col]] <- fine_labels
  so
}

test_that(".run_stcat_annotation errors when no SingleR gate columns exist", {
  skip_gate()
  skip_if_not_installed("Seurat")
  counts <- matrix(rpois(50, 2), nrow = 5,
                   dimnames = list(paste0("g", 1:5), paste0("c", 1:10)))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  expect_error(
    .run_stcat_annotation(so, assay = "RNA", python_script = NULL, prefix = NULL,
                          fallback_column = "SingleR_fine", cleanup = TRUE,
                          verbose = FALSE),
    "T cell gating requires SingleR annotations"
  )
})

test_that("zero T cells: STCAT is skipped and Prediction filled from fallback", {
  skip_gate()
  so <- .make_gate_so(main_labels = rep("B cells", 10),
                      fine_labels = paste0("Naive B cells ", 1:10))

  expect_warning(
    res <- .run_stcat_annotation(so, assay = "RNA", python_script = NULL,
                                 prefix = NULL, fallback_column = "SingleR_fine",
                                 cleanup = TRUE, verbose = FALSE),
    "No cells passed the T cell gate"
  )
  expect_identical(res@meta.data$Prediction, res@meta.data$SingleR_fine)
  expect_identical(res@meta.data$Cluster,    res@meta.data$SingleR_fine)
  expect_true(all(is.na(res@meta.data$Uncertainty_score)))
})

test_that("zero T cells + missing fallback column: warns and leaves NA", {
  skip_gate()
  so <- .make_gate_so(main_labels = rep("Monocytes", 10))
  w <- capture_warnings(
    res <- .run_stcat_annotation(so, assay = "RNA", python_script = NULL,
                                 prefix = NULL, fallback_column = "SingleR_fine",
                                 cleanup = TRUE, verbose = FALSE)
  )
  expect_true(any(grepl("No cells passed the T cell gate", w)))
  expect_true(any(grepl("Fallback column 'SingleR_fine' not found", w)))
  expect_true(all(is.na(res@meta.data$Prediction)))
})

test_that("zero T cells with prefix: prefixed columns are used and filled", {
  skip_gate()
  so <- .make_gate_so(main_labels = rep("NK cells", 10),
                      fine_labels = rep("Natural killer cells", 10),
                      prefix = "run2")
  expect_warning(
    res <- .run_stcat_annotation(so, assay = "RNA", python_script = NULL,
                                 prefix = "run2",
                                 fallback_column = "run2_SingleR_fine",
                                 cleanup = TRUE, verbose = FALSE),
    "No cells passed the T cell gate"
  )
  expect_identical(res@meta.data$run2_Prediction, res@meta.data$run2_SingleR_fine)
  expect_false("Prediction" %in% colnames(res@meta.data))
})

test_that("gate is the UNION of main and fine labels", {
  skip_gate()

  main <- c("T cells", "CD4+ T cells", "CD8+ T cells", "NK cells", "NK cells",
            rep("B cells", 5))
  fine <- c(rep("Natural killer cells", 3), "MAIT cells", "Th17 cells",
            rep("Naive B cells", 5))
  m_main <- .is_tcell_label(main)
  m_fine <- .is_tcell_label(fine)
  expect_identical(m_main | m_fine,
                   c(rep(TRUE, 5), rep(FALSE, 5)))

  expect_lt(sum(m_main), 5)
  expect_lt(sum(m_fine), 5)
})
