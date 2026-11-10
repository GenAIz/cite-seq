skip_vr <- function() {
  if (!exists(".rescue_missing_vocab", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
  skip_if_not_installed("Seurat")
}

.make_vr_so <- function(n_plt = 4, n_ery = 4, n_bg = 42, seed = 11) {
  set.seed(seed)
  plt_genes <- c("PPBP", "PF4", "NRGN", "GNG11", "TUBB1")
  ery_genes <- c("HBB", "HBA1", "HBA2", "ALAS2", "AHSP")
  filler    <- paste0("FILLER", 1:40)
  genes <- c(plt_genes, ery_genes, filler)
  n <- n_plt + n_ery + n_bg
  m <- matrix(rpois(length(genes) * n, 1), nrow = length(genes),
              dimnames = list(genes, paste0("c", seq_len(n))))

  m[c(plt_genes, ery_genes), ] <- rpois(10 * n, 0.05)
  m[plt_genes, seq_len(n_plt)]         <- m[plt_genes, seq_len(n_plt)] + 50
  m[ery_genes, n_plt + seq_len(n_ery)] <- m[ery_genes, n_plt + seq_len(n_ery)] + 50
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(m, "sparseMatrix")))
  so <- suppressWarnings(Seurat::NormalizeData(so, verbose = FALSE))
  so@meta.data$SingleR_main <- "Monocytes"
  so@meta.data$SingleR_fine <- "Classical monocytes"
  so
}

test_that("platelet- and erythrocyte-like cells are rescued; background is not", {
  skip_vr()
  so  <- .make_vr_so()
  res <- .rescue_missing_vocab(so, assay = "RNA", verbose = FALSE)
  md  <- res@meta.data
  expect_identical(unname(md$SingleR_main[1:4]),  rep("Platelets", 4))
  expect_identical(unname(md$SingleR_fine[1:4]),  rep("Platelets", 4))
  expect_identical(unname(md$SingleR_main[5:8]), rep("Erythrocytes", 4))
  expect_identical(unname(md$SingleR_main[9:50]), rep("Monocytes", 42))
  expect_identical(unname(md$SingleR_fine[9:50]), rep("Classical monocytes", 42))
})

test_that("audit columns record what fired and how strongly", {
  skip_vr()
  res <- .rescue_missing_vocab(.make_vr_so(), assay = "RNA", verbose = FALSE)
  md  <- res@meta.data
  expect_true(all(c("vocab_rescue", "vocab_rescue_score") %in% colnames(md)))
  expect_identical(unname(md$vocab_rescue[1:4]),  rep("Platelets", 4))
  expect_identical(unname(md$vocab_rescue[5:8]), rep("Erythrocytes", 4))
  expect_true(all(is.na(md$vocab_rescue[9:50])))
  expect_true(all(md$vocab_rescue_score[1:8] >= 2))
  expect_true(all(is.na(md$vocab_rescue_score[9:50])))
})

test_that("evidence thresholds hold: weak signal is not rescued", {
  skip_vr()
  so <- .make_vr_so()

  res <- .rescue_missing_vocab(so, assay = "RNA", min_z = 50, verbose = FALSE)
  expect_true(all(res@meta.data$SingleR_main == "Monocytes"))
  expect_true(all(is.na(res@meta.data$vocab_rescue)))
})

test_that("missing signature genes: rescue skips gracefully", {
  skip_vr()
  so <- .make_vr_so()

  sigs <- list("Platelets" = c("NOPE1", "NOPE2", "NOPE3", "NOPE4"))
  res  <- .rescue_missing_vocab(so, assay = "RNA", signatures = sigs, verbose = FALSE)
  expect_true(all(res@meta.data$SingleR_main == "Monocytes"))
  expect_false("vocab_rescue" %in% colnames(res@meta.data))
})

test_that("no SingleR columns: object returned unchanged", {
  skip_vr()
  so <- .make_vr_so()
  so@meta.data$SingleR_main <- NULL
  so@meta.data$SingleR_fine <- NULL
  res <- .rescue_missing_vocab(so, assay = "RNA", verbose = FALSE)
  expect_false("vocab_rescue" %in% colnames(res@meta.data))
})

test_that("prefix routes label and audit columns", {
  skip_vr()
  so <- .make_vr_so()
  names(so@meta.data)[names(so@meta.data) == "SingleR_main"] <- "run2_SingleR_main"
  names(so@meta.data)[names(so@meta.data) == "SingleR_fine"] <- "run2_SingleR_fine"
  res <- .rescue_missing_vocab(so, assay = "RNA", prefix = "run2", verbose = FALSE)
  md <- res@meta.data
  expect_identical(unname(md$run2_SingleR_main[1:4]), rep("Platelets", 4))
  expect_true("run2_vocab_rescue" %in% colnames(md))
  expect_false("vocab_rescue" %in% colnames(md))
})

test_that("rescued labels never enter the T cell gate", {
  skip_vr()
  if (!exists(".is_tcell_label", mode = "function")) skip("gate helper not sourced")
  expect_false(any(.is_tcell_label(c("Platelets", "Erythrocytes"))))
})
