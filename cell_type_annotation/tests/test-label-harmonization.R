skip_lh <- function() {
  if (!exists(".harmonize_prediction_labels", mode = "function"))
    skip("cell_type_annotation.R not sourced (Seurat unavailable at runner load)")
  skip_if_not_installed("Seurat")
}

.make_lh_so <- function(preds, prefix = NULL) {
  n <- length(preds)
  m <- matrix(rpois(5 * n, 2), nrow = 5,
              dimnames = list(paste0("g", 1:5), paste0("c", seq_len(n))))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = as(m, "sparseMatrix")))
  p <- function(x) if (!is.null(prefix)) paste0(prefix, "_", x) else x
  so@meta.data[[p("Prediction")]]   <- preds
  so@meta.data[[p("SingleR_fine")]] <- preds
  so
}

test_that("all eleven concept twins are renamed to STCAT vocabulary", {
  skip_lh()
  map_in  <- c("MAIT cells", "Naive CD4 T cells", "Naive CD8 T cells",
               "Central memory CD8 T cells", "Effector memory CD8 T cells",
               "Terminal effector CD8 T cells", "Terminal effector CD4 T cells",
               "T regulatory cells", "Th1 cells", "Th17 cells",
               "Vd2 gd T cells", "Non-Vd2 gd T cells")
  map_out <- c("MAIT", "CD4 Tn", "CD8 Tn", "CD8 Tcm", "CD8 Tem", "CD8 Temra",
               "CD4 Temra", "CD4 Treg", "CD4 Th1", "CD4 Th17", "Tgd", "Tgd")
  so  <- .make_lh_so(map_in)
  res <- .harmonize_prediction_labels(so, verbose = FALSE)
  expect_identical(unname(res@meta.data$Prediction), map_out)
})

test_that("no-twin Monaco labels and STCAT/other labels are untouched", {
  skip_lh()
  keep <- c("Th2 cells", "Follicular helper T cells", "Th1/Th17 cells",
            "CD4 Tn", "CD8 Temra", "MAIT", "Tgd", "Heterogeneous",
            "Classical monocytes", "Naive B cells", "Natural killer cells",
            "Platelets", "Erythrocytes", "CD4 Trm cell-death", NA)
  so  <- .make_lh_so(keep)
  res <- .harmonize_prediction_labels(so, verbose = FALSE)
  expect_identical(res@meta.data$Prediction, so@meta.data$Prediction)
})

test_that("SingleR_fine is never modified", {
  skip_lh()
  so  <- .make_lh_so(c("MAIT cells", "Vd2 gd T cells", "Th2 cells"))
  res <- .harmonize_prediction_labels(so, verbose = FALSE)
  expect_identical(res@meta.data$SingleR_fine,
                   c("MAIT cells", "Vd2 gd T cells", "Th2 cells"))
  expect_identical(unname(res@meta.data$Prediction), c("MAIT", "Tgd", "Th2 cells"))
})

test_that("a factor Prediction column is harmonized by LABEL, not by level code", {
  skip_lh()

  preds <- c("MAIT cells", "Th2 cells", "Vd2 gd T cells", "MAIT cells")
  so <- .make_lh_so(preds)
  so@meta.data$Prediction <- factor(preds)
  res <- .harmonize_prediction_labels(so, verbose = FALSE)
  out <- res@meta.data$Prediction
  expect_identical(as.character(out), c("MAIT", "Th2 cells", "Tgd", "MAIT"))

  expect_false(any(grepl("^[0-9]+$", as.character(out))))
})

test_that("prefix routes the Prediction column", {
  skip_lh()
  so  <- .make_lh_so(c("MAIT cells", "T regulatory cells"), prefix = "run2")
  res <- .harmonize_prediction_labels(so, prefix = "run2", verbose = FALSE)
  expect_identical(unname(res@meta.data$run2_Prediction), c("MAIT", "CD4 Treg"))
  expect_false("Prediction" %in% colnames(res@meta.data))
})

test_that("object without a Prediction column passes through", {
  skip_lh()
  so <- .make_lh_so(c("MAIT cells", "MAIT cells"))
  so@meta.data$Prediction <- NULL
  expect_identical(colnames(.harmonize_prediction_labels(so, verbose = FALSE)@meta.data),
                   colnames(so@meta.data))
})

test_that("the map targets only STCAT-vocabulary names (no identity or chained renames)", {
  skip_lh()
  map <- .prediction_harmonization_map()

  expect_length(intersect(names(map), unname(map)), 0)

  expect_identical(sort(names(map)[map == "Tgd"]),
                   c("Non-Vd2 gd T cells", "Vd2 gd T cells"))
})

test_that("None-T fallback sets Uncertainty_score to NA (provenance rule)", {
  skip_lh()

  cand <- c("cell_type_annotation.R", "../cell_type_annotation.R",
            "script/cell_type_annotation.R", "../script/cell_type_annotation.R")
  path <- cand[file.exists(cand)][1]
  skip_if(is.na(path), "cell_type_annotation.R not found from test wd")
  block <- paste(readLines(path), collapse = "\n")
  expect_true(grepl(
    "col_uncertainty\\]\\]\\[none_t_cells\\] <- NA", block),
    info = "None-T fallback must NA the uncertainty score")
})
