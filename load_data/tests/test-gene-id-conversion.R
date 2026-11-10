test_that("gene symbols pass through unchanged", {
  mat <- matrix(1:9, nrow = 3,
                dimnames = list(c("CD3E", "CD4", "CD8A"), paste0("C", 1:3)))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_identical(rownames(result), c("CD3E", "CD4", "CD8A"))
})

test_that("non-Ensembl IDs (mixed case) pass through unchanged", {

  genes <- c("ENS_FAKE", "NOTENSEMBL", "GENE001")
  mat   <- matrix(0, nrow = 3, ncol = 2, dimnames = list(genes, c("A","B")))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_identical(rownames(result), genes)
})

test_that("matrix with NULL rownames passes through unchanged", {
  mat    <- matrix(1:4, nrow = 2)
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_null(rownames(result))
})

test_that("empty gene vector passes through unchanged", {
  mat    <- matrix(numeric(0), nrow = 0, ncol = 5,
                   dimnames = list(character(0), paste0("C", 1:5)))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_equal(nrow(result), 0L)
})

test_that("unsupported species throws an informative error", {
  genes <- paste0("ENSG", sprintf("%011d", 1:5))
  mat   <- matrix(0, nrow = 5, ncol = 2, dimnames = list(genes, c("A","B")))
  expect_error(
    .convert_gene_ids(mat, species = "zebrafish"),
    regexp = "species"
  )
})

test_that("AnnotationDbi conversion is attempted for Ensembl IDs", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  genes  <- c("ENSG00000139618",
               "ENSG00000141736")
  mat    <- matrix(0, nrow = 2, ncol = 2, dimnames = list(genes, c("A","B")))
  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)

  expect_false(any(grepl("^ENSG", rownames(result))),
               info = "Expected symbols, got Ensembl IDs")
  expect_equal(length(rownames(result)), 2L)
})

test_that("duplicate symbols after conversion are made unique", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  genes <- c("SYMBOL_A", "SYMBOL_A")
  mat   <- matrix(0, nrow = 2, ncol = 2, dimnames = list(genes, c("A","B")))

  result <- .convert_gene_ids(mat, species = "human", verbose = FALSE)
  expect_equal(rownames(result), genes)
})
