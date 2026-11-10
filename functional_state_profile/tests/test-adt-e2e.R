skip_adt_e2e <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")
  if (!exists("classify_functional_states_ucell", mode = "function"))
    skip("classify_functional_state.R not sourced")
}

write_adt_markers <- function() {
  markers <- list(
    "_comment" = "test fixture - canonical names, resolved via aliases",
    "_aliases" = list(
      "PD-1"  = list("CD279", "PDCD1"),
      "TIM3"  = list("CD366", "HAVCR2"),
      "CD62L" = list("SELL"),
      "CD127" = list("IL7R")
    ),
    naiveish = list(core = list("CD62L", "CCR7", "CD45RA", "CD127")),
    exhaustedish = list(core = list("PD-1", "TIM3", "TIGIT", "CD39"))
  )
  f <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(markers, auto_unbox = TRUE, pretty = TRUE), f)
  f
}

make_adt_seurat <- function(n_cells = 40, seed = 7) {
  set.seed(seed)
  suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

  adt_features <- c("Ab.SELL",
                    "CCR7-TotalSeqC",
                    "CD45RA.AbSeq",
                    "Hu-IL7R",
                    "CD279-TotalSeqC",
                    "Ab.CD366",
                    "TIGIT.TotalSeqC",
                    "CD39-AB",
                    "CD34-TotalSeqC",
                    "IgG1-isotype-ctrl")
  cell_names <- paste0("Cell", sprintf("%03d", seq_len(n_cells)))

  adt_counts <- matrix(rpois(length(adt_features) * n_cells, lambda = 5),
                       nrow = length(adt_features), ncol = n_cells,
                       dimnames = list(adt_features, cell_names))

  adt_counts[5:8, 1:20]   <- adt_counts[5:8, 1:20]  + rpois(4 * 20, 60)
  adt_counts[1:4, 21:40]  <- adt_counts[1:4, 21:40] + rpois(4 * 20, 60)

  rna <- matrix(rpois(5 * n_cells, 1), nrow = 5,
                dimnames = list(paste0("GENE", 1:5), cell_names))
  obj <- CreateSeuratObject(counts = as(rna, "sparseMatrix"), project = "ADT_E2E")
  obj[["ADT"]] <- CreateAssayObject(counts = as(adt_counts, "sparseMatrix"))
  obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                       margin = 2, verbose = FALSE)
  obj
}

test_that("ADT assay scores end-to-end through alias expansion + token matching", {
  skip_adt_e2e()
  obj <- make_adt_seurat()
  f   <- write_adt_markers(); on.exit(unlink(f))

  out <- classify_functional_states_ucell(
    seurat_obj       = obj,
    marker_json_path = f,
    assay            = "ADT",
    maxRank          = 10,
    use_hvg          = FALSE,
    verbose          = FALSE
  )

  md <- out@meta.data
  expect_true(all(c("naiveish", "exhaustedish",
                    "dominant_state", "dominant_state_score",
                    "dominant_state_margin") %in% colnames(md)))

  expect_identical(as.character(md$dominant_state[1:20]),  rep("exhaustedish", 20))
  expect_identical(as.character(md$dominant_state[21:40]), rep("naiveish", 20))
})

test_that("an ADT marker absent from the panel does not sink the run; a fully absent state errors", {
  skip_adt_e2e()
  obj <- make_adt_seurat()

  markers <- list(
    "_aliases" = list("PD-1" = list("CD279")),
    okstate  = list(core = list("PD-1", "TIGIT")),
    gonestate = list(core = list("CD999", "NOTONPANEL"))
  )
  f <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(markers, auto_unbox = TRUE, pretty = TRUE), f)
  on.exit(unlink(f))

  expect_error(
    classify_functional_states_ucell(obj, f, assay = "ADT",
                                     maxRank = 10, verbose = FALSE),
    "gonestate")
})

test_that("the shipped ADT marker file resolves against a realistic decorated panel", {
  skip_adt_e2e()

  shipped <- local({
    fn <- "cell_state_markers_adt.json"
    cands <- c(fn, file.path("..", fn),
               file.path(dirname(dirname(FIXTURES)), fn),
               file.path(dirname(FIXTURES), "..", fn))
    hit <- cands[file.exists(cands)][1]
    if (is.na(hit)) stop("cannot locate shipped ADT marker file")
    normalizePath(hit)
  })
  d       <- jsonlite::fromJSON(shipped, simplifyVector = TRUE)
  states  <- d[!startsWith(names(d), "_")]
  all_mk  <- unique(unlist(lapply(states, function(s) c(s$core, s$additional)),
                           use.names = FALSE))

  suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
  set.seed(11)
  feats <- vapply(seq_along(all_mk), function(i) {
    if (i %% 2L) paste0("Ab.", gsub("[^A-Za-z0-9]+", ".", all_mk[i]))
    else paste0(gsub("[^A-Za-z0-9]+", ".", all_mk[i]), "-TotalSeqC")
  }, character(1))
  n_cells <- 30
  cells <- paste0("C", seq_len(n_cells))
  adt <- matrix(rpois(length(feats) * n_cells, 20), nrow = length(feats),
                dimnames = list(feats, cells))
  rna <- matrix(rpois(5 * n_cells, 1), nrow = 5,
                dimnames = list(paste0("G", 1:5), cells))
  obj <- CreateSeuratObject(counts = as(rna, "sparseMatrix"))
  obj[["ADT"]] <- CreateAssayObject(counts = as(adt, "sparseMatrix"))
  obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                       margin = 2, verbose = FALSE)

  expect_no_error({
    out <- classify_functional_states_ucell(obj, shipped, assay = "ADT",
                                            maxRank = length(feats),
                                            verbose = FALSE)
  })
  expect_true(all(names(states) %in% colnames(out@meta.data)))
})
