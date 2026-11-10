skip_proc <- function() {
  skip_if_not_installed("Seurat")
  if (!exists(".clean_adt_features", mode = "function"))
    skip("process_seurat_pipeline.R not sourced (Seurat unavailable at runner load)")
}

test_that(".clean_adt_features drops isotype controls and non-marker tags", {
  skip_proc()
  obj    <- make_cite_seurat()
  before <- nrow(obj[["ADT"]])
  cleaned <- .clean_adt_features(obj, drop_pattern = "isotype|^Set1-Sample",
                                 verbose = FALSE)
  feats <- rownames(cleaned[["ADT"]])
  expect_false(any(grepl("isotype", feats, ignore.case = TRUE)))
  expect_false(any(grepl("^Set1-Sample", feats)))
  expect_true(all(c("Ab.CD3", "Ab.CD4", "Ab.CD8") %in% feats))
  expect_lt(length(feats), before)
})

test_that(".clean_adt_features keeps the assay unchanged if < 2 features would remain", {
  skip_proc()
  obj <- make_cite_seurat()
  cleaned <- .clean_adt_features(obj, drop_pattern = ".", verbose = FALSE)
  expect_equal(nrow(cleaned[["ADT"]]), nrow(obj[["ADT"]]))
})

test_that(".clean_adt_features honours keep_pattern", {
  skip_proc()
  obj <- make_cite_seurat()
  cleaned <- .clean_adt_features(obj, keep_pattern = "^Ab\\.",
                                 drop_pattern = "isotype", verbose = FALSE)
  feats <- rownames(cleaned[["ADT"]])
  expect_true(all(grepl("^Ab\\.", feats)))
  expect_false(any(grepl("^Set1-Sample", feats)))
  expect_false(any(grepl("isotype", feats, ignore.case = TRUE)))
})

test_that(".clean_adt_features is a no-op without an ADT assay", {
  skip_proc()
  obj <- make_cite_seurat()
  obj[["ADT"]] <- NULL
  expect_identical(.clean_adt_features(obj, verbose = FALSE), obj)
})

test_that(".filter_genes_by_min_cells removes rarely-detected genes", {
  skip_proc()
  obj    <- make_cite_seurat()
  before <- nrow(obj)
  filt   <- .filter_genes_by_min_cells(obj, min_cells = 3, verbose = FALSE)
  expect_lt(nrow(filt), before)
  expect_true("Gene1" %in% rownames(filt))
})

test_that(".filter_genes_by_min_cells with min_cells <= 0 is a no-op", {
  skip_proc()
  obj <- make_cite_seurat()
  expect_equal(nrow(.filter_genes_by_min_cells(obj, min_cells = 0, verbose = FALSE)),
               nrow(obj))
})

test_that(".drop_small_batches removes undersized batches", {
  skip_proc()
  obj     <- make_cite_seurat()
  dropped <- .drop_small_batches(obj, batch_column = "orig.ident",
                                 min_cells_per_sample = 10, verbose = FALSE)
  expect_false("BatchTiny" %in% dropped$orig.ident)
  expect_true("BatchA" %in% dropped$orig.ident)
  expect_lt(ncol(dropped), ncol(obj))
})

test_that(".drop_small_batches never drops every cell", {
  skip_proc()
  obj  <- make_cite_seurat()
  kept <- .drop_small_batches(obj, batch_column = "orig.ident",
                              min_cells_per_sample = 1e6, verbose = FALSE)
  expect_equal(ncol(kept), ncol(obj))
})

test_that(".drop_small_batches no-ops on a missing batch column", {
  skip_proc()
  obj <- make_cite_seurat()
  out <- .drop_small_batches(obj, batch_column = "NoSuchCol",
                             min_cells_per_sample = 10, verbose = FALSE)
  expect_equal(ncol(out), ncol(obj))
})

test_that(".to_disk_backed moves RNA counts on-disk and leaves ADT in memory", {
  skip_proc()
  skip_if_not_installed("BPCells")
  obj <- make_cite_seurat()
  out <- .to_disk_backed(obj, cache_dir = tempfile("bpc_"),
                         assays = "RNA", verbose = FALSE)
  expect_s4_class(out, "Seurat")
  rna_counts <- SeuratObject::LayerData(out, assay = "RNA", layer = "counts")
  adt_counts <- SeuratObject::LayerData(out, assay = "ADT", layer = "counts")
  expect_true(inherits(rna_counts, "IterableMatrix"))
  expect_false(inherits(adt_counts, "IterableMatrix"))
  expect_equal(ncol(out), ncol(obj))
})

test_that(".detect_rna_batch returns the expected result structure", {
  skip_proc()
  obj <- make_cite_seurat()
  suppressWarnings({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 10, verbose = FALSE)
  })
  res <- .detect_rna_batch(obj, rna_reduction = "pca", batch_column = "orig.ident",
                           n_batches = 2, lisi_threshold = 0.7,
                           silhouette_threshold = 0.15, silhouette_n_cells = 50,
                           verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("median_lisi", "mean_silhouette", "lisi_flag",
                    "silhouette_flag", "correction_recommended") %in% names(res)))
  expect_type(res$correction_recommended, "logical")
})

test_that(".adt_degenerate_markers flags globally dead / constant markers", {
  if (!exists(".adt_degenerate_markers", mode = "function"))
    skip("process_seurat_pipeline.R not sourced")
  set.seed(1)
  n <- 200
  m <- cbind(
    good     = rnorm(n, 5, 2),
    allzero  = rep(0, n),
    constant = rep(3.5, n)
  )
  bad <- .adt_degenerate_markers(m)
  expect_true(all(c("allzero", "constant") %in% bad))
  expect_false("good" %in% bad)
})

test_that(".adt_degenerate_markers flags within-batch constant / low-complexity", {
  if (!exists(".adt_degenerate_markers", mode = "function"))
    skip("process_seurat_pipeline.R not sourced")
  set.seed(2)
  n    <- 120
  samp <- rep(c("A", "B"), each = n / 2)
  m <- cbind(
    global_ok = rnorm(n, 5, 2),
    flat_in_B = c(rnorm(n / 2, 5, 2), rep(4, n / 2)),
    lowdist   = sample(0:2, n, replace = TRUE)
  )
  bad <- .adt_degenerate_markers(m, samples = samp, min_distinct = 10L)
  expect_true(all(c("flat_in_B", "lowdist") %in% bad))
  expect_false("global_ok" %in% bad)
})

test_that(".normalize_adt applies CLR to integer ADT counts", {
  skip_proc()
  if (!exists(".normalize_adt", mode = "function")) skip("not sourced")
  obj  <- make_cite_seurat()
  out  <- .normalize_adt(obj, verbose = FALSE)
  cnts <- as.matrix(GetAssayData(out, assay = "ADT", layer = "counts"))
  dat  <- as.matrix(GetAssayData(out, assay = "ADT", layer = "data"))
  expect_false(isTRUE(all.equal(cnts, dat)))
})

test_that(".normalize_adt skips CLR when ADT counts are non-integer", {
  skip_proc()
  if (!exists(".normalize_adt", mode = "function")) skip("not sourced")
  obj <- make_cite_seurat()
  cn  <- GetAssayData(obj, assay = "ADT", layer = "counts")
  cn@x <- cn@x + 0.5
  obj <- SetAssayData(obj, assay = "ADT", layer = "counts", new.data = cn)
  out  <- .normalize_adt(obj, verbose = FALSE)
  cnts <- as.matrix(GetAssayData(out, assay = "ADT", layer = "counts"))
  dat  <- as.matrix(GetAssayData(out, assay = "ADT", layer = "data"))
  expect_equal(cnts, dat)
})

test_that(".clean_adt_features drops zero-variance markers (not just by pattern)", {
  skip_proc()
  obj <- make_cite_seurat()
  cn  <- GetAssayData(obj, assay = "ADT", layer = "counts")
  cn["Ab.CD3", ] <- 0
  obj <- SetAssayData(obj, assay = "ADT", layer = "counts", new.data = cn)

  cleaned <- .clean_adt_features(obj, drop_pattern = "isotype", verbose = FALSE)
  expect_false("Ab.CD3" %in% rownames(cleaned[["ADT"]]))
  expect_true("Ab.CD4" %in% rownames(cleaned[["ADT"]]))
})

test_that(".adt_degenerate_markers ignores batches too small to assess", {
  skip_proc()
  set.seed(3)

  mat <- matrix(rpois(43 * 4, 6), nrow = 43,
                dimnames = list(NULL, paste0("M", 1:4)))
  samples <- c(rep("Big", 40), rep("Tiny", 3))
  bad <- .adt_degenerate_markers(mat, samples = samples)
  expect_false(setequal(bad, colnames(mat)))
})

test_that(".detect_adt_batch returns the expected result structure", {
  skip_proc()
  if (!exists(".detect_adt_batch", mode = "function")) skip("not sourced")
  obj <- make_cite_seurat(n_cells = 80)
  obj$Pool <- rep(c("P1", "P2"), each = 40)
  obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                       margin = 2, verbose = FALSE)
  res <- .detect_adt_batch(obj, batch_column = "Pool", n_batches = 2,
                           wasserstein_threshold = 0.3, verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("flagged_proteins", "pct_proteins_flagged", "mean_wasserstein",
                    "correction_recommended", "per_protein_distance") %in% names(res)))
  expect_type(res$correction_recommended, "logical")
})

test_that(".adt_match_features resolves vendor-decorated marker names", {
  skip_proc()
  feats <- c("Ab.CD3", "Hu-CD19", "CD56_TotalSeqC", "CD14", "MOPC.21")
  expect_equal(.adt_match_features(feats, "CD3"),  "Ab.CD3")
  expect_equal(.adt_match_features(feats, "CD19"), "Hu-CD19")
  expect_equal(.adt_match_features(feats, "CD56"), "CD56_TotalSeqC")
  expect_equal(.adt_match_features(feats, "CD14"), "CD14")
})

test_that(".adt_match_features does not confuse digit-extended markers", {
  skip_proc()
  feats <- c("Ab.CD34", "Ab.CD38", "Ab.CD3", "CD141", "CD14", "CD200", "CD20")
  expect_equal(.adt_match_features(feats, "CD3"),  "Ab.CD3")
  expect_equal(.adt_match_features(feats, "CD14"), "CD14")
  expect_equal(.adt_match_features(feats, "CD20"), "CD20")
})

test_that(".adt_match_features returns empty for absent markers", {
  skip_proc()
  expect_length(.adt_match_features(c("Ab.CD3", "Ab.CD4"), "CD19"), 0)
})

test_that(".adt_crosslineage_doublets returns a named logical over all cells", {
  skip_proc()
  obj <- make_cite_seurat()
  res <- .adt_crosslineage_doublets(obj, assay = "ADT", verbose = FALSE)
  expect_type(res, "logical")
  expect_length(res, ncol(obj))
  expect_equal(names(res), colnames(obj))
})

test_that(".adt_crosslineage_doublets flags a constructed cross-lineage cell", {
  skip_proc()
  obj <- make_cite_seurat()
  cts <- as.matrix(SeuratObject::LayerData(obj, assay = "ADT", layer = "counts"))

  cts["Ab.CD3",  1] <- 5000
  cts["Ab.CD19", 1] <- 5000
  obj[["ADT"]] <- SeuratObject::CreateAssayObject(counts = as(cts, "sparseMatrix"))
  res <- .adt_crosslineage_doublets(obj, assay = "ADT", verbose = FALSE)
  expect_true(res[[1]])
})

test_that(".adt_crosslineage_doublets is a no-op when the ADT assay is absent", {
  skip_proc()
  obj <- make_cite_seurat()
  obj[["ADT"]] <- NULL
  res <- .adt_crosslineage_doublets(obj, assay = "ADT", verbose = FALSE)
  expect_false(any(res))
  expect_length(res, ncol(obj))
})

test_that(".detect_doublets adds the expected metadata columns", {
  skip_proc()
  obj <- make_cite_seurat()
  out <- .detect_doublets(obj, sample_column = "orig.ident",
                          min_cells_per_sample = 10L, verbose = FALSE)
  for (col in c("doublet_score", "doublet_class_rna", "adt_crosslineage",
                "doublet_class", "doublet_method"))
    expect_true(col %in% colnames(out@meta.data), info = col)
  expect_true(all(out$doublet_class %in% c("doublet", "singlet")))
  expect_false(is.null(out@misc$doublet_qc))
})

test_that(".detect_doublets falls back to orig.ident for a bad sample column", {
  skip_proc()
  obj <- make_cite_seurat()
  out <- .detect_doublets(obj, sample_column = "does_not_exist",
                          min_cells_per_sample = 10L, verbose = FALSE)
  expect_equal(out@misc$doublet_qc$sample_column, "orig.ident")
})

test_that(".remove_doublets drops flagged cells and keeps the summary", {
  skip_proc()
  obj <- make_cite_seurat()
  obj$doublet_class <- c("doublet", "doublet", rep("singlet", ncol(obj) - 2))
  obj@misc$doublet_qc <- list(n_total_doublet = 2)
  out <- .remove_doublets(obj, verbose = FALSE)
  expect_equal(ncol(out), ncol(obj) - 2)
  expect_true(all(out$doublet_class == "singlet"))
  expect_equal(out@misc$doublet_qc$n_total_doublet, 2)
})

test_that(".remove_doublets refuses to empty the object", {
  skip_proc()
  obj <- make_cite_seurat()
  obj$doublet_class <- rep("doublet", ncol(obj))
  out <- .remove_doublets(obj, verbose = FALSE)
  expect_equal(ncol(out), ncol(obj))
})

test_that(".remove_doublets is a no-op without the doublet_class column", {
  skip_proc()
  obj <- make_cite_seurat()
  expect_equal(ncol(.remove_doublets(obj, verbose = FALSE)), ncol(obj))
})

test_that(".adt_match_features does not match genes that merely contain a marker", {
  skip_proc()
  genes <- c("ABCD3", "HACD3", "NUDCD3", "PLCD3", "PTCD3")
  expect_length(.adt_match_features(genes, "CD3"), 0)
})

test_that(".adt_match_features resolves separator-decorated antibody names", {
  skip_proc()
  panel <- c("CD3", "Ab.CD3", "Hu-CD3", "CD3-TotalSeqC", "anti-human CD3",
             "Ab.CD34", "Ab.CD38")
  expect_setequal(.adt_match_features(panel, "CD3"),
                  c("CD3", "Ab.CD3", "Hu-CD3", "CD3-TotalSeqC", "anti-human CD3"))
})

test_that("ADT lineage groups contain no gene-symbol aliases", {
  skip_proc()
  toks <- unlist(.adt_lineage_groups(), use.names = FALSE)

  expect_true(all(grepl("^CD[0-9]+[A-Z]?$", toks)))
  expect_false(any(c("NCAM1", "MS4A1", "CD3E", "CD3D", "CD3G") %in% toks))
})

test_that("lineage groups still cover the mutually exclusive immune lineages", {
  skip_proc()
  g <- .adt_lineage_groups()
  expect_setequal(names(g),
                  c("T_cell", "B_cell", "NK_cell", "Myeloid", "Platelet"))
})

test_that(".detect_nonmarker_tags matches common multiplexing conventions", {
  skip_proc()
  tags <- c("HTO", "HTO1", "HTO_2", "HTO-3", "Hashtag4", "hash-tag5",
            "HT7", "BC2", "TotalSeq-C0251 Hashtag 1", "Set1-Sample1",
            "Sample2", "CMO301", "LMO3", "MULTI-A1")
  expect_setequal(.detect_nonmarker_tags(tags), tags)
})

test_that(".detect_nonmarker_tags spares genuine markers from real panels", {
  skip_proc()

  markers <- c("Ab.CD3-UCHT1", "Ab.CD14-M5E2", "Ab.CD141", "Ab.CD303",
               "Ab.CD1c", "Ab.CD226-11A8", "Ab.Isotype-RTK2071", "Ab.HLA.DR",
               "Ab.integrin.b7", "Ab.TCR.Va7.2", "AB_CD56_NCAM",
               "AB_CD274_B7_H1_PD_L1", "AB_humanCD11b", "CD19", "CD34")
  expect_length(.detect_nonmarker_tags(markers), 0)
})

test_that(".clean_adt_features removes multiplexing tags with no config at all", {
  skip_proc()
  obj <- make_cite_seurat()
  cleaned <- .clean_adt_features(obj, verbose = FALSE)
  feats <- rownames(cleaned[["ADT"]])
  expect_false(any(grepl("^Set1-Sample", feats)))
  expect_true(all(c("Ab.CD3", "Ab.CD4", "Ab.CD8") %in% feats))
})

test_that("multiplexing-tag removal is additive to a user drop_pattern", {
  skip_proc()
  obj <- make_cite_seurat()
  cleaned <- .clean_adt_features(obj, drop_pattern = "isotype", verbose = FALSE)
  feats <- rownames(cleaned[["ADT"]])
  expect_false(any(grepl("^Set1-Sample", feats)))
  expect_false(any(grepl("isotype", feats, ignore.case = TRUE)))
})

test_that("multiplexing-tag removal is skipped if it would gut the panel", {
  skip_proc()
  obj <- make_cite_seurat()
  cts <- as.matrix(SeuratObject::LayerData(obj, assay = "ADT", layer = "counts"))
  rownames(cts) <- paste0("HTO_", seq_len(nrow(cts)))
  obj[["ADT"]] <- SeuratObject::CreateAssayObject(counts = as(cts, "sparseMatrix"))
  cleaned <- .clean_adt_features(obj, verbose = FALSE)
  expect_equal(nrow(cleaned[["ADT"]]), nrow(obj[["ADT"]]))
})

test_that(".adt_match_features resolves markers whose own name has punctuation", {
  feats <- c("Ab.HLA.DR", "Ab.HLA.ABC", "Ab.CD3-UCHT1", "Ab.TCR.AB", "HLADR", "Ab.CD4")

  expect_equal(sort(.adt_match_features(feats, "HLA-DR")),
               sort(c("Ab.HLA.DR", "HLADR")))
  expect_equal(.adt_match_features(feats, "TCR-AB"), "Ab.TCR.AB")

  expect_false("Ab.HLA.ABC" %in% .adt_match_features(feats, "HLA-DR"))
  expect_equal(.adt_match_features(feats, "CD3"), "Ab.CD3-UCHT1")
  expect_equal(.adt_match_features(c("ABCD3", "Ab.CD3"), "CD3"), "Ab.CD3")
})
