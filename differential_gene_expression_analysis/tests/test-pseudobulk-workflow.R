test_that("limma stratified workflow returns a data frame with expected columns", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj      = obj,
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    ident.1         = "Treated",
    ident.2         = "Control",
    verbose         = FALSE
  )

  expect_s3_class(result, "data.frame")

  expect_gt(nrow(result), 0)
  for (col in c("gene", "cluster", "timepoint", "contrast",
                "logFC", "P.Value", "adj.P.Val")) {
    expect_true(col %in% colnames(result),
                info = sprintf("column '%s' missing", col))
  }

  expect_identical(unique(result$contrast), "Treated - Control")
  expect_identical(unique(as.character(result$timepoint)), "T1")

  top20 <- head(result$gene[order(result$P.Value)], 20)
  expect_gte(sum(top20 %in% paste0("GENE", sprintf("%04d", 1:20))), 12)
})

test_that("edgeR workflow returns results with expected columns", {
  skip_if_not_installed("edgeR")

  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj      = obj,
    test.use        = "edgeR",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    ident.1         = "Treated",
    ident.2         = "Control",
    verbose         = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  for (col in c("gene", "logFC", "PValue", "FDR")) {
    expect_true(col %in% colnames(result),
                info = sprintf("column '%s' missing", col))
  }
  expect_identical(unique(result$contrast), "Treated - Control")
})

test_that("direct comparison mode works when timepoints_to_analyze is NULL", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj          = obj,
    test.use            = "limma",
    pseudobulking       = TRUE,
    design_strategy     = "stratified",
    cluster_col         = "CellType",
    sample_col          = "SampleID",
    condition_col       = "Condition",
    timepoint_col       = "Timepoint",
    timepoints_to_analyze = NULL,
    ident.1             = "Treated",
    ident.2             = "Control",
    verbose             = FALSE
  )

  expect_s3_class(result, "data.frame")

  expect_gt(nrow(result), 0)
  expect_identical(unique(result$contrast), "Treated - Control")

  expect_false("timepoint" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))
})

test_that("interaction design produces per-timepoint contrasts", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj <- make_dge_seurat(n_cells = 120, n_samples = 6,
                          n_conditions = 2, n_timepoints = 2, n_clusters = 2)

  result <- perform_dge_analysis(
    seurat_obj      = obj,
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "interaction",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    participant_col = "ParticipantID",
    ident.1         = "Treated",
    ident.2         = "Control",
    verbose         = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  expect_setequal(unique(result$contrast),
                  c("T1_Treated - T1_Control", "T2_Treated - T2_Control"))
})

test_that("multiple input objects produces results with object column", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj1 <- make_mini_dge_seurat(seed = 42)
  obj2 <- make_mini_dge_seurat(seed = 99)

  result <- perform_dge_analysis(
    seurat_obj      = list(Batch1 = obj1, Batch2 = obj2),
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    ident.1         = "Treated",
    ident.2         = "Control",
    verbose         = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("object" %in% colnames(result))

  expect_setequal(unique(result$object), c("Batch1", "Batch2"))
})

test_that("an invalid timepoints_to_analyze value is a hard validation error", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj <- make_mini_dge_seurat()
  expect_error(
    perform_dge_analysis(
      seurat_obj      = obj,
      test.use        = "limma",
      pseudobulking   = TRUE,
      design_strategy = "stratified",
      cluster_col     = "CellType",
      sample_col      = "SampleID",
      condition_col   = "Condition",
      timepoint_col   = "Timepoint",
      timepoints_to_analyze = "T99",
      verbose         = FALSE
    ),
    regexp = "Invalid timepoints_to_analyze"
  )
})

test_that("a run where no contrast is estimable returns an empty data frame", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  obj <- make_mini_dge_seurat()

  result <- perform_dge_analysis(
    seurat_obj      = obj,
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    min_samples_per_group = 999,
    verbose         = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 0L)
})
