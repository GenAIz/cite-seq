three_cond_dge <- function(seed = 42) {
  make_dge_seurat(n_cells = 240, n_genes = 60, n_samples = 12,
                  n_timepoints = 2, n_clusters = 2, seed = seed,
                  condition_names = c("Control", "TreatA", "TreatB"))
}

test_that("default contrasts produce all pairwise comparisons", {
  skip_if_not_installed("limma")

  result <- perform_dge_analysis(
    seurat_obj      = three_cond_dge(),
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    verbose         = FALSE
  )

  expect_gt(nrow(result), 0)

  expect_setequal(unique(result$contrast),
                  c("Control - TreatA", "Control - TreatB", "TreatA - TreatB"))
})

test_that("reference_condition produces contrasts against the reference only", {
  skip_if_not_installed("limma")

  result <- perform_dge_analysis(
    seurat_obj          = three_cond_dge(),
    test.use            = "limma",
    pseudobulking       = TRUE,
    design_strategy     = "stratified",
    cluster_col         = "CellType",
    sample_col          = "SampleID",
    condition_col       = "Condition",
    timepoint_col       = "Timepoint",
    timepoints_to_analyze = "T1",
    reference_condition = "Control",
    verbose             = FALSE
  )

  expect_gt(nrow(result), 0)

  expect_setequal(unique(result$contrast),
                  c("TreatA - Control", "TreatB - Control"))
})

test_that("ident.1/ident.2 produces a single specific contrast", {
  skip_if_not_installed("limma")

  result <- perform_dge_analysis(
    seurat_obj      = three_cond_dge(),
    test.use        = "limma",
    pseudobulking   = TRUE,
    design_strategy = "stratified",
    cluster_col     = "CellType",
    sample_col      = "SampleID",
    condition_col   = "Condition",
    timepoint_col   = "Timepoint",
    timepoints_to_analyze = "T1",
    ident.1         = "TreatA",
    ident.2         = "Control",
    verbose         = FALSE
  )

  expect_gt(nrow(result), 0)

  expect_identical(unique(result$contrast), "TreatA - Control")
})
