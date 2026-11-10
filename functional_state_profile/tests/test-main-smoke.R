.smoke_sf_root <- local({
  cands <- c("main.R", file.path("..", "main.R"),
             file.path(dirname(dirname(FIXTURES)), "main.R"))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("cannot locate main.R from test wd")
  dirname(normalizePath(hit))
})

skip_main <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("UCell")
  skip_if_not_installed("BiocParallel")
  skip_if_not_installed("optparse")
  skip_if_not_installed("jsonlite")
}

test_that("main.R runs a config end-to-end and exports dynamic state columns", {
  skip_main()

  sf_root <- .smoke_sf_root
  in_dir  <- file.path(tempdir(), "fsp_smoke_in")
  out_dir <- file.path(tempdir(), "fsp_smoke_out")
  unlink(c(in_dir, out_dir), recursive = TRUE)
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE), add = TRUE)

  saveRDS(make_fsp_seurat(), file.path(in_dir, "data", "annotated_seurat.rds"))
  file.copy(TEST_MARKERS_JSON, file.path(in_dir, "test_markers.json"))

  cfg <- list(
    input_rds        = "data/annotated_seurat.rds",
    output_rds       = "smoke_seurat.rds",
    output_csv       = "smoke_scores.csv",
    marker_json_path = "test_markers.json",
    use_hvg          = FALSE,
    maxRank          = 100,
    min_score        = 0.0,
    run_validation   = TRUE,
    verbose          = FALSE
  )
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 8),
             cfg_path)

  old <- setwd(sf_root); on.exit(setwd(old), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", out_dir)),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(res, "status") %||% 0L
  expect_identical(status, 0L)
  if (status != 0L) print(tail(res, 30))

  rds_out <- list.files(file.path(out_dir, "output"),
                        pattern = "^smoke_seurat_.*\\.rds$", full.names = TRUE)
  csv_out <- list.files(file.path(out_dir, "output"),
                        pattern = "^smoke_scores_.*\\.csv$", full.names = TRUE)
  expect_length(rds_out, 1)
  expect_length(csv_out, 1)

  scores <- read.csv(csv_out[1], check.names = FALSE)

  for (col in c("cell_barcode", "dominant_state", "dominant_state_score",
                "dominant_state_margin", "exhaustion", "naive-like",
                "cytotoxicity")) {
    expect_true(col %in% colnames(scores),
                info = sprintf("column '%s' missing from score CSV", col))
  }
  expect_equal(nrow(scores), 50)
  expect_true(all(scores$dominant_state %in%
                  c("exhaustion", "naive-like", "cytotoxicity")))

  expect_false(file.exists(file.path(out_dir, "error.json")))
})

test_that("main.R rejects an out-of-range min_score at config validation", {
  skip_main()

  sf_root <- .smoke_sf_root
  in_dir  <- file.path(tempdir(), "fsp_smoke_badcfg")
  unlink(in_dir, recursive = TRUE); dir.create(in_dir, recursive = TRUE)
  on.exit(unlink(in_dir, recursive = TRUE), add = TRUE)

  dir.create(file.path(in_dir, "data"))
  saveRDS(1L, file.path(in_dir, "data", "annotated_seurat.rds"))

  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(
    list(input_rds = "data/annotated_seurat.rds", min_score = 1.5),
    auto_unbox = TRUE), cfg_path)

  old <- setwd(sf_root); on.exit(setwd(old), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", tempdir())),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(attr(res, "status") %||% 0L, 1L)
  hit <- any(grepl("min_score", res))
  if (!hit) print(res)
  expect_true(hit)
})

test_that("main.R reports a missing marker_json_path at config validation", {
  skip_main()

  sf_root <- .smoke_sf_root
  in_dir  <- file.path(tempdir(), "fsp_smoke_badmarkers")
  unlink(in_dir, recursive = TRUE); dir.create(in_dir, recursive = TRUE)
  on.exit(unlink(in_dir, recursive = TRUE), add = TRUE)

  dir.create(file.path(in_dir, "data"))
  saveRDS(1L, file.path(in_dir, "data", "annotated_seurat.rds"))

  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(
    list(input_rds        = "data/annotated_seurat.rds",
         marker_json_path = "markers/does_not_exist.json"),
    auto_unbox = TRUE), cfg_path)

  old <- setwd(sf_root); on.exit(setwd(old), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", tempdir())),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(attr(res, "status") %||% 0L, 1L)
  hit <- any(grepl("marker_json_path", res))
  if (!hit) print(res)
  expect_true(hit)
})

test_that("input auto-discovery finds an .rds NESTED under the data port", {

  skip_main()

  sf_root <- .smoke_sf_root
  in_dir  <- file.path(tempdir(), "fsp_nested_in")
  out_dir <- file.path(tempdir(), "fsp_nested_out")
  unlink(c(in_dir, out_dir), recursive = TRUE)
  dir.create(file.path(in_dir, "data", "output"), recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE), add = TRUE)

  saveRDS(make_fsp_seurat(),
          file.path(in_dir, "data", "output", "annotated_seurat.rds"))
  file.copy(TEST_MARKERS_JSON, file.path(in_dir, "test_markers.json"))

  cfg <- list(

    output_rds       = "nested_seurat.rds",
    output_csv       = "nested_scores.csv",
    marker_json_path = "test_markers.json",
    use_hvg          = FALSE,
    maxRank          = 100,
    run_validation   = FALSE,
    verbose          = FALSE
  )
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 8),
             cfg_path)

  old <- setwd(sf_root); on.exit(setwd(old), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", out_dir)),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(res, "status") %||% 0L
  if (status != 0L) print(tail(res, 30))
  expect_identical(status, 0L)

  expect_true(any(grepl("using newest .rds on the data port", res, fixed = TRUE)))
  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^nested_seurat_.*\\.rds$"), 1)
})

test_that("a marker JSON on the attachment port is discovered by name", {
  skip_main()

  sf_root <- .smoke_sf_root
  in_dir  <- file.path(tempdir(), "fsp_smoke_att_in")
  out_dir <- file.path(tempdir(), "fsp_smoke_att_out")
  unlink(c(in_dir, out_dir), recursive = TRUE)
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  dir.create(file.path(in_dir, "attachment"), recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE), add = TRUE)

  saveRDS(make_fsp_seurat(), file.path(in_dir, "data", "annotated_seurat.rds"))
  file.copy(TEST_MARKERS_JSON, file.path(in_dir, "attachment", "test_markers.json"))

  cfg <- list(
    input_rds        = "data/annotated_seurat.rds",
    output_rds       = "att_seurat.rds",
    output_csv       = "att_scores.csv",
    marker_json_path = "test_markers.json",
    use_hvg          = FALSE,
    maxRank          = 100,
    run_validation   = FALSE,
    verbose          = FALSE
  )
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 8),
             cfg_path)

  old <- setwd(sf_root); on.exit(setwd(old), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env    = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", out_dir)),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(res, "status") %||% 0L
  expect_identical(status, 0L)
  if (status != 0L) print(tail(res, 30))

  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^att_seurat_.*\\.rds$"), 1)
  expect_length(list.files(file.path(out_dir, "output"),
                           pattern = "^att_scores_.*\\.csv$"), 1)
  expect_false(file.exists(file.path(out_dir, "error.json")))
})
