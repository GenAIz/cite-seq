.ORCH_PATH <- local({
  cands <- c("orchestrator.R", file.path("..", "orchestrator.R"))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("cannot locate orchestrator.R from test wd")
  normalizePath(hit)
})
.ORCH_ROOT  <- dirname(.ORCH_PATH)
.orch_exprs <- parse(.ORCH_PATH)

.orch_rhs <- function(name) {
  for (e in .orch_exprs) {
    if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
        is.name(e[[2]]) && identical(as.character(e[[2]]), name)) {
      return(e[[3]])
    }
  }
  stop("no top-level definition of '", name, "' in orchestrator.R")
}

.orch_env <- function(config, run_timestamp = "20260101_000000") {
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(a, b) if (!is.null(a)) a else b
  env$config <- config
  env$run_timestamp <- run_timestamp
  for (fn in c("step_output_path", "step_resume_path", "load_previous_output")) {
    assign(fn, eval(.orch_rhs(fn), env), envir = env)
  }
  env
}

.dispatch_steps <- names(as.list(.orch_rhs("step_dispatch")))[-1]

test_that("the dispatch table registers exactly the documented steps", {
  expect_setequal(.dispatch_steps,
                  c("load_data", "add_metadata", "process_seurat",
                    "cell_type_annotation", "aggregate_ml", "export_loupe"))
})

test_that("every dispatch entry maps to a run function defined in orchestrator.R", {
  entries <- as.list(.orch_rhs("step_dispatch"))[-1]
  for (step in names(entries)) {
    fn_name <- as.character(entries[[step]])
    expect_identical(fn_name, paste0("run_", step),
                     info = sprintf("dispatch naming convention broken for '%s'", step))
    rhs <- .orch_rhs(fn_name)
    expect_true(is.call(rhs) && identical(as.character(rhs[[1]]), "function"),
                info = sprintf("'%s' is not a function definition", fn_name))
  }
})

test_that("default_steps are the canonical five, all dispatchable, load_data first", {
  default_steps <- eval(.orch_rhs("default_steps"), baseenv())
  expect_identical(default_steps,
                   c("load_data", "add_metadata", "process_seurat",
                     "cell_type_annotation", "aggregate_ml"))
  expect_true(all(default_steps %in% .dispatch_steps))
  expect_identical(default_steps[1], "load_data")

  expect_false("export_loupe" %in% default_steps)
})

test_that("resume from a step filters to remaining steps", {
  steps <- c("load_data", "process_seurat", "add_metadata", "aggregate_ml")
  resume_index <- which(steps == "add_metadata")
  expect_equal(steps[resume_index:length(steps)], c("add_metadata", "aggregate_ml"))
})

test_that("resume respects custom step order", {
  steps <- c("load_data", "cell_type_annotation", "add_metadata", "aggregate_ml")
  resume_index <- which(steps == "add_metadata")
  expect_equal(steps[resume_index:length(steps)], c("add_metadata", "aggregate_ml"))
})

test_that("resume from step not in list yields empty index", {
  steps <- c("load_data", "process_seurat")
  expect_length(which(steps == "add_metadata"), 0L)
})

test_that("step_output_path generates timestamped filenames", {
  tmp_dir <- tempfile("orch_out_"); dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  env <- .orch_env(config = list(output_files = list(load_data = "my_output.rds"),
                                 paths = list(output_dir = tmp_dir)),
                   run_timestamp = "20260531_120000")

  expect_match(env$step_output_path("load_data", "load_data.rds"),
               "my_output_20260531_120000\\.rds$")
  expect_match(env$step_output_path("add_metadata", "add_metadata.rds"),
               "add_metadata_20260531_120000\\.rds$")
})

test_that("step_resume_path finds most recent timestamped file", {
  tmp_dir <- tempfile("orch_resume_"); dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  writeLines("old", file.path(tmp_dir, "add_metadata_20260101_010000.rds"))
  Sys.sleep(0.1)
  writeLines("new", file.path(tmp_dir, "add_metadata_20260531_120000.rds"))

  env <- .orch_env(config = list(output_files = list(),
                                 paths = list(output_dir = tmp_dir)))
  expect_match(basename(env$step_resume_path("add_metadata", "add_metadata.rds")),
               "20260531_120000")
})

test_that("step_resume_path falls back to non-timestamped file", {
  tmp_dir <- tempfile("orch_fallback_"); dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  env <- .orch_env(config = list(output_files = list(),
                                 paths = list(output_dir = tmp_dir)))
  expect_match(basename(env$step_resume_path("add_metadata", "add_metadata.rds")),
               "^add_metadata\\.rds$")
})

.plant <- function(dir, ...) {
  for (f in c(...)) writeLines("x", file.path(dir, f))
}
.resume_env <- function(tmp_dir) {
  .orch_env(config = list(output_files = list(),
                          paths = list(output_dir = tmp_dir)))
}

test_that("a downstream artefact is never handed to an upstream step", {
  tmp <- tempfile("orch_ds_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))

  .plant(tmp, "cell_type_annotation_20260101_000000.rds")
  env <- .resume_env(tmp)
  expect_null(env$load_previous_output(1, "process_seurat"))
})

test_that("the nearest upstream artefact wins over more distant ones", {
  tmp <- tempfile("orch_near_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  .plant(tmp, "load_data_20260101_000000.rds",
              "process_seurat_20260101_000000.rds")
  env <- .resume_env(tmp)
  hit <- env$load_previous_output(1, "cell_type_annotation")
  expect_match(basename(hit), "^process_seurat_")
})

test_that("this run's prior steps are searched before older artefacts", {
  tmp <- tempfile("orch_prior_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  .plant(tmp, "load_data_20260101_000000.rds")
  env <- .resume_env(tmp)
  hit <- env$load_previous_output(2, c("load_data", "process_seurat"))
  expect_match(basename(hit), "^load_data_")
})

test_that("the first step of the pipeline has nothing to resume from", {
  tmp <- tempfile("orch_first_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  .plant(tmp, "cell_type_annotation_20260101_000000.rds")
  env <- .resume_env(tmp)
  expect_null(env$load_previous_output(1, "load_data"))
})

test_that("add_metadata placed late in the run resumes from the annotated object", {

  tmp <- tempfile("orch_pos_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  .plant(tmp, "load_data_20260101_000000.rds",
              "cell_type_annotation_20260101_000000.rds")
  env <- .resume_env(tmp)
  hit <- env$load_previous_output(
    4, c("load_data", "process_seurat", "cell_type_annotation", "add_metadata"))
  expect_match(basename(hit), "^cell_type_annotation_")
})

.run_orchestrator_dryrun <- function(steps) {
  in_dir  <- tempfile("orch_cfg_");  dir.create(in_dir)
  out_dir <- file.path(in_dir, "out"); log_dir <- file.path(in_dir, "logs")
  cfg <- list(
    metadata    = list(project_name = "orch_test"),
    paths       = list(output_dir = out_dir, logs_dir = log_dir,
                       input_data = list()),
    execution   = list(mode = "step-by-step", steps = steps, n_cores = 1),
    load_data   = list(),
    processing  = list(),
    aggregation = list()
  )
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), cfg_path)

  old <- setwd(.ORCH_ROOT); on.exit({setwd(old); unlink(in_dir, recursive = TRUE)})
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("orchestrator.R", "--config", shQuote(cfg_path), "--dry-run"),
    env = paste0("SF_LOG_DIR=", log_dir),
    stdout = TRUE, stderr = TRUE
  ))
  list(status = attr(res, "status") %||% 0L, output = paste(res, collapse = "\n"))
}

skip_orch_subprocess <- function() {
  skip_if_not_installed("optparse")
  skip_if_not_installed("jsonlite")

  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
}

test_that("a valid step list passes --dry-run", {
  skip_orch_subprocess()
  r <- .run_orchestrator_dryrun(c("load_data", "add_metadata"))
  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)
  expect_match(r$output, "Dry-run")
})

test_that("unknown step names are rejected by the orchestrator", {
  skip_orch_subprocess()
  r <- .run_orchestrator_dryrun(c("load_data", "bogus_step"))
  expect_identical(r$status, 1L)
  expect_match(r$output, "Unknown steps")
})

test_that("load_data at non-first position is rejected by the orchestrator", {
  skip_orch_subprocess()
  r <- .run_orchestrator_dryrun(c("process_seurat", "load_data"))
  expect_identical(r$status, 1L)
  expect_match(r$output, "load_data must be the first step")
})
