# ============================================================================
# test-portability.R
# Self-containment guard (portability.R): a foreign disk-backed layer whose
# backing cannot be read must be caught with INPUT_NOT_SELF_CONTAINED
# instead of dying deep inside the analysis (production failure:
# functional-state-profile inside UCell with "there is no package called
# 'BPCells'"). Outputs are bundles (disk_backed.R): any backed layer in the
# saved .rds must have its backing directory under <name>_matrices/.
# ============================================================================

.pt_sf_root <- local({
  cands <- c("portability.R", file.path("..", "portability.R"),
             file.path(dirname(dirname(FIXTURES)), "portability.R"))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("cannot locate portability.R from test wd")
  dirname(normalizePath(hit))
})
if (!exists("sf_materialize_layers", mode = "function")) {
  source(file.path(.pt_sf_root, "error_codes.R"))
  source(file.path(.pt_sf_root, "portability.R"))
}
if (!exists("sf_write_backed", mode = "function")) {
  source(file.path(.pt_sf_root, "disk_backed.R"))
}

fake_disk_matrix <- function() {
  # Class attribute mimicking a serialized BPCells matrix whose package is
  # NOT installed: detection must work from the class attribute alone,
  # without loading any namespace.
  structure(list(), class = structure("IterableMatrix", package = "BPCells"))
}

test_that(".sf_matrix_portable accepts base and Matrix classes only", {
  skip_if_not_installed("Matrix")
  expect_true(.sf_matrix_portable(matrix(1:4, 2)))
  expect_true(.sf_matrix_portable(
    Matrix::Matrix(matrix(c(0, 1, 0, 2), 2), sparse = TRUE)))
  expect_false(.sf_matrix_portable(fake_disk_matrix()))
})

test_that(".sf_assay_layers reads assays via direct slot access", {
  skip_if_not_installed("Seurat")
  obj <- make_tiny_seurat(n_cells = 10, n_genes = 20, add_adt = FALSE)
  layers <- .sf_assay_layers(obj@assays$RNA)
  expect_gte(length(layers), 1)
  expect_true(all(vapply(layers, .sf_matrix_portable, logical(1))))
})

test_that("sf_materialize_layers passes a self-contained object through unchanged", {
  skip_if_not_installed("Seurat")
  obj <- make_tiny_seurat(n_cells = 10, n_genes = 20, add_adt = FALSE)
  out <- sf_materialize_layers(obj, verbose = FALSE)
  expect_s4_class(out, "Seurat")
  before <- lapply(.sf_assay_layers(obj@assays$RNA), class)
  after <- lapply(.sf_assay_layers(out@assays$RNA), class)
  expect_identical(before, after)
})

test_that("non-Seurat input is returned unchanged", {
  x <- list(1, 2)
  expect_identical(sf_materialize_layers(x, verbose = FALSE), x)
})

# --- Subprocess test: the full main.R contract for a disk-backed input -------

.run_ld_main_pt <- function(cfg, in_dir, out_dir) {
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg_path <- file.path(in_dir, "config.json")
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), cfg_path)
  old <- setwd(.pt_sf_root); on.exit(setwd(old))
  rscript <- file.path(R.home("bin"), "Rscript")
  res <- suppressWarnings(system2(
    rscript, c("main.R", "--config", shQuote(cfg_path)),
    env = c(paste0("SF_INPUT_PATH=", in_dir),
               paste0("SF_OUTPUT_PATH=", out_dir)),
    stdout = TRUE, stderr = TRUE
  ))
  list(status = attr(res, "status") %||% 0L,
       output = paste(res, collapse = "\n"),
       error_json = file.path(out_dir, "error.json"))
}

test_that("a disk-backed .rds input fails INPUT_NOT_SELF_CONTAINED at the entry point", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("optparse")
  skip_if_not_installed("Seurat")
  in_dir <- tempfile("pt_in_"); out_dir <- tempfile("pt_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)

  # Inject a fake disk-backed layer: the class attribute of a serialized
  # BPCells matrix, with no BPCells needed to build or detect it.
  obj <- make_tiny_seurat(n_cells = 10, n_genes = 20, add_adt = FALSE)
  assay <- obj@assays$RNA
  fake <- structure(list(),
                     class = structure("IterableMatrix", package = "BPCells"))
  sn <- methods::slotNames(class(assay))
  if ("layers" %in% sn) {
    lay <- methods::slot(assay, "layers"); lay[["counts"]] <- fake
    methods::slot(assay, "layers", check = FALSE) <- lay
  } else {
    methods::slot(assay, "counts", check = FALSE) <- fake
  }
  obj@assays$RNA <- assay
  suppressWarnings(saveRDS(obj, file.path(in_dir, "data", "disk_backed.rds")))

  r <- .run_ld_main_pt(list(load_data_path = "data/disk_backed.rds",
                            project_name = "PortabilityTest", verbose = FALSE),
                       in_dir, out_dir)
  expect_identical(r$status, 1L)
  expect_true(file.exists(r$error_json),
              info = paste("no error.json; output was:\n", r$output))
  ej <- jsonlite::fromJSON(r$error_json)
  expect_identical(ej$error_code, "INPUT_NOT_SELF_CONTAINED")
  expect_identical(ej$step, "save-output")
  expect_true(any(grepl("RNA/counts", ej$details)))
  expect_true(any(grepl("dgCMatrix", ej$details)))
})

test_that("an .h5ad input yields a bundle output whose backed layers travel with it", {
  # The h5ad path opens X with BPCells::open_matrix_anndata_hdf5(), so the
  # Seurat object returned by load_data() carries IterableMatrix layers
  # pointing at the INPUT .h5ad. The save re-writes every such layer into
  # the output's own <name>_matrices/ directory. This test pins the
  # contract: any backed layer in the saved .rds has its backing directory
  # next to the .rds, never a path into the input.
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("optparse")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("BPCells")
  skip_if_no_py_anndata()
  in_dir <- tempfile("pt_h5_in_"); out_dir <- tempfile("pt_h5_out_")
  on.exit(unlink(c(in_dir, out_dir), recursive = TRUE))
  dir.create(file.path(in_dir, "data"), recursive = TRUE)
  .write_test_h5ad(file.path(in_dir, "data", "tiny.h5ad"),
                   n_cells = 10L, n_genes = 6L)

  r <- .run_ld_main_pt(list(load_data_path = "data/tiny.h5ad",
                            project_name = "PortabilityH5ad",
                            verbose = FALSE),
                       in_dir, out_dir)
  if (r$status != 0L) print(r$output)
  expect_identical(r$status, 0L)

  out_rds <- list.files(file.path(out_dir, "output"),
                        pattern = "\\.rds$", full.names = TRUE)
  expect_length(out_rds, 1)
  obj <- readRDS(out_rds[1])
  expect_s4_class(obj, "Seurat")
  root <- .sf_backing_root(out_rds[1])
  for (a in names(obj@assays)) {
    layers <- .sf_assay_layers(obj@assays[[a]])
    for (ln in names(layers)) {
      m <- layers[[ln]]
      if (.sf_matrix_portable(m)) next
      expect_identical(attr(class(m), "package"), "BPCells",
                       info = sprintf("assay %s/%s: unexpected class %s",
                                      a, ln, class(m)[1]))
      expect_true(dir.exists(.sf_backing_dir(root, a, ln)),
                  info = sprintf("no backing directory for %s/%s", a, ln))
    }
  }
  reloaded <- sf_load_backed(out_rds[1], verbose = FALSE)
  expect_s4_class(reloaded, "Seurat")
  expect_identical(dim(reloaded), dim(obj))
})
