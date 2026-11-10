#' Backing-directory root for a bundle .rds
#'
#' @param rds_path Path to a boundary `.rds` file.
#' @return Path of the sibling `<rds basename>_matrices` directory that holds
#'   the BPCells backing matrices for that `.rds`.
.sf_backing_root <- function(rds_path)
  paste0(sub("\\.rds$", "", rds_path, ignore.case = TRUE), "_matrices")

#' Backing directory for one assay layer
#'
#' @param root Backing root returned by [.sf_backing_root()].
#' @param assay Assay name (e.g. `"RNA"`).
#' @param layer Layer name (e.g. `"counts"`).
#' @return Directory under `root` holding that layer's BPCells matrix.
#'   Characters outside `[A-Za-z0-9_.-]` are replaced with `_` so the name is
#'   portable across filesystems.
.sf_backing_dir <- function(root, assay, layer)
  file.path(root, gsub("[^A-Za-z0-9_.-]", "_", paste0(assay, ".", layer)))

#' Test whether a matrix layer is BPCells-backed, from the class attribute alone.
#'
#' @param m A matrix layer of any class.
#' @return `TRUE` when the layer's class comes from the BPCells package.
.sf_bpcells_layer <- function(m) {
  pkg <- attr(class(m), "package")
  !is.null(pkg) && identical(pkg, "BPCells")
}

#' Write a Seurat object as a boundary bundle
#'
#' Writes `obj` to `rds_path` with every sparse assay layer stored as a
#' BPCells matrix directory under `<rds basename>_matrices/`, so the saved
#' `.rds` carries only small stubs and downstream smart functions can
#' memory-map the matrix values. Dense base-matrix layers stay inside the
#' `.rds`. Slot-based (v3/v4) assays are materialized to in-memory matrices
#' instead. A disk-backed input layer whose backing cannot be read fails via
#' `sf_fail("INPUT_NOT_SELF_CONTAINED")`. When BPCells is not installed the
#' function falls back to a fully materialized plain `.rds`.
#'
#' @param obj Seurat object to write; a non-Seurat object is written as-is
#'   after the materialize guard.
#' @param rds_path Output `.rds` path; the backing directory is created next
#'   to it.
#' @param verbose Logical; print one line per layer written.
#' @return Invisibly, the (possibly modified) object that was serialized.
sf_write_backed <- function(obj, rds_path, verbose = TRUE) {
  backable <- inherits(obj, "Seurat") &&
    requireNamespace("BPCells", quietly = TRUE)
  if (!backable) {
    obj <- sf_materialize_layers(obj, input_path = NULL, verbose = verbose)
  } else {
    root <- .sf_backing_root(rds_path)
    for (a in names(obj@assays)) {
      assay <- obj@assays[[a]]
      sn <- methods::slotNames(class(assay))
      if (!("layers" %in% sn)) {
        changed <- FALSE
        for (sl in intersect(c("counts", "data", "scale.data"), sn)) {
          m <- methods::slot(assay, sl)
          dm <- dim(m)
          if (is.null(m) || length(dm) != 2 || any(dm == 0)) next
          if (.sf_matrix_portable(m)) next
          mem <- tryCatch(methods::as(m, "dgCMatrix"), error = function(e) {
            sf_fail("INPUT_NOT_SELF_CONTAINED", .SF_PORTABILITY_MESSAGE,
                    details = c(
                      sprintf("disk-backed layer - %s/%s: %s", a, sl, class(m)[1]),
                      sprintf("materialization failed: %s", conditionMessage(e)),
                      .SF_PORTABILITY_REMEDIATION))
          })
          methods::slot(assay, sl, check = FALSE) <- mem
          changed <- TRUE
        }
        if (changed) obj@assays[[a]] <- assay
        next
      }
      lay <- methods::slot(assay, "layers")
      if (is.null(names(lay))) next
      changed <- FALSE
      for (ln in names(lay)) {
        m <- lay[[ln]]
        if (is.null(m) || is.matrix(m)) next
        if (!.sf_bpcells_layer(m) && !.sf_matrix_portable(m)) {
          m <- tryCatch(methods::as(m, "dgCMatrix"), error = function(e) {
            sf_fail("INPUT_NOT_SELF_CONTAINED", .SF_PORTABILITY_MESSAGE,
                    details = c(
                      sprintf("disk-backed layer - %s/%s: %s", a, ln, class(m)[1]),
                      sprintf("materialization failed: %s", conditionMessage(e)),
                      .SF_PORTABILITY_REMEDIATION))
          })
        } else if (!.sf_bpcells_layer(m) && !inherits(m, "dgCMatrix")) {
          m2 <- tryCatch(methods::as(m, "dgCMatrix"), error = function(e) NULL)
          if (is.null(m2)) next
          m <- m2
        }
        d <- .sf_backing_dir(root, a, ln)
        if (dir.exists(d)) unlink(d, recursive = TRUE)
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        lay[[ln]] <- tryCatch(
          BPCells::write_matrix_dir(m, dir = d, compress = TRUE,
                                    overwrite = TRUE),
          error = function(e) {
            if (.sf_bpcells_layer(m))
              sf_fail("INPUT_NOT_SELF_CONTAINED", .SF_PORTABILITY_MESSAGE,
                      details = c(
                        sprintf("disk-backed layer - %s/%s: %s",
                                a, ln, class(m)[1]),
                        sprintf("backing read failed: %s", conditionMessage(e)),
                        .SF_PORTABILITY_REMEDIATION))
            stop(e)
          })
        changed <- TRUE
        if (verbose) cat(sprintf("  Backed %s/%s -> %s/%s\n",
                                 a, ln, basename(root), basename(d)))
      }
      if (changed) {
        methods::slot(assay, "layers", check = FALSE) <- lay
        obj@assays[[a]] <- assay
      }
    }
  }
  con <- gzfile(rds_path, open = "wb", compression = 1)
  closed <- FALSE
  on.exit(if (!closed) try(close(con), silent = TRUE), add = TRUE)
  saveRDS(obj, file = con)
  close(con)
  closed <- TRUE
  invisible(obj)
}

#' Read a boundary .rds, reattaching its backed layers
#'
#' Reads `rds_path` and, when a `<rds basename>_matrices/` directory sits
#' next to it, reattaches each BPCells layer to its backing directory by
#' convention (embedded paths from the producing container are ignored).
#' Layers whose name matches `keep_backed` stay disk-backed; all others are
#' materialized to in-memory `dgCMatrix`. A plain `.rds` with no backing
#' directory goes through `sf_materialize_layers()` unchanged, preserving the
#' original self-containment guard for foreign uploads.
#'
#' @param rds_path Path to the input `.rds`.
#' @param keep_backed Regular expression matched against layer names; matching
#'   layers stay disk-backed. `NULL` or `""` materializes every layer.
#' @param verbose Logical; print one line per reattached or materialized
#'   layer.
#' @return The loaded Seurat object. Exits via `sf_fail()` when the file
#'   cannot be read (`INPUT_LOAD_FAILED`), a backing directory is missing or
#'   unreadable (`INPUT_NOT_SELF_CONTAINED`), or BPCells is not installed in
#'   the image.
sf_load_backed <- function(rds_path, keep_backed = "^counts", verbose = TRUE) {
  obj <- tryCatch(
    readRDS(rds_path),
    error = function(e) {
      sf_fail("INPUT_LOAD_FAILED",
              "The input Seurat .rds file could not be read.",
              details = c(sprintf("path: %s", rds_path), e$message))
    })
  root <- .sf_backing_root(rds_path)
  if (!inherits(obj, "Seurat") || !dir.exists(root))
    return(sf_materialize_layers(obj, input_path = rds_path, verbose = verbose))
  if (!requireNamespace("BPCells", quietly = TRUE))
    sf_fail("INPUT_NOT_SELF_CONTAINED",
            paste("The input carries BPCells-backed matrix layers but this",
                  "image does not have the BPCells package installed."),
            details = sprintf("backing directory: %s", root))
  for (a in names(obj@assays)) {
    assay <- obj@assays[[a]]
    if (!("layers" %in% methods::slotNames(class(assay)))) next
    lay <- methods::slot(assay, "layers")
    if (is.null(names(lay))) next
    changed <- FALSE
    for (ln in names(lay)) {
      m <- lay[[ln]]
      if (is.null(m) || !.sf_bpcells_layer(m)) next
      d <- .sf_backing_dir(root, a, ln)
      m2 <- tryCatch({
        x <- BPCells::open_matrix_dir(d)
        stopifnot(identical(dim(x), dim(m)))
        if (all(dim(x) > 0)) invisible(methods::as(x[1, 1], "dgCMatrix"))
        x
      }, error = function(e) {
        sf_fail("INPUT_NOT_SELF_CONTAINED",
                paste("A disk-backed matrix layer could not be reattached to",
                      "its backing directory next to the input .rds."),
                details = c(
                  sprintf("layer: %s/%s", a, ln),
                  sprintf("expected backing directory: %s", d),
                  conditionMessage(e),
                  paste("the upstream dataset must transfer the",
                        basename(root), "directory alongside the .rds")))
      })
      keep <- !is.null(keep_backed) && nzchar(keep_backed) &&
        grepl(keep_backed, ln)
      if (keep) {
        lay[[ln]] <- m2
        if (verbose) cat(sprintf("  Reattached %s/%s (stays disk-backed)\n", a, ln))
      } else {
        lay[[ln]] <- methods::as(m2, "dgCMatrix")
        if (verbose) cat(sprintf("  Loaded %s/%s into memory from backing\n", a, ln))
      }
      changed <- TRUE
    }
    if (changed) {
      methods::slot(assay, "layers", check = FALSE) <- lay
      obj@assays[[a]] <- assay
    }
  }
  obj
}

#' Copy a boundary bundle
#'
#' Copies a boundary `.rds` together with its backing directory, renaming the
#' backing root to match the destination `.rds` so convention-based
#' reattachment keeps working downstream. A plain `.rds` (no backing
#' directory) is copied alone.
#'
#' @param src_rds Source `.rds` path.
#' @param dst_rds Destination `.rds` path.
#' @return `TRUE` on success, `FALSE` on any copy failure.
sf_copy_bundle <- function(src_rds, dst_rds) {
  ok <- tryCatch(file.copy(src_rds, dst_rds, overwrite = TRUE),
                 error = function(e) FALSE)
  if (!isTRUE(ok)) return(FALSE)
  src_root <- .sf_backing_root(src_rds)
  if (!dir.exists(src_root)) return(TRUE)
  dst_root <- .sf_backing_root(dst_rds)
  if (dir.exists(dst_root)) unlink(dst_root, recursive = TRUE)
  if (!dir.create(dst_root, recursive = TRUE, showWarnings = FALSE)) return(FALSE)
  entries <- list.files(src_root, full.names = TRUE)
  if (length(entries) == 0) return(TRUE)
  isTRUE(all(file.copy(entries, dst_root, recursive = TRUE)))
}
