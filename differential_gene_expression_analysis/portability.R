#' TRUE when a layer is a plain in-memory matrix (base or Matrix package) readable in any container.
.sf_matrix_portable <- function(m) {
  if (is.matrix(m)) return(TRUE)
  pkg <- attr(class(m), "package")
  !is.null(pkg) && identical(pkg, "Matrix")
}

#' Named list of an assay's matrix layers via direct slot access (Seurat v5 Assay5 or v3/v4 Assay).
.sf_assay_layers <- function(assay) {
  sn <- methods::slotNames(class(assay))
  if ("layers" %in% sn) {
    lay <- methods::slot(assay, "layers")
    Filter(function(m) !is.null(m), as.list(lay))
  } else {
    out <- list()
    for (sl in intersect(c("counts", "data", "scale.data"), sn)) {
      m <- methods::slot(assay, sl)
      d <- dim(m)
      if (!is.null(m) && length(d) == 2 && all(d > 0)) out[[sl]] <- m
    }
    out
  }
}

.SF_PORTABILITY_MESSAGE <- paste(
  "The input .rds is not self-contained: it carries disk-backed matrix",
  "layers whose backing files are not part of the upload.")

.SF_PORTABILITY_REMEDIATION <- c(
  paste("The matrix values live in backing files on the computer where the",
        "object was created; those files are not inside the .rds and cannot",
        "travel with it."),
  paste("Re-export the object with in-memory matrices before uploading -",
        "in the session that created it, convert each disk-backed layer,",
        "e.g. obj[[\"RNA\"]]$counts <- as(obj[[\"RNA\"]]$counts,",
        "\"dgCMatrix\") - or provide the original .h5ad/.h5mu instead."))

#' Materialize every assay layer to an in-memory matrix, or exit via sf_fail() if one cannot be; non-Seurat input is returned unchanged.
sf_materialize_layers <- function(obj, input_path = NULL, verbose = TRUE) {
  if (!inherits(obj, "Seurat")) return(obj)
  for (a in names(obj@assays)) {
    assay  <- obj@assays[[a]]
    layers <- tryCatch(.sf_assay_layers(assay), error = function(e) list())
    for (ln in names(layers)) {
      m <- layers[[ln]]
      if (.sf_matrix_portable(m)) next
      cls <- class(m)[1]
      mem <- tryCatch(methods::as(m, "dgCMatrix"), error = function(e) {
        sf_fail("INPUT_NOT_SELF_CONTAINED", .SF_PORTABILITY_MESSAGE,
                details = c(
                  if (!is.null(input_path)) sprintf("input: %s", input_path),
                  sprintf("disk-backed layer - %s/%s: %s", a, ln, cls),
                  sprintf("materialization failed: %s", conditionMessage(e)),
                  .SF_PORTABILITY_REMEDIATION))
      })
      if (verbose) cat(sprintf(
        "  Materialized %s/%s (%s -> dgCMatrix) so the object is self-contained\n",
        a, ln, cls))
      sn <- methods::slotNames(class(assay))
      if ("layers" %in% sn) {
        lay <- methods::slot(assay, "layers")
        lay[[ln]] <- mem
        methods::slot(assay, "layers", check = FALSE) <- lay
      } else {
        methods::slot(assay, ln, check = FALSE) <- mem
      }
    }
    obj@assays[[a]] <- assay
  }
  obj
}
