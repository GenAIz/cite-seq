.loupe_reencode_barcodes <- function(bc) {
  seqs <- sub("-.*$", "", bc)
  sufs <- sub("^[^-]*", "", bc)
  new_bc <- paste0(seqs, "-", as.integer(factor(sufs)))
  if (anyDuplicated(new_bc)) return(NULL)
  new_bc
}

.loupe_cluster_columns <- function(meta, max_groups = 200L) {
  keep <- character(0)
  for (cc in colnames(meta)) {
    v <- meta[[cc]]
    cat_like <- is.factor(v) || is.character(v) ||
                (is.numeric(v) && all(v == round(v), na.rm = TRUE))
    if (!cat_like) next
    nu <- length(unique(v[!is.na(v)]))
    if (nu >= 2L && nu <= max_groups) keep <- c(keep, cc)
  }
  keep
}

export_seurat_to_loupe <- function(seurat_obj,
                                   output_dir,
                                   output_name,
                                   assays = c("RNA", "ADT"),
                                   executable_path = NULL,
                                   verbose = TRUE) {
  if (!requireNamespace("loupeR", quietly = TRUE))
    stop("Package 'loupeR' is required for the export_loupe step. ",
         "Install it from 10x Genomics and run loupeR::setup() once.")

  ow <- options(warn = 1); on.exit(options(ow), add = TRUE)
  obj <- seurat_obj
  present <- intersect(assays, names(obj@assays))
  if (length(present) == 0) {
    warning("None of the requested assays are present; nothing to export.")
    return(character(0))
  }

  for (a in present)
    for (lyr in SeuratObject::Layers(obj[[a]])) {
      if (!grepl("^counts", lyr)) next
      m <- SeuratObject::LayerData(obj, assay = a, layer = lyr)
      if (inherits(m, "IterableMatrix")) {
        if (verbose) cat(sprintf(" Materializing %s/%s into memory...\n", a, lyr))
        SeuratObject::LayerData(obj, assay = a, layer = lyr) <- as(m, "dgCMatrix")
      }
    }

  obj <- tryCatch(
    Seurat::DietSeurat(obj, assays = present, layers = "counts",
                       dimreducs = names(obj@reductions)),
    error = function(e) {
      if (verbose) cat(sprintf(" DietSeurat skipped (%s); continuing.\n", conditionMessage(e)))
      obj
    })

  bc <- colnames(obj)
  new_bc <- .loupe_reencode_barcodes(bc)
  if (is.null(new_bc)) {
    if (verbose) cat(" WARN: re-encoded barcodes are not unique; leaving barcodes unchanged.\n")
  } else if (!identical(new_bc, bc)) {
    obj <- tryCatch(
      SeuratObject::RenameCells(obj, new.names = new_bc),
      error = function(e) {
        if (verbose) cat(sprintf(" WARN: could not rename barcodes (%s); leaving unchanged.\n",
                                 conditionMessage(e)))
        obj
      })
    if (verbose)
      cat(sprintf(" Re-encoded %d barcodes to 10x format (<seq>-<gem group>, %d gem groups).\n",
                  length(bc), nlevels(factor(sub("^[^-]*", "", bc)))))
  }

  cluster_cols <- .loupe_cluster_columns(obj@meta.data, max_groups = 200L)
  for (cc in cluster_cols)
    if (!is.factor(obj@meta.data[[cc]]))
      obj@meta.data[[cc]] <- as.factor(as.character(obj@meta.data[[cc]]))
  if (verbose)
    cat(sprintf(" Using %d metadata column(s) as Loupe clusters (<= 200 groups each).\n",
                length(cluster_cols)))
  for (r in names(obj@reductions))
    if (ncol(SeuratObject::Embeddings(obj, r)) != 2) obj[[r]] <- NULL

  created <- character(0)
  for (a in present) {
    cnt <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = "counts"),
                    error = function(e) NULL)
    if (is.null(cnt)) { if (verbose) cat(sprintf(" SKIP %s: no counts layer.\n", a)); next }
    vals <- if (methods::is(cnt, "dgCMatrix")) cnt@x else as.numeric(cnt)
    finite <- vals[is.finite(vals)]
    is_int <- length(finite) == 0 || all(abs(finite - round(finite)) < 1e-8)
    if (!is_int) {
      if (verbose)
        cat(sprintf(paste0(" SKIP %s: Loupe requires integer counts, but %s counts are ",
                           "non-integer (min %.3f, max %.3f) - likely already normalized.\n"),
                    a, a, min(finite), max(finite)))
      next
    }

    SeuratObject::DefaultAssay(obj) <- a
    nm <- paste0(output_name, "_", a)
    if (verbose) cat(sprintf(" Writing Loupe file for %s assay: %s.cloupe\n", a, nm))
    ok <- tryCatch(
      loupeR::create_loupe_from_seurat(
        obj, output_dir = output_dir, output_name = nm,
        metadata_cols = cluster_cols, force = TRUE, executable_path = executable_path),
      error = function(e) { cat(sprintf(" ERROR exporting %s: %s\n", a, conditionMessage(e))); FALSE })
    if (isTRUE(ok)) {
      created <- c(created, file.path(output_dir, paste0(nm, ".cloupe")))
      if (verbose) cat(sprintf(" OK wrote %s.cloupe\n", nm))
    } else if (verbose) {
      cat(sprintf(" FAILED: loupeR did not write %s.cloupe (returned %s) - see warnings above.\n",
                  a, as.character(ok)))
    }
  }

  if (length(created) == 0)
    warning("export_seurat_to_loupe produced no .cloupe files.")
  created
}
