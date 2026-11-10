suppressPackageStartupMessages(library(reticulate))

hdr <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")
ok <- function(x) cat(" [ok] ", x, "\n")
bad <- function(x) cat(" [FAIL] ", x, "\n")

hdr("1. Interpreter and package versions")
cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
cat(" reticulate python :", if (is.null(cfg)) "<none>" else cfg$python, "\n")
cat(" R anndata pkg :", requireNamespace("anndata", quietly = TRUE), "\n")
if (requireNamespace("anndata", quietly = TRUE))
  cat(" R anndata version :", as.character(utils::packageVersion("anndata")), "\n")
cat(" BPCells version :", as.character(utils::packageVersion("BPCells")), "\n")
py_run_string("
import anndata, h5py, scipy, numpy
print(' python anndata :', anndata.__version__)
print(' python h5py :', h5py.__version__)
print(' python scipy :', scipy.__version__)
")

hdr("2. Write the fixture exactly as helper-fixtures.R does")
tmp <- file.path(tempdir(), "fixture.h5ad")
py_run_string(sprintf("
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
cells = ['Cell%%d' %% i for i in range(1, 11)]
genes = ['Gene%%d' %% i for i in range(1, 6)]
X = sp.csr_matrix(np.random.poisson(5, (10, 5)).astype('float32'))
obs = pd.DataFrame({'batch': ['A']*5 + ['B']*5}, index=cells)
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
a.write_h5ad(r'%s')
print(' wrote', r'%s')
", tmp, tmp))

dump_py <- "
import h5py
def dump(path, label):
    print('\\n --- %s ---' % label)
    with h5py.File(path, 'r') as f:
        X = f['X']
        print(' X :', type(X).__name__, dict(X.attrs) if X.attrs else {})
        if isinstance(X, h5py.Group):
            for k in X: print(' X/%-8s dtype=%s shape=%s' % (k, X[k].dtype, X[k].shape))
        o = f['obs']
        idx = o.attrs.get('_index', b'_index')
        idx = idx.decode() if isinstance(idx, bytes) else idx
        print(' obs._index :', repr(idx))
        d = o[idx]
        if isinstance(d, h5py.Group):
            # anndata >= 0.12 encodes string arrays as a GROUP. BPCells opens
            # obs/_index as a DATASET, hence '(Dataset) Inappropriate type'.
            print(' ->', idx, ': GROUP <-- BPCells cannot read this',
                  dict(d.attrs))
            for k in d:
                print(' %s/%s dtype=%s shape=%s' % (idx, k, d[k].dtype, d[k].shape))
        else:
            t = d.id.get_type()
            print(' ->', idx, ': DATASET vlen=%s cset=%s size=%s' %
                  (t.is_variable_str(), t.get_cset(), t.get_size()))
        v = f['var']
        vidx = v.attrs.get('_index', b'_index')
        vidx = vidx.decode() if isinstance(vidx, bytes) else vidx
        vd = v[vidx]
        print(' var._index :', repr(vidx),
              'GROUP' if isinstance(vd, h5py.Group) else
              'DATASET vlen=%s' % vd.id.get_type().is_variable_str())
"
py_run_string(dump_py)

REAL <- "/raw/Yoshida_et_al_2022_airway.h5ad"
hdr("3. On-disk structure: fixture vs real file")
py_run_string(sprintf("dump(r'%s', 'FIXTURE (written here)')", tmp))
if (file.exists(REAL)) {
  py_run_string(sprintf("dump(r'%s', 'REAL Yoshida h5ad')", REAL))
} else {
  cat("\n (real file not mounted at ", REAL, " - skipping)\n", sep = "")
}

hdr("4. BPCells::open_matrix_anndata_hdf5")
try_bp <- function(path, label) {
  if (!file.exists(path)) { cat(" --", label, ": not mounted, skipped\n"); return(invisible()) }
  res <- tryCatch({
    m <- BPCells::open_matrix_anndata_hdf5(path = path)
    ok(sprintf("%-26s dims = %d x %d", label, nrow(m), ncol(m)))
    TRUE
  }, error = function(e) {
    bad(sprintf("%-26s %s", label, conditionMessage(e))); FALSE
  })
  invisible(res)
}
try_bp(tmp, "FIXTURE")
try_bp(REAL, "REAL Yoshida")

hdr("5. .read_h5ad_metadata() - verbose, both backends")
src <- if (file.exists("/app/script/load_data_updated.R"))
  "/app/script/load_data_updated.R" else "script/load_data_updated.R"
cat(" sourcing", src, "\n")
suppressWarnings(suppressPackageStartupMessages(source(src)))

cat("\n -- with the R anndata package available (as installed) --\n")
r1 <- tryCatch(.read_h5ad_metadata(tmp, verbose = TRUE),
               error = function(e) { bad(conditionMessage(e)); list(obs = NULL) })
cat(" obs recovered:", !is.null(r1$obs),
    if (!is.null(r1$obs)) sprintf(" (%d x %d, cols: %s)", nrow(r1$obs), ncol(r1$obs),
                                  paste(colnames(r1$obs), collapse = ",")) else "", "\n")

hdr("6. Isolate the two backends")
cat(" -- R anndata::read_h5ad directly --\n")
if (requireNamespace("anndata", quietly = TRUE)) {
  tryCatch({
    a <- anndata::read_h5ad(tmp)
    cat(" class(a$obs) :", paste(class(a$obs), collapse = "/"), "\n")
    cat(" nrow(a$obs) :", paste(nrow(a$obs), collapse = ","), "\n")
    cat(" is.data.frame:", is.data.frame(a$obs), "\n")
  }, error = function(e) bad(paste("anndata::read_h5ad ->", conditionMessage(e))))
} else cat(" R anndata not installed\n")

cat("\n -- reticulate + py_to_r directly --\n")
tryCatch({
  adp <- reticulate::import("anndata", convert = FALSE)
  a2 <- adp$read_h5ad(tmp)
  df <- reticulate::py_to_r(a2$obs)
  cat(" class(df) :", paste(class(df), collapse = "/"), "\n")
  cat(" dim(df) :", paste(dim(df), collapse = " x "), "\n")
  cat(" colnames(df) :", paste(colnames(df), collapse = ","), "\n")
}, error = function(e) bad(paste("reticulate ->", conditionMessage(e))))

hdr("done")
