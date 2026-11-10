# Cell Type Annotation Smart Function

## Context

The Cell Type Annotation Smart Function assigns cell type labels to a processed Seurat object. It supports two complementary annotation methods - SingleR and STCAT - which can be run independently or together.

The function performs one step with two optional sub-methods:

1. **SingleR** - A reference-based method that uses bulk RNA-seq reference datasets (Monaco, HPCA, Blueprint/ENCODE, DICE, or ImmGen) to label individual cells. Labels are assigned at both broad (`main`) and fine-grained (`fine`) levels. Uses a **pseudo-bulk aggregation** strategy: expression is averaged per cluster (n_genes x n_clusters dense matrix) instead of converting the full n_genes x n_cells matrix to SingleCellExperiment, keeping peak RAM proportional to n_clusters rather than n_cells. This makes SingleR tractable on datasets with >500k cells. Runs with `SerialParam` to avoid fork-based OOM inside Docker.

2. **STCAT** - A Python-based deep learning model for immune cell type annotation. Invoked as a subprocess via `python/run_stcat.py` using the Python interpreter bundled in the image. The R side streams the expression matrix from disk to h5ad via `BPCells::write_matrix_anndata_hdf5()` (no `dgCMatrix` materialisation in R RAM). The Python script monkey-patches `multiprocessing.Process` to run STCAT's parallel CD4/CD8 subtyping sequentially, eliminating the ~3x fork-copy RAM spike. scipy sparse dtype compatibility (`uint32` -> `int32`) is corrected automatically. Produces a `Prediction` column in the Seurat metadata.

When `method = "both"`, the two methods run as one coordinated sequence rather than side by side. This module and the vendored copy in `alternate_workflow` are kept byte-identical apart from their header comment, so a smart-function run and a full-pipeline run produce the same annotation:

1. **SingleR first**, at cluster level (main + fine).
2. **Per-cell fine refinement** - fine labels are re-annotated cell by cell in bounded chunks (`singler_fine_chunk_size`), so a mixed cluster is no longer blanketed with a single label. Memory scales with the chunk, not the dataset.
3. **Vocabulary rescue** - platelets and erythrocytes, which the Monaco reference has no labels for (SingleR files them under monocytes), are relabelled from canonical markers when the evidence is decisive. Recorded in `vocab_rescue`.
4. **T cell gate** - only cells whose `SingleR_main` **or** `SingleR_fine` label references T cells are sent to STCAT; the rest keep their SingleR label. STCAT is a T cell tool, so this both cuts its runtime roughly in half and stops it being asked about cells it cannot classify. STCAT's neighbour graph is built from a T-subset re-embedding (subset HVGs -> PCA -> Harmony when `stcat_batch_column` has more than one level), so its sub-clustering resolves T substructure instead of separating by batch.
5. **Label harmonization** - Monaco fallback labels that duplicate STCAT concepts (`MAIT cells` -> `MAIT`, `Naive CD4 T cells` -> `CD4 Tn`, ...) are renamed to the STCAT vocabulary, so one cell type has one name. `SingleR_fine` itself is left reference-faithful.

## Main Commands

Single command with **input** from `$SF_INPUT_PATH` and **output** to `$SF_OUTPUT_PATH`:

- **run** - Annotate cell types in a processed Seurat object -> writes an annotated Seurat `.rds` to `$SF_OUTPUT_PATH/output/` (the collected output-port directory).

## Prerequisites

- R 4.5.1+ or Docker
- For STCAT: the Python stack bundled in the image (`STCAT==1.0.8`, `scanpy`, `anndata`, `scikit-learn`, `umap-learn`). The interpreter is supplied by the image via `RETICULATE_PYTHON`; it is not a configuration setting.

## Environment (Docker / Smart Function)

| Variable | Purpose |
|--------------------|----------------------------------------------------------------------------------|
| `SF_INPUT_PATH` | Directory for **inputs** (Seurat `.rds` and `configs/` folder). |
| `SF_OUTPUT_PATH` | Directory for **outputs** (deliverables go to its `output/` port subdirectory). |
| `SF_STATUS_FILE` | If set, writes `SUCCESS` on completion (used by Genaiz platform). |
| `ANNOTATION_CONFIG`| Config filename inside `$SF_INPUT_PATH/configs/` (default: `config.json`). |

## Commands

### run

Loads a processed Seurat RDS, runs cell type annotation using SingleR and/or STCAT, and saves the annotated object.

```bash
Rscript main.R --config configs/config.json
```

- **Input:**
  - Seurat RDS: path specified in `input_rds` (relative to `$SF_INPUT_PATH`); when unset or not found, the newest `.rds` on the `data` input port (`$SF_INPUT_PATH/data/`) is used, so chained runs with timestamped upstream filenames work
  - Config: `$SF_INPUT_PATH/configs/$ANNOTATION_CONFIG`
- **Output:** `$SF_OUTPUT_PATH/output/<output_name>_YYYYMMDD_HHMMSS.rds` (timestamp automatically appended for compliance; the `output/` port subdirectory is what the platform collects into the output dataset; the dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers - the two must always travel together)

## Output Columns Added to Seurat Metadata

| Column | Source | Description |
|--------|--------|-------------|
| `SingleR_main` | SingleR | Broad cell type label (cluster level) |
| `SingleR_main_score` | SingleR | Confidence of the main label |
| `SingleR_fine` | SingleR | Fine-grained cell type label, assigned **per cell**. Reference-faithful Monaco names (e.g. `Vd2 gd T cells`), never harmonized - the place to look when `Prediction` collapses a distinction |
| `SingleR_fine_score` | SingleR | Confidence of the fine label |
| `vocab_rescue` | Markers | `Platelets` / `Erythrocytes` for rescued cells, `NA` otherwise. Also the flag for platelet-monocyte aggregates, which some analyses prefer to exclude |
| `Prediction` | STCAT + SingleR | Final cell type: STCAT's call for gated T cells, the fallback label otherwise, harmonized to one vocabulary |
| `Uncertainty_score` | STCAT | STCAT uncertainty **and** a provenance flag: non-`NA` means the label came from STCAT, `NA` means it came from SingleR |

Working columns used during annotation (`SingleR_*_delta`, `SingleR_fine_cluster*`, `Cluster`, `vocab_rescue_score`, `seurat_clusters`) are dropped before the object is written.

## Config Files

- **`configs/config.example.json`** - Example configuration. Copy to `config.json` and fill in. `input_rds` is relative to `$SF_INPUT_PATH`.

## Integrations

### Docker

```bash
make docker

docker run --rm \
  -e SF_INPUT_PATH=/inputs \
  -e SF_OUTPUT_PATH=/outputs \
  -v /host/inputs:/inputs \
  -v /host/outputs:/outputs \
  cite-seq-cell-type-annotation:latest
```

### Genaiz

```bash
make genaiz
```

Test runs can be made using:

```bash
make genaiz-test
```

The above uses the configuration found at `run/in/configs/config.json`, which can be set up using `configs/config.example.json` as a reference.

## Error Reporting

Any failure writes **`error.json`** - to `$SF_OUTPUT_PATH` and mirrored into the collected `$SF_OUTPUT_PATH/output/` port directory, so it reaches the platform's files API - and exits with status 1, so the front-end can display a specific, user-readable message. A successful run writes **no** `error.json`; the `SF_STATUS_FILE` `SUCCESS` behaviour is unchanged.

| Code | Meaning |
|------|---------|
| `CONFIG_MISSING` | No `--config` given, or the config file does not exist |
| `CONFIG_PARSE_FAILED` | Config file is not valid JSON |
| `CONFIG_VALIDATION_FAILED` | Config keys invalid - note `method` accepts SingleR/STCAT/both (any case) and `reference` accepts monaco/hpca/blueprint_encode/dice/immgen |
| `INPUT_LOAD_FAILED` | The input `.rds` could not be read |
| `INPUT_NOT_SELF_CONTAINED` | The input `.rds` carries disk-backed matrix layers (e.g. BPCells) whose backing files live on the machine that created the object and are not part of the upload, so the matrix values are unreadable in any other container. Re-export with in-memory matrices, or provide the original `.h5ad`/`.h5mu` |
| `INPUT_NOT_SEURAT` | The input `.rds` does not contain a Seurat object |
| `ASSAY_NOT_FOUND` | The configured assay is not in the input object |
| `CLUSTER_COLUMN_INVALID` | The cluster column is missing or contains NA values (details list candidates) |
| `UNKNOWN_REFERENCE` | The reference dataset name is not recognised |
| `REFERENCE_LOAD_FAILED` | The celldex reference could not be loaded. All five references are baked into the image as an offline cache at build time (platform containers have no internet), so this usually means the image predates the baked cache |
| `NO_COMMON_GENES` | Input and reference share no gene names (often Ensembl-ID input) |
| `SINGLER_LABELS_MISSING` | STCAT-only run on an object with no `SingleR_main`/`SingleR_fine`: the T cell gate has nothing to gate on - run with `method` `both` |
| `ANNOTATION_FAILED` | Pseudo-bulk aggregation or SingleR itself failed |
| `PYTHON_UNAVAILABLE` | The Python environment for STCAT could not be configured |
| `STCAT_SCRIPT_MISSING` | The STCAT Python script was not found |
| `STCAT_FAILED` | STCAT ran but failed (details carry the Python output) |
| `OUTPUT_WRITE_FAILED` | The output `.rds` could not be written |
| `PIPELINE_FAILED` | Any other unexpected runtime error (backstop; raw R error in `details`) |

The helper lives in `error_codes.R`, vendored byte-identically into every smart function and compared by `make drift`.
