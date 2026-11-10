# Load Data Smart Function

## Context

The Load Data Smart Function is the entry point of the CITE-seq pipeline. It ingests one pre-assembled multi-modal object and produces a unified Seurat object ready for downstream processing.

The function performs one step:

**Data Loading** - Reads a single `.rds` (Seurat object), `.h5ad` (AnnData) or `.h5mu` (MuData) file; each already carries RNA and ADT together, which are separated automatically. For `.h5ad` files, Ensembl gene IDs are converted to gene symbols using `AnnotationDbi::mapIds()` with org.Hs.eg.db (human) or org.Mm.eg.db (mouse) - no Azimuth dependency. Gene names that are already symbols are passed through unchanged. For `.h5mu` files, the muon Python library is used via reticulate. Because `.h5ad`/`.h5mu` reading depends on Python imported in-process, the module probes that backend first and fails loudly (`H5AD_PYTHON_BACKEND_UNAVAILABLE`) rather than loading a CITE-seq file as RNA-only.

Integrating external sample metadata is **not** part of this module - it lives in the standalone [add_metadata](../add_metadata/README.md) smart function, which can run at any point after loading.

## Main Commands

Single command with **input** from `$SF_INPUT_PATH` and **output** to `$SF_OUTPUT_PATH`:

- **run** - Load the input object -> writes a Seurat `.rds` file to `$SF_OUTPUT_PATH/output/` (the collected output-port directory).

## Prerequisites

- R 4.5.1+ or Docker
- For `.h5mu` / `.h5ad` files: the Python stack bundled in the image (`muon`, `anndata`, `scanpy`). The interpreter is supplied by the image via `RETICULATE_PYTHON`; it is not a configuration setting.

## Environment (Docker / Smart Function)

| Variable | Purpose |
|--------------------|----------------------------------------------------------------------------|
| `SF_INPUT_PATH` | Directory for **inputs** (raw data files and `configs/` folder). |
| `SF_OUTPUT_PATH` | Directory for **outputs** (deliverables go to its `output/` port subdirectory). |
| `SF_STATUS_FILE` | If set, writes `SUCCESS` on completion (used by Genaiz platform). |
| `LOAD_DATA_CONFIG` | Config filename inside `$SF_INPUT_PATH/configs/` (default: `config.json`). |

## Commands

### run

Loads the input object (RNA + ADT are split automatically) and saves the result as a compressed Seurat RDS file.

```bash
Rscript main.R --config configs/config.json
```

- **Input:**
  - Pre-processed object: path specified in `load_data_path` (relative to `$SF_INPUT_PATH`)
  - Config: `$SF_INPUT_PATH/configs/$LOAD_DATA_CONFIG`
- **Output:** `$SF_OUTPUT_PATH/output/<output_name>_YYYYMMDD_HHMMSS.rds` (timestamp automatically appended for compliance; the `output/` port subdirectory is what the platform collects into the output dataset; the dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers - the two must always travel together)

## Supported Input Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| AnnData | `.h5ad` | Ensembl IDs converted to symbols; ADT split from `AB_` features or `obsm` |
| MuData | `.h5mu` | RNA + ADT loaded via muon (Python) |
| Seurat RDS | `.rds` | One file, must contain a Seurat object; bare matrices and directories are rejected |

Raw formats (`.csv`, `.csv.gz`, `.tar.gz`, 10X `.h5`) are rejected. This module loads an object that has already been pre-processed; RNA and ADT are separated automatically when both are present in the file.

## Config Files

- **`configs/config.example.json`** - Example configuration. Copy to `config.json` and fill in. All file paths are relative to `$SF_INPUT_PATH`.


## Error Reporting

Any failure writes **`error.json`** - to `$SF_OUTPUT_PATH` and mirrored into the collected `$SF_OUTPUT_PATH/output/` port directory, so it reaches the platform's files API - and exits with status 1, so the front-end can display a specific, user-readable message. A successful run writes **no** `error.json`; the `SF_STATUS_FILE` `SUCCESS` behaviour is unchanged.

| Code | Meaning |
|------|---------|
| `CONFIG_MISSING` | No `--config` given, or the config file does not exist |
| `CONFIG_PARSE_FAILED` | Config file is not valid JSON |
| `CONFIG_VALIDATION_FAILED` | Config keys are invalid, e.g. `load_data_path` missing or pointing nowhere |
| `INPUT_LOAD_FAILED` | The input file could not be read |
| `INPUT_NOT_SELF_CONTAINED` | The input `.rds` carries disk-backed matrix layers (e.g. BPCells) whose backing files live on the machine that created the object and are not part of the upload, so the matrix values are unreadable in any other container. Re-export with in-memory matrices, or provide the original `.h5ad`/`.h5mu` |
| `INPUT_NOT_SEURAT` | The `.rds` holds a bare matrix/list, not a pre-processed Seurat object |
| `UNSUPPORTED_FORMAT` | Input is not one of `.rds`, `.h5ad`, `.h5mu` (raw formats are rejected by design) |
| `H5AD_PYTHON_BACKEND_UNAVAILABLE` | For `.h5ad`/`.h5mu` input: `h5py`/`anndata` could not be imported in-process, so the file cannot be read completely. The run stops rather than silently dropping the ADT data |
| `H5AD_PARSE_FAILED` | The `.h5ad` could not be parsed into a Seurat object |
| `H5AD_SQUARE_MATRIX` | The `.h5ad` X matrix is square (cells == genes), which the reader cannot process safely |
| `H5MU_PARSE_FAILED` | The `.h5mu` could not be parsed (e.g. no ADT/PROTEIN modality found) |
| `UNSUPPORTED_SPECIES` | `species_name` is not `human`/`mouse` |
| `GENE_ID_CONVERSION_FAILED` | A required annotation package for gene ID conversion is missing |
| `OUTPUT_WRITE_FAILED` | The output `.rds` could not be written |
| `PIPELINE_FAILED` | Any other unexpected runtime error (backstop; raw R error in `details`) |

The helper lives in `error_codes.R`, vendored byte-identically into every smart function and compared by `make drift`.

## Integrations

### Docker

```bash
make docker

docker run --rm \
  -e SF_INPUT_PATH=/inputs \
  -e SF_OUTPUT_PATH=/outputs \
  -v /host/inputs:/inputs \
  -v /host/outputs:/outputs \
  cite-seq-load-data:latest
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
