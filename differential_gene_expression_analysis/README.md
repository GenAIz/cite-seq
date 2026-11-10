# Differential Gene Expression Analysis Smart Function

## Context

The Differential Gene Expression Analysis Smart Function identifies genes that are significantly up- or down-regulated between biological conditions in an annotated CITE-seq dataset. It supports both single-cell native and pseudobulk testing strategies, allowing the method to be matched to the experimental design.

The function supports two broad testing approaches:

1. **Single-cell native testing** (`pseudobulking = false`) - Tests are run directly on individual cells using Seurat's `FindMarkers` framework. Supported methods are MAST, wilcox, bimod, roc, t, negbinom, poisson and LR. Each cluster is tested against the rest by default; `ident_1`/`ident_2` run a single two-group comparison instead, and `ident_1` alone compares that group against all other cells. Suitable for exploratory analysis and datasets without strong sample-level batch effects.

2. **Pseudobulk testing** (`pseudobulking = true`) - Aggregates counts per sample per cluster before testing, then applies bulk RNA-seq methods (DESeq2, limma, or edgeR). This approach is statistically more appropriate for multi-sample studies and reduces inflation of false positives caused by treating cells as independent replicates.

Results are returned as a flat data frame with gene-level statistics (log fold change, p-value, adjusted p-value) across all cell type clusters.

## Main Commands

Single command with **input** from `$SF_INPUT_PATH` and **output** to `$SF_OUTPUT_PATH`:

- **run** - Perform DGE analysis on an annotated Seurat object -> writes a timestamped results `.csv` to `$SF_OUTPUT_PATH/output/` (the collected output-port directory).

The DGE *results* are a table, but so that this node can sit mid-chain, it also **passes its input Seurat object through unchanged**: every run writes a byte-exact copy of the input `.rds` (`dge_passthrough_YYYYMMDD_HHMMSS.rds`), plus its `_matrices/` directory when present, alongside the results CSV. A node chained downstream of DGE therefore receives exactly the object this analysis ran on - the DGE results themselves travel only in the CSV.

## Prerequisites

- R 4.5.1+ or Docker
- No Python dependencies

## Environment (Docker / Smart Function)

| Variable | Purpose |
|----------------|--------------------------------------------------------------------------------|
| `SF_INPUT_PATH` | Directory for **inputs** (annotated Seurat `.rds` and `configs/` folder). |
| `SF_OUTPUT_PATH`| Directory for **outputs** (deliverables go to its `output/` port subdirectory). |
| `SF_STATUS_FILE`| If set, writes `SUCCESS` on completion (used by Genaiz platform). |
| `DGE_CONFIG` | Config filename inside `$SF_INPUT_PATH/configs/` (default: `config.json`). |

## Commands

### run

Loads an annotated Seurat RDS, runs differential gene expression analysis, and saves the results table.

```bash
Rscript main.R --config configs/config.json
```

- **Input:**
  - Seurat RDS: path specified in `input_rds` (relative to `$SF_INPUT_PATH`); when unset or not found, the newest `.rds` on the `data` input port (`$SF_INPUT_PATH/data/`) is used, so chained runs with timestamped upstream filenames work
  - Config: `$SF_INPUT_PATH/configs/$DGE_CONFIG`
- **Output:** `$SF_OUTPUT_PATH/output/<output_name>_YYYYMMDD_HHMMSS.csv` - the results table as a flat CSV (timestamp automatically appended for compliance) - plus `$SF_OUTPUT_PATH/output/dge_passthrough_YYYYMMDD_HHMMSS.rds`, a byte-exact copy of the input object (with its `_matrices/` directory when present) so downstream nodes can chain from this one. The retired `output_rds`/`output_csv` keys are rejected at validation.

## Supported Test Methods

| Method | Mode | Package | Best For |
|--------|------|---------|----------|
| `MAST` | Single-cell | MAST | Exploratory analysis, small sample sizes |
| `wilcox` | Single-cell | Seurat | Fast, non-parametric comparisons |
| `bimod` | Single-cell | Seurat | Likelihood-ratio test for bimodal expression |
| `roc` | Single-cell | Seurat | Classifier power per gene (AUC) |
| `t` | Single-cell | Seurat | Student's t-test on expression |
| `negbinom` | Single-cell | Seurat | Negative binomial GLM on UMI counts |
| `poisson` | Single-cell | Seurat | Poisson GLM on UMI counts |
| `LR` | Single-cell | Seurat | Logistic regression, supports latent variables |
| `DESeq2` | Pseudobulk | DESeq2 | Multi-sample studies (recommended) |
| `limma` | Pseudobulk | limma/voom | Continuous covariates, complex designs |
| `edgeR` | Pseudobulk | edgeR | Small sample sizes with pseudobulking |

## Config Files

- **`configs/config.example.json`** - Example configuration including a pseudobulk DESeq2 example. Copy to `config.json` and fill in.

## Integrations

### Docker

```bash
make docker

docker run --rm \
  -e SF_INPUT_PATH=/inputs \
  -e SF_OUTPUT_PATH=/outputs \
  -v /host/inputs:/inputs \
  -v /host/outputs:/outputs \
  cite-seq-dge:latest
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
| `CONFIG_VALIDATION_FAILED` | Config keys invalid - `test_use` accepts MAST/wilcox/bimod/roc/t/negbinom/poisson/LR/DESeq2/limma/edgeR (any case), `design_strategy` accepts stratified/interaction, `pseudobulk_mode` accepts cluster_sample/sample_only/celltype_condition; `cluster_col` is always required, `sample_col`/`condition_col` are required when pseudobulking, and `ident_2` without `ident_1` is rejected; retired `output_rds`/`output_csv` keys are rejected here, before any compute |
| `INPUT_LOAD_FAILED` | The input `.rds` could not be read |
| `INPUT_NOT_SELF_CONTAINED` | The input `.rds` carries disk-backed matrix layers (e.g. BPCells) whose backing files live on the machine that created the object and are not part of the upload, so the matrix values are unreadable in any other container. Re-export with in-memory matrices, or provide the original `.h5ad`/`.h5mu` |
| `COLUMN_NOT_FOUND` | A configured metadata column (cluster/sample/condition/timepoint) is missing from the object |
| `INVALID_CONDITIONS` | `ident.1`/`ident.2`/`conditions_to_analyze`/`reference_condition` name values absent from the condition column |
| `INVALID_TIMEPOINTS` | `timepoints_to_analyze` names values absent from the timepoint column |
| `NO_DGE_RESULTS` | The analysis ran but every cluster/contrast was skipped (e.g. group sizes below `min_samples_per_group`) - surfaced as an error instead of a silently empty CSV |
| `OUTPUT_WRITE_FAILED` | The results CSV could not be written |
| `PIPELINE_FAILED` | Any other unexpected runtime error (backstop; raw R error in `details`) |

The helper lives in `error_codes.R`, vendored byte-identically into every smart function and compared by `make drift`.
