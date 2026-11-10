# Functional State Profile Smart Function

## Context

The Functional State Profile Smart Function scores T cells against a panel of biologically defined functional states (15 in the bundled RNA marker file, 8 surface states in the ADT file) using UCell, a rank-based gene signature scoring method. UCell computes a Mann-Whitney U statistic per cell per gene set, making scores robust to dataset size and sequencing depth differences.

Each cell receives a continuous score (0-1) for each functional state, and a `dominant_state` label is assigned based on the highest-scoring state (optionally `unassigned` below a configurable `min_score`; a `dominant_state_margin` column records top-minus-runner-up). An optional validation step checks score distributions and flags unexpected results.

The functional states scored by the bundled RNA marker file are:

| State | Biological Meaning |
|-------|--------------------|
| `stress-response` | Heat shock and cellular stress adaptation (HSP70/HSP90 families) |
| `anergy` | Hyporesponsive T cells with impaired TCR signalling (NR4A/EGR families) |
| `cell-death` | Apoptotic and cell death pathways |
| `naive-like` | Undifferentiated or recently activated cells |
| `memory` | Antigen-experienced compartment (TCM/TEM/TEMRA/TRM/TSCM) |
| `virtual-memory` | Memory-phenotype cells arising without cognate antigen |
| `adhesion` | Tissue residency and migratory adhesion programmes |
| `senescence` | Irreversible growth arrest; terminally differentiated |
| `exhaustion` | Chronic antigen exposure with progressive loss of effector function |
| `regulating` | Immunosuppressive / regulatory T cell programme |
| `cytotoxicity` | Effector killing capacity (granzymes, perforin) |
| `proliferation` | Actively cycling cells |
| `quiescence` | Resting, metabolically inactive cells |
| `IFN-response` | Interferon-stimulated gene expression |
| `activated` | Recently activated effector state |

Gene signatures for each state (core genes and additional genes) are defined in `cell_state_markers.json`, which ships with the image and is documented with literature references.

## Main Commands

Single command with **input** from `$SF_INPUT_PATH` and **output** to `$SF_OUTPUT_PATH`:

- **run** - Score cells against every functional state in the marker JSON and assign dominant state labels -> writes Seurat `.rds` and per-cell scores `.csv` to `$SF_OUTPUT_PATH/output/` (the collected output-port directory).

## Prerequisites

- R 4.5.1+ or Docker
- No Python dependencies

## Environment (Docker / Smart Function)

| Variable | Purpose |
|------------------|----------------------------------------------------------------------------------|
| `SF_INPUT_PATH` | Directory for **inputs** (annotated Seurat `.rds`, optional marker JSON, `configs/`). |
| `SF_OUTPUT_PATH` | Directory for **outputs** (deliverables go to its `output/` port subdirectory). |
| `SF_STATUS_FILE` | If set, writes `SUCCESS` on completion (used by Genaiz platform). |
| `FSP_CONFIG` | Config filename inside `$SF_INPUT_PATH/configs/` (default: `config.json`). |

## Commands

### run

Loads an annotated Seurat RDS, computes UCell scores for every functional state in the marker JSON, assigns dominant state labels, and saves both the enriched Seurat object and a per-cell score CSV.

```bash
Rscript main.R --config configs/config.json
```

- **Input:**
  - Seurat RDS: path specified in `input_rds` (relative to `$SF_INPUT_PATH`); when unset or not found, the newest `.rds` on the `data` input port (`$SF_INPUT_PATH/data/`) is used, so chained runs with timestamped upstream filenames work
  - Marker JSON: specified in `marker_json_path`; resolved from the input directory, then from the `attachment` input port (same name first, newest `.json` otherwise), with the bundled `cell_state_markers.json` as the final fallback
  - Config: `$SF_INPUT_PATH/configs/$FSP_CONFIG`
- **Output:**
  - `$SF_OUTPUT_PATH/output/<output_rds>_YYYYMMDD_HHMMSS.rds` - Seurat object with UCell score columns and `dominant_state`; the dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers - the two must always travel together
  - `$SF_OUTPUT_PATH/output/<output_csv>_YYYYMMDD_HHMMSS.csv` - Per-cell scores for all states plus `dominant_state`, `dominant_state_score`, `dominant_state_margin`

## Output Columns Added to Seurat Metadata

| Column | Type | Description |
|--------|------|-------------|
| `dominant_state` | String | Name of the highest-scoring functional state for each cell (`unassigned` when `min_score` is set and no state reaches it) |
| `dominant_state_score` | Float | UCell score of the dominant state |
| `dominant_state_margin` | Float | Top score minus runner-up score (small = near-tie) |
| `stress-response` | Float | UCell score (0-1) |
| `anergy` | Float | UCell score (0-1) |
| `cell-death` | Float | UCell score (0-1) |
| `naive-like` | Float | UCell score (0-1) |
| `memory` | Float | UCell score (0-1) |
| `virtual-memory` | Float | UCell score (0-1) |
| `adhesion` | Float | UCell score (0-1) |
| `senescence` | Float | UCell score (0-1) |
| `exhaustion` | Float | UCell score (0-1) |
| `regulating` | Float | UCell score (0-1) |
| `cytotoxicity` | Float | UCell score (0-1) |
| `proliferation` | Float | UCell score (0-1) |
| `quiescence` | Float | UCell score (0-1) |
| `IFN-response` | Float | UCell score (0-1) |
| `activated` | Float | UCell score (0-1) |


## Error Reporting

Any failure writes **`error.json`** - to `$SF_OUTPUT_PATH` and mirrored into the collected `$SF_OUTPUT_PATH/output/` port directory, so it reaches the platform's files API - and exits with status 1, so the front-end can display a specific, user-readable message. A successful run writes **no** `error.json` (presence of the file = failure); the `SF_STATUS_FILE` `SUCCESS` behaviour is unchanged.

```json
{
  "status": "error",
  "smart_function": "functional_state_profile",
  "error_code": "MARKER_STATE_EMPTY",
  "message": "One or more functional states have no marker genes present in the input data. ...",
  "details": ["No marker genes survive the gene universe for state(s): memory. ..."],
  "step": "classification",
  "timestamp": "2026-08-10T21:40:12-0600"
}
```

| Code | Meaning |
|------|---------|
| `CONFIG_MISSING` | No `--config` given, or the config file does not exist |
| `CONFIG_PARSE_FAILED` | Config file is not valid JSON |
| `CONFIG_VALIDATION_FAILED` | Config keys are invalid (details list every problem) |
| `INPUT_LOAD_FAILED` | The input `.rds` could not be read |
| `INPUT_NOT_SELF_CONTAINED` | The input `.rds` carries disk-backed matrix layers (e.g. BPCells) whose backing files live on the machine that created the object and are not part of the upload, so the matrix values are unreadable in any other container. Re-export with in-memory matrices, or provide the original `.h5ad`/`.h5mu` |
| `INPUT_NOT_SEURAT` | The input `.rds` does not contain a Seurat object |
| `MARKER_JSON_MISSING` | Marker definitions JSON not found |
| `MARKER_STATE_EMPTY` | A state has no marker genes present in the data (often Ensembl-ID input or over-restrictive HVG filter) |
| `ASSAY_NOT_FOUND` | The configured assay is not in the input object |
| `SCORING_FAILED` | UCell scoring failed (missing package or scoring error) |
| `OUTPUT_WRITE_FAILED` | Output `.rds`/`.csv` could not be written |
| `PIPELINE_FAILED` | Any other unexpected runtime error (backstop; raw R error in `details`) |

The helper lives in `error_codes.R`, which is vendored byte-identically into every smart function and compared by `make drift`.

## Config Files

- **`configs/config.example.json`** - Example configuration. Copy to `config.json` and fill in.
- **`cell_state_markers.json`** - RNA gene signatures for all 15 states with literature references and per-state `curation_notes`. Bundled in the image. Can be overridden by placing a custom JSON at `$SF_INPUT_PATH/<marker_json_path>`.
- **`cell_state_markers_adt.json`** - Surface (ADT) marker signatures in canonical protein names with an `_aliases` map; run with `assay = "ADT"`. Vendor-decorated panel names are resolved by whole-token matching.

## Integrations

### Docker

```bash
make docker

docker run --rm \
  -e SF_INPUT_PATH=/inputs \
  -e SF_OUTPUT_PATH=/outputs \
  -v /host/inputs:/inputs \
  -v /host/outputs:/outputs \
  cite-seq-functional-state-profile:latest
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
