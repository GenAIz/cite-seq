# Add Metadata

Integrates external sample/patient metadata from Excel, CSV, or TSV files into a Seurat object by performing a left join on a specified key column.

This module can be run at **any point** in the pipeline after data loading - it only requires a Seurat `.rds` file as input and a metadata file to join.

## Quick Start

```bash
# Local
Rscript main.R --config configs/config.example.json

# Docker
docker build -t cite-seq-add-metadata .
docker run -v $(pwd)/run/in:/app/run/in -v $(pwd)/run/out:/app/run/out \
  cite-seq-add-metadata

# Tests
Rscript tests/testthat.R
```

## Supported Metadata Formats

- Excel (.xlsx, .xls)
- CSV (.csv)
- TSV / tab-delimited (.tsv, .txt)

## Key Features

- Left join preserves all cells in the Seurat object
- Cells without matching metadata receive NA values
- Optional match status column for QC
- Optional backup of original metadata
- Join statistics stored in `@misc$metadata_join_stats`
- Compatible with Seurat v4 and v5
- When `input_rds` is unset or the configured path is missing, the newest `.rds` on the `data` input port (`$SF_INPUT_PATH/data/`) is used automatically, so chained pipeline runs work with timestamped upstream filenames
- Metadata files uploaded to the `attachment` input port are discovered automatically: a configured name whose path does not exist is matched there, and with `metadata_file` unset the newest supported file on the port is used
- Outputs are written to `$SF_OUTPUT_PATH/output/`, the port subdirectory the platform collects; the dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers - the two must always travel together

## Configuration

See [config.md](config.md) for full parameter documentation.

## Error Reporting

Any failure writes **`error.json`** - to `$SF_OUTPUT_PATH` and mirrored into the collected `$SF_OUTPUT_PATH/output/` port directory, so it reaches the platform's files API - and exits with status 1, so the front-end can display a specific, user-readable message. A successful run writes **no** `error.json`; the `SF_STATUS_FILE` `SUCCESS` behaviour is unchanged.

| Code | Meaning |
|------|---------|
| `CONFIG_MISSING` | No `--config` given, or the config file does not exist |
| `CONFIG_PARSE_FAILED` | Config file is not valid JSON |
| `CONFIG_VALIDATION_FAILED` | Required keys missing/invalid (details list every problem) |
| `INPUT_LOAD_FAILED` | The input Seurat `.rds` could not be read |
| `INPUT_NOT_SELF_CONTAINED` | The input `.rds` carries disk-backed matrix layers (e.g. BPCells) whose backing files live on the machine that created the object and are not part of the upload, so the matrix values are unreadable in any other container. Re-export with in-memory matrices, or provide the original `.h5ad`/`.h5mu` |
| `METADATA_FILE_MISSING` | The metadata file does not exist |
| `METADATA_LOAD_FAILED` | The metadata file exists but could not be read/parsed |
| `UNSUPPORTED_METADATA_FORMAT` | Metadata is not xlsx/xls/csv/tsv/txt |
| `EXCEL_SUPPORT_MISSING` | Excel input but readxl is not installed in the image |
| `OBJECT_KEY_NOT_FOUND` | `object_key` column absent from the object metadata |
| `METADATA_KEY_NOT_FOUND` | `metadata_key` column absent from the metadata file |
| `NO_METADATA_MATCHES` | No cells matched any metadata row (key mismatch) |
| `METADATA_JOIN_CORRUPT` | Join would change the cell count (non-unique key); run stopped to protect the object |
| `OUTPUT_WRITE_FAILED` | The output `.rds` could not be written |
| `PIPELINE_FAILED` | Any other unexpected runtime error (backstop; raw R error in `details`) |

The helper lives in `error_codes.R`, vendored byte-identically into every smart function and compared by `make drift`.
