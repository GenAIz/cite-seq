# Add Metadata Configuration Parameters

The following parameters are required to run the Add Metadata Smart Function within the CITE-seq Solution.

## Workspace Configuration Parameters

### Config File Name

* Type: String
* Description: Name of the configuration file to load from the `configs/` input port.
* Location: `METADATA_CONFIG` environment variable
* Configuration: Defaults to `config.json` if not set.

### Input Path

* Type: String
* Description: Directory from which input files (input object, metadata file, and configs) are read.
* Location: `SF_INPUT_PATH` environment variable
* Configuration: Set by the Genaiz platform at runtime. Defaults to `run/in` when running locally.

### Output Path

* Type: String
* Description: Directory to which the output files are written.
* Location: `output/` subdirectory of the `SF_OUTPUT_PATH` environment variable (the platform collects the output port dataset from there)
* Configuration: Set by the Genaiz platform at runtime. Defaults to `run/out` when running locally.

---

## Workflow Configuration Parameters

### input_rds (optional)

* Type: String
* Description: Path (relative to the input directory) to the `.rds` object to add metadata to.
* Note: when unset or missing, the newest `.rds` found anywhere under the `data` input port is used automatically.

### metadata_file

* Type: String
* Required: only when the `attachment` input port holds no metadata file.
* Description: Path (relative to the input directory) to the external metadata file (.xlsx, .xls, .csv, .tsv, .txt).
* Note: when the path does not exist, a file with the same name is looked up on the `attachment` input port; when unset, the newest metadata file on that port is used.
* Note: any other format, including `.rds`, fails validation.

### object_key (required)

* Type: String
* Description: Column name in the object's metadata to use for joining (e.g., "orig.ident", "MULTI_ID").
* Note: the former `seurat_key` spelling is still accepted.

### metadata_key (required)

* Type: String
* Description: Column name in the external metadata file to match with object_key.

### sheet

* Type: String or Integer
* Default: 1
* Description: For Excel files, the sheet name or number to read.

### output_name

* Type: String
* Default: "metadata.rds"
* Description: Base name for the output `.rds` file; a run timestamp is appended automatically.
* Note: the output dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers; the two must always travel together.

### verbose

* Type: Boolean
* Default: true
* Description: Print progress messages during execution.

## Advanced parameters

These are not in the example config. They keep their defaults unless you add them:

| Key | Default | Effect |
|---|---|---|
| `col_types` | `"text"` | Column types passed to `readxl::read_excel()`. Excel only. Allowed values (validated): `"text"`, `"guess"`, `"list"`. |
| `subset_to_existing` | `true` | Keep only cells whose key matches a metadata row. |
| `overwrite_existing` | `false` | On a column-name clash, `false` preserves the object's copy. |
| `add_match_status` | `true` | Adds a column recording whether each cell matched. |
| `backup_metadata` | `false` | Stashes the object's original metadata in its `misc` slot first. |
| `na_values` | `"", "NA", "N/A", "na", "n/a"` | Strings read as `NA`. |
