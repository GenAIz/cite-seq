# Functional State Profile Configuration Parameters

The following parameters are required to run the Functional State Profile Smart Function within the CITE-seq Solution.

## Workspace Configuration Parameters

The following parameters are tied to the Genaiz workspace.

### Config File Name

* Type: String
* Description: Name of the configuration file to load from the `configs/` input port.
* Location: `FSP_CONFIG` environment variable
* Configuration:
  + Workspace entry: Manual entry at creation.
  + Defaults to `config.json` if not set.

### Input Path

* Type: String
* Description: Directory from which input files (input object, optional marker JSON, and configs) are read.
* Location: `SF_INPUT_PATH` environment variable
* Configuration:
  + Set by the Genaiz platform at runtime.
  + Defaults to `run/in` when running locally.
  + Non-mutable at runtime.

### Output Path

* Type: String
* Description: Directory to which the output files are written.
* Location: `output/` subdirectory of the `SF_OUTPUT_PATH` environment variable (the platform collects the output port dataset from there)
* Configuration:
  + Set by the Genaiz platform at runtime.
  + Defaults to `run/out` when running locally.
  + Non-mutable at runtime.

---

## Workflow Configuration Parameters

The following parameters are defined inside `configs/config.json` and control execution of the Functional State Profile Smart Function.

### Input RDS

* Required: No
* Type: String
* Description: Path (relative to the input directory) to the annotated `.rds` object.
* Note: when unset or missing, the newest `.rds` found anywhere under the `data` input port is used automatically.
* Location: `input_rds` inside `configs/config.json`
* Example: `"data/annotated_object.rds"`

### Output RDS

* Required: Yes
* Type: String
* Description: Filename for the output `.rds` object with the UCell score columns added; a run timestamp is appended automatically.
* Note: the output dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers; the two must always travel together.
* Location: `output_rds` inside `configs/config.json`
* Default value: `"functional_state_object.rds"`

### Output CSV

* Required: Yes
* Type: String
* Description: Filename for the per-cell scores CSV (one row per cell, one column per state, plus the dominant-state columns).
* Location: `output_csv` inside `configs/config.json`
* Default value: `"functional_state_scores.csv"`

### Marker JSON Path

* Required: No
* Type: String
* Description: Path (relative to the input directory) to the gene signature definitions JSON.
* Note: when the path does not exist, the `attachment` input port is searched - first for the same name, then for the newest `.json`; when unset, the newest `.json` on that port is used; the bundled `cell_state_markers.json` inside the image is the final fallback.
* Location: `marker_json_path` inside `configs/config.json`
* Configuration:
  + To use custom gene signatures, provide a JSON following the same schema as the bundled file (keys: `core`, `additional`, `description`, `references`).
  + Default value: `"cell_state_markers.json"` (bundled with the image)

### Use Core Only

* Required: No
* Type: Boolean
* Description: If `true`, only each signature's `core` genes are scored, excluding `additional` genes.
* Location: `use_core_only` inside `configs/config.json`
* Default value: `false`

### Use Highly Variable Genes

* Required: No
* Type: Boolean
* Description: If `true`, restricts scoring to the object's highly variable genes.
* Location: `use_hvg` inside `configs/config.json`
* Default value: `false`

### Assay

* Required: No
* Type: String
* Description: Assay from which expression values are read for scoring.
* Location: `assay` inside `configs/config.json`
* Default value: `"RNA"`

### Max Rank

* Required: No
* Type: Integer
* Description: Maximum gene rank considered by the UCell statistic; genes ranked lower do not contribute.
* Location: `maxRank` inside `configs/config.json`
* Default value: `1500`

### Minimum Score

* Required: No
* Type: Number (0-1)
* Description: Minimum score the winning state must reach; cells below it are labelled `unassigned`.
* Note: when `null`, every cell gets its top-scoring state; the score and margin columns are written either way.
* Location: `min_score` inside `configs/config.json`
* Default value: unset (`null`) - no thresholding

### Number of Cores

* Required: No
* Type: Integer
* Description: Number of CPU cores to use for parallel UCell scoring via BiocParallel.
* Location: `ncores` inside `configs/config.json`
* Default value: `1`

### Run Validation

* Required: Yes
* Type: Boolean
* Description: If `true`, runs score-distribution checks after scoring and flags unexpected results.
* Location: `run_validation` inside `configs/config.json`
* Default value: `true`

### Verbose

* Required: No
* Type: Boolean
* Description: If `true`, prints detailed progress messages during scoring.
* Location: `verbose` inside `configs/config.json`
* Default value: `true`
