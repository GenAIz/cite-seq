# Differential Gene Expression Analysis Configuration Parameters

The following parameters are required to run the Differential Gene Expression Analysis Smart Function within the CITE-seq Solution.

## Workspace Configuration Parameters

The following parameters are tied to the Genaiz workspace.

### Config File Name

* Type: String
* Description: Name of the configuration file to load from the `configs/` input port.
* Location: `DGE_CONFIG` environment variable
* Configuration:
  + Workspace entry: Manual entry at creation.
  + Defaults to `config.json` if not set.

### Input Path

* Type: String
* Description: Directory from which input files (input object and configs) are read.
* Location: `SF_INPUT_PATH` environment variable
* Configuration:
  + Set by the Genaiz platform at runtime.
  + Defaults to `run/in` when running locally.
  + Non-mutable at runtime.

### Output Path

* Type: String
* Description: Directory to which the DGE results files are written.
* Location: `output/` subdirectory of the `SF_OUTPUT_PATH` environment variable (the platform collects the output port dataset from there)
* Configuration:
  + Set by the Genaiz platform at runtime.
  + Defaults to `run/out` when running locally.
  + Non-mutable at runtime.

---

## Workflow Configuration Parameters

The following parameters are defined inside `configs/config.json` and control execution of the DGE Smart Function.

### Input RDS

* Required: No
* Type: String
* Description: Path (relative to the input directory) to the annotated `.rds` object.
* Note: when unset or missing, the newest `.rds` found anywhere under the `data` input port is used automatically.
* Location: `input_rds` inside `configs/config.json`
* Example: `"data/annotated_object.rds"`
* Note: the object must carry a log-normalized `data` layer; without one the run fails with `INPUT_MISSING_NORMALIZED_DATA` before compute.

### Output CSV

* Required: No
* Type: String
* Description: Filename for the results CSV; a run timestamp is appended automatically.
* Location: `output_name` inside `configs/config.json`
* Default value: `"dge_results.csv"`
* Note: the module also writes `dge_passthrough_<timestamp>.rds` (a copy of the input object, plus its `_matrices/` directory when present) so downstream nodes can chain from it.
* Note: the former `output_rds` and `output_csv` keys are retired and rejected at validation; use `output_name`.

---

## Testing Strategy Parameters

### Test Method

* Required: No
* Type: String
* Description: Statistical test to use for differential expression (case-insensitive).
* Location: `test_use` inside `configs/config.json`
* Configuration:
  + Single-cell native: `"MAST"`, `"wilcox"`, `"bimod"`, `"roc"`, `"t"`, `"negbinom"`, `"poisson"`, `"LR"`
  + Pseudobulk (bulk RNA-seq models; require `pseudobulking: true`): `"DESeq2"`, `"limma"`, `"edgeR"`
  + Any other value fails validation (`CONFIG_VALIDATION_FAILED`), as does a bulk method without `pseudobulking: true`.
  + Default value: `"MAST"`

### Pseudobulking

* Required: No
* Type: Boolean
* Description: If `true`, aggregates counts per sample per cluster before testing (recommended for multi-sample studies).
* Location: `pseudobulking` inside `configs/config.json`
* Default value: `false`

### Pseudobulk Mode

* Required: No
* Type: String
* Description: How to form pseudobulk units when `pseudobulking` is `true`.
* Location: `pseudobulk_mode` inside `configs/config.json`
* Configuration:
  + Accepted values: `"cluster_sample"` (default), `"sample_only"`, `"celltype_condition"`; any other value fails validation.
  + Note: the current implementation always aggregates one pseudobulk per cluster per sample (per timepoint in the stratified/interaction designs), i.e. the `"cluster_sample"` behaviour. The other two values are accepted by validation but do not change the aggregation yet.

### Design Strategy

* Required: No
* Type: String
* Description: Statistical design used for the pseudobulk model formulas.
* Location: `design_strategy` inside `configs/config.json`
* Configuration:
  + `"stratified"` - Condition effect estimated within clusters (default).
  + `"interaction"` - Includes a cluster-by-condition interaction term.
  + Can be changed between workflow runs.

---

## Column Mapping Parameters

### Cluster Column

* Required: **Yes** - validation fails (`CONFIG_VALIDATION_FAILED`) when it is missing or empty
* Type: String
* Description: Metadata column with cell type or cluster labels; DGE runs separately per unique value.
* Location: `cluster_col` inside `configs/config.json`
* Example: `"SingleR_main"`

### Sample Column

* Required: Only when `pseudobulking` is `true` (validation fails without it); otherwise optional
* Type: String
* Description: Metadata column identifying biological samples or donors.
* Location: `sample_col` inside `configs/config.json`
* Default value: `"Participant ID"`

### Condition Column

* Required: Only when `pseudobulking` is `true` (validation fails without it); otherwise optional
* Type: String
* Description: Metadata column identifying the biological condition to test (e.g., treatment vs. control).
* Note: cell-level runs use it only to resolve `ident_1`/`ident_2` groups.
* Location: `condition_col` inside `configs/config.json`
* Example: `"Terminal Group"`

### Timepoint Column

* Required: No
* Type: String
* Description: Metadata column identifying timepoints, used for longitudinal designs.
* Location: `timepoint_col` inside `configs/config.json`
* Default value: `"Timepoint (months)"`

### Participant Column

* Required: No
* Type: String
* Description: Metadata column identifying participants, used for repeated-measures blocking in pseudobulk designs.
* Location: `participant_col` inside `configs/config.json`
* Default value: same as `sample_col`

### Pseudobulk Assay

* Required: No
* Type: String
* Description: Assay from which counts are aggregated for pseudobulk testing.
* Location: `pseudobulk_assay` inside `configs/config.json`
* Default value: `"RNA"`

---

## Filtering and Statistical Parameters

### Reference Condition

* Required: No
* Type: String or null
* Description: Condition level used as the reference (denominator) in fold changes; `null` uses the first level alphabetically.
* Location: `reference_condition` inside `configs/config.json`
* Default value: `null`

### Minimum Samples Per Group

* Required: No
* Type: Integer
* Description: Minimum samples per condition group for a cluster to be tested; smaller clusters are skipped.
* Location: `min_samples_per_group` inside `configs/config.json`
* Default value: `2`

### Minimum Cells Per Pseudobulk

* Required: No
* Type: Integer
* Description: Minimum cells required to form a pseudobulk unit; smaller combinations are excluded.
* Location: `min_cells_per_pseudobulk` inside `configs/config.json`
* Default value: `10`

### Log Fold Change Threshold

* Required: No
* Type: Float
* Description: Minimum absolute log2 fold change to report; `0` returns all genes, and negative values fail validation.
* Location: `lfc_threshold` inside `configs/config.json`
* Default value: `0`

### P-value Adjustment Method

* Required: No
* Type: String
* Description: Method for multiple testing correction applied to raw p-values.
* Location: `p_adjust_method` inside `configs/config.json`
* Configuration:
  + Accepted values (validated): `"BH"`, `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`, `"BY"`, `"fdr"`, `"none"`. Any other value fails validation.
  + Default value: `"BH"` (Benjamini-Hochberg)

---

## Subset and Comparison Parameters

### Timepoints to Analyze

* Required: No
* Type: Array of strings or numbers
* Description: Restricts pseudobulk runs to these `timepoint_col` values; `null` analyzes all, unknown values fail with `INVALID_TIMEPOINTS`.
* Note: ignored by cell-level methods.
* Location: `timepoints_to_analyze` inside `configs/config.json`
* Default value: `null` (all timepoints)

### Conditions to Analyze

* Required: No
* Type: Array of strings
* Description: Restricts pseudobulk runs to these `condition_col` values; `null` includes all, unknown values fail with `INVALID_CONDITIONS`.
* Note: ignored by cell-level methods.
* Location: `conditions_to_analyze` inside `configs/config.json`
* Default value: `null` (all conditions)

### Ident 1

* Required: No
* Type: String
* Description: First group of the comparison; a value of `condition_col`.
* Note: also accepted spelled `ident.1`.
* Note: set alone, cell-level methods compare this group against all other cells; pseudobulk runs ignore it unless `ident_2` is also set.
* Note: pseudobulk runs validate the values against `condition_col`; cell-level runs also accept `cluster_col` values.
* Note: with both set, cell-level methods run one comparison of the two groups across all cells instead of each cluster against the rest.
* Note: with both set in pseudobulk runs, only the `ident_1` vs `ident_2` contrast is computed (per timepoint in the stratified design) instead of all pairwise contrasts.
* Location: `ident_1` inside `configs/config.json`
* Default value: `null`

### Ident 2

* Required: No
* Type: String
* Description: Second group of the comparison; a value of `condition_col`.
* Note: also accepted spelled `ident.2`; requires `ident_1` (validation fails without it).
* Location: `ident_2` inside `configs/config.json`
* Default value: `null`

### Contrasts

* Required: No
* Type: Object (named list)
* Description: Custom contrast formulas for the pseudobulk models, e.g. `{"GroupA_vs_GroupB": "GroupA - GroupB"}`; `null` auto-generates them.
* Location: `contrasts` inside `configs/config.json`
* Default value: `null` (auto-generated)

---

### Verbose

* Required: No
* Type: Boolean
* Description: If `true`, prints detailed progress messages during analysis.
* Location: `verbose` inside `configs/config.json`
* Default value: `true`

## Common use cases

Each recipe below shows only the keys to set on top of the required
`cluster_col`; recipes that compare condition groups also need `condition_col`.
Everything else keeps its default.

### Marker genes: each cluster vs the rest (default)

Leave `ident_1`/`ident_2` unset and use a cell-level method:

```json
{
  "test_use": "MAST"
}
```

### Two-group comparison

Set `ident_1` and `ident_2` to two values of `condition_col`; one comparison
of the two groups is run across all cells:

```json
{
  "test_use": "wilcox",
  "ident_1": "Treated",
  "ident_2": "Control"
}
```

### One group vs the rest

Set only `ident_1`; cell-level methods compare that group against all other
cells:

```json
{
  "test_use": "wilcox",
  "ident_1": "Treated"
}
```

### Pseudobulk analysis (DESeq2, limma, or edgeR)

Set these four keys together to aggregate counts per cluster per sample and
fit a bulk RNA-seq model (they were previously shown as a
`_pseudobulk_example` block in the example config):

```json
{
  "test_use": "DESeq2",
  "pseudobulking": true,
  "pseudobulk_mode": "cluster_sample",
  "design_strategy": "interaction"
}
```

Pseudobulking aggregates counts per cluster x sample before testing, which is
the statistically appropriate mode when samples, not cells, are the unit of
replication. `"limma"` and `"edgeR"` are drop-in replacements for
`"DESeq2"` here.

### Pseudobulk, one contrast only

Add `ident_1`/`ident_2` to a pseudobulk run to compute just that contrast
instead of all pairwise contrasts:

```json
{
  "test_use": "limma",
  "pseudobulking": true,
  "design_strategy": "stratified",
  "ident_1": "Treated",
  "ident_2": "Control"
}
```

### Restricting pseudobulk runs to timepoints or conditions

Limit a pseudobulk run with `timepoints_to_analyze` and/or
`conditions_to_analyze` (cell-level methods ignore both):

```json
{
  "test_use": "DESeq2",
  "pseudobulking": true,
  "timepoints_to_analyze": [0, 6],
  "conditions_to_analyze": ["Treated", "Control"]
}
```

### Custom pseudobulk contrasts

Replace the auto-generated contrasts with named formulas over
`condition_col` values:

```json
{
  "test_use": "DESeq2",
  "pseudobulking": true,
  "contrasts": { "Treated_vs_Control": "Treated - Control" }
}
```
