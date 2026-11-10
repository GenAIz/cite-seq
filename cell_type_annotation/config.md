# Cell Type Annotation Configuration Parameters

The following parameters are required to run the Cell Type Annotation Smart Function within the CITE-seq Solution.

## Workspace Configuration Parameters

The following parameters are tied to the Genaiz workspace.

### Config File Name

* Type: String
* Description: Name of the configuration file to load from the `configs/` input port.
* Location: `ANNOTATION_CONFIG` environment variable
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
* Description: Directory to which the output files are written.
* Location: `output/` subdirectory of the `SF_OUTPUT_PATH` environment variable (the platform collects the output port dataset from there)
* Configuration:
  + Set by the Genaiz platform at runtime.
  + Defaults to `run/out` when running locally.
  + Non-mutable at runtime.

---

## Workflow Configuration Parameters

The following parameters are defined inside `configs/config.json` and control execution of the Cell Type Annotation Smart Function.

### Input RDS

* Required: No
* Type: String
* Description: Path (relative to the input directory) to the processed `.rds` object to annotate.
* Note: when unset or missing, the newest `.rds` found anywhere under the `data` input port is used automatically.
* Location: `input_rds` inside `configs/config.json`
* Example: `"data/processed_object.rds"`

### Output Name

* Required: No
* Type: String
* Description: Filename for the annotated `.rds` output; a run timestamp is appended automatically.
* Note: the output dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers; the two must always travel together.
* Location: `output_name` inside `configs/config.json`
* Default value: `"annotated_object.rds"`

### Method

* Required: No
* Type: String
* Description: Annotation method(s) to run.
* Location: `method` inside `configs/config.json`
* Configuration:
  + `"singler"` - Run only SingleR.
  + `"stcat"` - Run only STCAT.
  + `"both"` - Run both SingleR and STCAT (default).
  + Accepted spellings (validated): `"singler"`/`"SingleR"`, `"stcat"`/`"STCAT"`, `"both"`. Any other value fails validation.

### Reference

* Required: No
* Type: String
* Description: Reference dataset used by SingleR for label transfer.
* Location: `reference` inside `configs/config.json`
* Configuration:
  + `"monaco"` - Monaco Immune Reference (recommended for T/B/NK cells).
  + `"hpca"` - Human Primary Cell Atlas reference.
  + `"blueprint_encode"` - Blueprint/ENCODE reference.
  + `"dice"` - DICE immune reference.
  + `"immgen"` - ImmGen reference (mouse immune cells).
  + Default value: `"monaco"`

### Cluster Column

* Required: No
* Type: String
* Description: Metadata column containing cluster assignments, used to aggregate SingleR scores per cluster.
* Location: `cluster_column` inside `configs/config.json`
* Default value: `"RNA_clusters"` (the column written by the pipeline's clustering step)
* Note: set it to `"seurat_clusters"` if the object was clustered outside this pipeline.

### STCAT Batch Column

* Required: No
* Type: String
* Description: Metadata column identifying batches.
* Note: with more than one batch level, the T-subset embedding sent to STCAT is Harmony-corrected on this column; a single level skips Harmony.
* Location: `stcat_batch_column` inside `configs/config.json`
* Default value: `"orig.ident"`

### Verbose

* Required: No
* Type: Boolean
* Description: If `true`, prints detailed progress messages during annotation.
* Location: `verbose` inside `configs/config.json`
* Default value: `true`

## Advanced parameters

These are not in the example config. They keep their defaults unless you add them:

| Key | Default | Effect |
|---|---|---|
| `assay` | `"RNA"` | Assay to annotate. Change only for sketch-processed objects (`"sketch"`). |
| `singler_label_type` | `["main", "fine"]` | Which SingleR label granularities to compute. SingleR always offers both. |
| `singler_workers` | `4` | Parallel workers for SingleR (BiocParallel). |
| `singler_fine_chunk_size` | `20000` | Cells per chunk in the per-cell fine annotation pass. |
| `stcat_python_script` | `null` | Override path to the STCAT runner. When null the copy shipped in the image is used. |
| `prefix` | `null` | Optional prefix for the added annotation columns. |
| `fallback_column` | `"SingleR_fine"` | Column used for cells STCAT cannot classify. |
| `cleanup_temp` | `true` | Remove the STCAT temp exchange directory afterwards. |

## What the module produces

Annotation runs as a sequence, matching the alternate workflow exactly:

1. **SingleR** assigns main and fine labels per cluster (pseudo-bulk).
2. **Per-cell fine refinement** re-annotates fine labels cell by cell in
   chunks, so a mixed cluster is no longer blanketed with a single label.
3. **Vocabulary rescue** relabels platelets and erythrocytes from canonical
   markers - types the Monaco reference cannot name, which SingleR otherwise
   files under monocytes. Overrides are recorded in `vocab_rescue`.
4. **T cell gate**: only cells whose `SingleR_main` OR `SingleR_fine` label
   references T cells are sent to STCAT; everything else keeps the fallback
   label. STCAT receives a T-subset re-embedding (subset HVGs -> PCA -> Harmony
   when `stcat_batch_column` has more than one level).
5. **Label harmonization** renames Monaco fallback labels that duplicate STCAT
   concepts to the STCAT vocabulary, so one cell type has one name.

Metadata columns on the returned object:

| Column | Meaning |
|---|---|
| `SingleR_main`, `SingleR_main_score` | Cluster-level main label and confidence. |
| `SingleR_fine`, `SingleR_fine_score` | Per-cell fine label and confidence. Reference-faithful Monaco names (e.g. `Vd2 gd T cells`), never harmonized. |
| `vocab_rescue` | `"Platelets"` / `"Erythrocytes"` for rescued cells, `NA` otherwise. Also flags platelet-monocyte aggregates, which some analyses prefer to exclude. |
| `Prediction` | Final cell type: STCAT's call for gated T cells, the fallback label otherwise. |
| `Uncertainty_score` | STCAT uncertainty **and** a provenance flag: non-`NA` means the label came from STCAT, `NA` means it came from SingleR. |

Working columns used during annotation (`SingleR_*_delta`,
`SingleR_fine_cluster*`, `Cluster`, `vocab_rescue_score`, `seurat_clusters`)
are removed before the object is written. The column configured as
`cluster_column` is always kept, even when it is `seurat_clusters`.

