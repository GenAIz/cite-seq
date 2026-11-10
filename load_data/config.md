# Load Data Configuration Parameters

The following parameters are required to run the Load Data Smart Function within the CITE-seq Solution.

> **Note:** Metadata integration (adding external sample/patient metadata to the
> object) has been moved to the standalone **add_metadata** module.
> It can be run at any point in the pipeline after data loading.

## Workspace Configuration Parameters

### Config File Name

* Type: String
* Description: Name of the configuration file to load from the `configs/` input port.
* Location: `LOAD_DATA_CONFIG` environment variable
* Configuration: Defaults to `config.json` if not set.

### Input Path

* Type: String
* Description: Directory from which input files and configs are read.
* Location: `SF_INPUT_PATH` environment variable
* Configuration: Set by the Genaiz platform at runtime. Defaults to `run/in` when running locally.

### Output Path

* Type: String
* Description: Directory to which the output files are written.
* Location: `output/` subdirectory of the `SF_OUTPUT_PATH` environment variable (the platform collects the output port dataset from there)
* Configuration: Set by the Genaiz platform at runtime. Defaults to `run/out` when running locally.

---

## Workflow Configuration Parameters

### load_data_path (required)

* Type: String
* Description: Path (relative to the input directory) to ONE pre-processed object (.rds, .h5ad, or .h5mu).
* Note: one file, not a directory; raw formats (.csv, .tar.gz, 10x .h5) and bare matrices are rejected - the input must already be processed.


### project_name

* Type: String
* Default: "CITEseqProject"
* Description: Project name assigned to the output object.

### species_name

* Type: String
* Default: "human"
* Description: Species for Ensembl to gene symbol conversion (`"human"` or `"mouse"`).


### output_name

* Type: String
* Default: "loaded_data.rds"
* Description: Base name for the output `.rds` file; a run timestamp is appended automatically.
* Note: the output dataset holds the `.rds` plus a `<name>_matrices/` directory of its matrix layers; the two must always travel together.

### verbose

* Type: Boolean
* Default: true
* Description: Print progress messages during execution.

---

## h5ad Data Handling

When loading `.h5ad` files, the module automatically:

1. **Maps the file's declared layers onto the output object** without changing
   any values. Raw counts are looked for in `layers/counts`, `layers/raw`, then
   `raw/X`. When one is found, it becomes the object's `counts` layer and the
   file's `X` (when non-integer, i.e. processed expression of the same shape)
   becomes the `data` layer. When the file ships only a non-integer `X`, that
   matrix fills both layers; when `X` is integer-valued raw counts and no
   processed matrix exists, only `counts` is filled and the `data` layer stays
   empty. Every run logs the mapping it chose (`Layer mapping: counts <- ...`).
2. **Extracts cell metadata** (`adata.obs`) and attaches it to `@meta.data`.
3. **Extracts embeddings** (`adata.obsm` - UMAP, PCA) and adds them as reductions on the object.
4. **Detects ADT/protein data** in three locations:
   - Features with `AB_` prefix interleaved in `X` (e.g., COMBAT). The split
     applies to both matrices, so the ADT assay receives its raw counts and,
     when available, its processed values as the ADT `data` slot.
   - Sparse matrix in `uns/antibody_X` or `uns/protein_X` (e.g., Stephenson)
   - Dense matrix in `obsm/protein` or `obsm/ADT`

## h5mu Data Handling

When loading `.h5mu` files, the `RNA` and `PROTEIN`/`ADT` modalities are read
separately and the same layer mapping is applied per modality: each modality's
integer `layers/counts` (or `layers/raw`) becomes the assay's `counts`, and its
`X` becomes the assay's `data` when it holds processed (non-integer) values.
Global `obs` metadata and the RNA modality's `obsm` embeddings are carried onto
the output object. No normalization or other computation is performed - the
module only relabels what the file already contains.
