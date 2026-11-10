# T cell Workflow - Comprehensive Guide

**Version:** 1.1 | **Status:** Production Ready | **Setup Time:** 5-10 minutes

---

## Quick Overview

This is a complete CITE-seq analysis pipeline that:
1. Loads RNA/ADT data (H5AD, H5, RDS, CSV, H5MU formats)
2. Adds metadata (optional)
3. Processes & clusters cells with Seurat
4. Annotates cell types (SingleR + STCAT)
5. Aggregates data for machine learning

**Expected runtime:** 30 minutes (10k cells) to 5+ hours (500k cells)

---

## Before You Start

### Requirements
- R 4.0+ installed
- Conda or Docker
- 50+ GB disk space (first time)
- Your CITE-seq data files

### Check Your System
```bash
R --version # R installed?
conda --version # Conda installed?
docker --version # Docker installed? (optional)
df -h # Enough disk space?
```

---

## Part 1: Setup (5 minutes)

### Step 1: Clone the Repository
```bash
git clone <your-gitlab-repo-url>
cd cite-seq/alternate_workflow
```

### Step 2: Create Config File
```bash
# Copy the template config
cp config_standard.json my_config.json
```

Edit `my_config.json` to set your data paths and parameters (see Part 3 below). The directory structure (`script/`, `python/`, `data/`, `results/`, `logs/`) is managed by Docker volume mounts - no manual setup required.

### Step 3: Verify Setup
```bash
# Run verification (10-point check)
bash verify.sh

# Expected output:
# [OK] Passed: 10/10
# Setup is complete!
```

---

## Part 2: Prepare Your Data

### Data Directory Structure
```
data/raw/
+-- rna_counts.h5ad # Required
+-- adt_counts.h5ad # Optional
+-- metadata.xlsx # Optional
```

### Supported Formats
- **RNA/ADT:** H5AD, H5, RDS, CSV, H5MU
- **Metadata:** Excel (.xlsx), CSV, TSV, RDS

### Copy Your Files
```bash
cp /path/to/rna_counts.h5ad data/raw/
cp /path/to/adt_counts.h5ad data/raw/
cp /path/to/metadata.xlsx data/raw/
```

---

## Part 3: Configuration (5 minutes)

### Create Config File
```bash
# Copy template
cp config_standard.json my_config.json

# Edit with your settings
nano my_config.json
```

Note: `my_config.json` is excluded from git (it contains machine-specific paths).

### Essential Settings to Change

```json
{
  "metadata": {
    "project_name": "MyProject"
  },
  "paths": {
    "input_data": {
      "rna_path": "/absolute/path/to/rna_counts.h5ad",
      "adt_path": "/absolute/path/to/adt_counts.h5ad",
      "metadata_file": "/absolute/path/to/metadata.xlsx"
    },
    "output_dir": "/results",
    "logs_dir": "/logs"
  },
  "output_files": {
    "load_data": "MyProject_load_data.rds",
    "add_metadata": "MyProject_add_metadata.rds",
    "process_seurat": "MyProject_process_seurat.rds",
    "cell_type_annotation": "MyProject_cell_type_annotation.rds",
    "aggregate_ml": "MyProject_aggregate_ml.xlsx"
  },
  "load_data": {
    "species_name": "human"
  },
  "execution": {
    "mode": "step-by-step",
    "n_cores": 8
  }
}
```

The `output_files` section is optional - the example above shows renames. If omitted, the defaults are `load_data.rds`, `add_metadata.rds`, `process_seurat.rds`, `cell_type_annotation.rds` and `aggregate_ml.xlsx` (the run timestamp is always appended, and resume mode looks for these same names). Only include the keys you want to rename.

### Get Absolute Paths
```bash
# Use pwd to get full path
pwd # Shows: /home/username/cite-seq/alternate_workflow

# Full path example:
# /home/username/cite-seq/alternate_workflow/raw/rna_counts.h5ad
# NOT: ~/raw/rna_counts.h5ad (tilde doesn't work)
# NOT: data/raw/rna_counts.h5ad (relative paths don't work)
```

### Validate Config
```bash
# Test without running pipeline
Rscript orchestrator.R --config my_config.json --dry-run

# Expected output:
# [OK] All validation checks passed!
```

---

## Part 4: Run the Pipeline

### Option A: Local (Command Line)

#### Create Python Environment
```bash
# Install dependencies
conda env create -f environment.yml

# Activate
conda activate tcell-env
```

#### Run Pipeline
```bash
# Test first (dry-run)
Rscript orchestrator.R --config my_config.json --dry-run

# Run full pipeline
Rscript orchestrator.R --config my_config.json

# Monitor progress
tail -f logs/tcell_workflow_*.log
```

#### Output Files

A **run timestamp** (`YYYYMMDD_HHMMSS`) is **always** appended to every output filename for compliance - regardless of mode. The base name comes from `output_files` in your config, or uses the defaults shown below.

**step-by-step mode** - saves a timestamped RDS after every step that was run:
```
results/
+-- load_data_20260423_180001.rds
+-- add_metadata_20260423_180001.rds
+-- process_seurat_20260423_180001.rds
+-- cell_type_annotation_20260423_180001.rds
+-- aggregate_ml_20260423_180001.xlsx (aggregate_ml: XLSX only, no RDS)
```

**chained mode** - saves only after cell type annotation and ML aggregation:
```
results/
+-- annotated_seurat_20260423_180001.rds
+-- aggregate_ml_20260423_180001.xlsx
```

**logs** (both modes):
```
logs/
+-- tcell_workflow_TIMESTAMP.log <- includes a full snapshot of the config used
```

When resuming a pipeline (`--resume <step>`), the orchestrator automatically finds the most recent timestamped intermediate file for each prior step.

---

### Option B: Docker (Recommended for reproducibility)

#### Build Docker Image
```bash
# Method 1: Using Makefile (easiest)
make build

# Method 2: Manual Docker command
docker build -t t_cell_workflow:latest .

# Verify build
docker images | grep t_cell_workflow
# Should show: t_cell_workflow ~2.5GB
```

#### Run Pipeline in Docker
```bash
# Dry-run test
make test CONFIG=my_config.json

# Full pipeline (default: 8 CPUs)
make run CONFIG=my_config.json

# Full pipeline with custom CPU count
make run CONFIG=my_config.json CPUS=4

# Interactive shell
make shell

# Cleanup
make clean
```

**What Docker does:**
- Installs R 4.5.1
- Installs all system dependencies
- Creates the `tcell-env` Python environment from `environment.yml`
- Installs all R and Python packages
- Copies `script/`, `python/`, and `orchestrator.R` into the image
- Runs completely isolated from your system

#### CPU Configuration

The Makefile allocates **8 CPUs by default**. If your machine has fewer CPUs (like 4), you must override this:

```bash
# Run with 4 CPUs (most common for modern machines)
make run CPUS=4 CONFIG=my_config.json

# Run with 2 CPUs (for limited hardware)
make run CPUS=2
```

**Check your available CPUs:**

macOS/Linux:
```bash
nproc
```

Windows/WSL:
```bash
wsl -e nproc
```

Docker Desktop Settings:
- macOS: Docker menu -> Preferences -> Resources -> CPUs
- Windows: Docker Desktop Settings -> Resources -> CPUs
- Look for the "CPUs" slider

#### Other Makefile Variables

**CONFIG** - Config file to use. Accepts a filename (with optional `CONFIG_DIR`) or a full absolute path:
```bash
# Filename only - looks in the current directory
make run CONFIG=my_config.json

# Filename + directory - looks in CONFIG_DIR
make run CONFIG=my_config.json CONFIG_DIR=/home/user/dev/citeseq_data_processing

# Full absolute path - CONFIG_DIR is ignored
make run CONFIG=/home/user/dev/citeseq_data_processing/my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/home/user/dev/citeseq_data_processing/results
```

**CONFIG_DIR** - Directory containing the config file (default: current directory). Use when CONFIG is just a filename but lives elsewhere.

**DATA_DIR** - The host directory mounted as `/raw` inside the container (default: current directory). The whole directory is mounted, not a `data/` subfolder. See [Docker Path Resolution](#docker-path-resolution-and-data_dir) below for a full explanation.
```bash
# On your local machine (data lives next to the Makefile)
make run CONFIG=my_config.json

# On a research/HPC machine (data lives somewhere else entirely)
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq
```

**OUTPUT_DIR** - The host directory mounted as `/results` inside the container. Results are written here (default: current directory).
```bash
# Write results to a dedicated location on a research machine
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq/results
```

**LOGS_DIR** - The host directory for log files, mounted as `/logs` inside the container (default: `OUTPUT_DIR/logs`). Set this explicitly if you want logs somewhere other than a `logs/` subfolder inside your output directory.
```bash
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq/results LOGS_DIR=/mnt/citeseq/logs
```

**Combine all variables:**
```bash
make run CPUS=4 CONFIG=/mnt/citeseq/data/my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq/results LOGS_DIR=/mnt/citeseq/logs
```

#### Memory Configuration

Default memory allocation: **90 GB** (set by `MEM` in the `Makefile`). Override at runtime with `MEM=`:

```bash
make run CONFIG=my_config.json MEM=48g # decrease on smaller machines
make run CONFIG=my_config.json MEM=8g # reduce to force OOM faster (debugging)
```

#### File Ownership and User Permissions

By default Docker runs as root, meaning files written to mounted volumes (`results/`, `logs/`) end up owned by root on your host machine - making them hard to delete or edit without `sudo`.

The Makefile automatically passes your host UID and GID into the container so output files are always owned by you:

```makefile
-e HOST_UID=$(shell id -u)
-e HOST_GID=$(shell id -g)
```

The `entrypoint.sh` uses `gosu` to create a matching user inside the container and run the pipeline as that user. No action required - this happens automatically on every `make run`.

#### Docker Compose Alternative

If you prefer using `docker-compose`, create a `docker-compose.yml`:

```yaml
version: '3.8'

services:
  workflow:
    image: t_cell_workflow:latest
    volumes:
      - ./data:/data
      - .:/output
      - ./my_config.json:/app/my_config.json:ro
    environment:
      - CONFIG=my_config.json
    cpus: '4'
    mem_limit: 32g
    shm_size: 8gb
```

Then run:
```bash
docker-compose run --rm workflow run my_config.json
```

---

#### Docker Path Resolution and DATA_DIR

This is the most common source of confusion when running the pipeline on different machines. Understanding it will save you significant debugging time.

**The core rule: config files always use container-internal paths.**

When Docker runs the pipeline, your files are not accessed directly at their host locations. Instead, Docker *mounts* specific paths from your host into the container. The Makefile defines three mounts:

```makefile
-v $(DATA_DIR)/data:/data # your data directory -> /data inside the container
-v $(OUTPUT_DIR):/output # your output directory -> /output inside the container
-v $(PWD)/$(CONFIG):/app/$(CONFIG):ro # your config file only -> /app/<config> (read-only)
```

This means:
- `$(DATA_DIR)` on your host machine becomes `/raw` inside the container
- `$(OUTPUT_DIR)` becomes `/results`, and `$(LOGS_DIR)` becomes `/logs` - this is where results and logs are written
- **Only your config file** is mounted into `/app` - all pipeline code (`orchestrator.R`, `script/`, `python/`, `entrypoint.sh`) is baked into the image at build time and is never exposed as a host-side mount
- Your config file must always reference `/raw/...`, `/results/...` and `/logs/...` paths - never your host machine's paths

**Example: local machine**

Your data and project live at `~/cite-seq/alternate_workflow/`. You run:
```bash
make run CONFIG=my_config.json
# DATA_DIR and OUTPUT_DIR both default to $(PWD) = ~/cite-seq/alternate_workflow
# Mounts:
# ~/cite-seq/alternate_workflow/data -> /data (input)
# ~/cite-seq/alternate_workflow -> /output (results go to ./results/ and ./logs/)
# ~/cite-seq/alternate_workflow/my_config.json -> /app/my_config.json (config only, read-only)
```

Your `my_config.json` must say:
```json
"rna_path": "/raw/combat.h5ad",
"output_dir": "/results",
"logs_dir": "/logs"
```
Not host paths like `/home/username/cite-seq/alternate_workflow/data/raw/...`.

**Example: research machine with separate output location**

Your data lives at `/mnt/citeseq/data/raw/` and you want results written to `/mnt/citeseq/results/`. You run:
```bash
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq
# Mounts:
# /mnt/citeseq/data -> /data (input)
# /mnt/citeseq -> /output (results + logs)
# ~/cite-seq/alternate_workflow/my_config.json -> /app/my_config.json (config only, read-only)
```

Your `my_config.json` is identical - no changes needed:
```json
"rna_path": "/raw/combat.h5ad",
"output_dir": "/results",
"logs_dir": "/logs"
```

**The config file never changes between machines. Only `DATA_DIR` and `OUTPUT_DIR` change. The config is portable; the mount points adapt to the machine.**

**What happens if you put the host path in the config?**

You'll see an error like:
```
Error: RNA file not found: /mnt/citeseq/raw/combat.h5ad
```

That path doesn't exist *inside the container* - only `/raw/combat.h5ad` does. The validation module (`validation_module.R`) calls `file.exists()` from inside the container, so it can only see container-internal paths.

**Summary table:**

| Scenario | Command | `rna_path` | `output_dir` | Code visible to user? |
|---|---|---|---|---|
| Local | `make run CONFIG=my_config.json` | `/raw/combat.h5ad` | `/results` | No - config only |
| Research machine | `make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq` | `/raw/combat.h5ad` | `/results` | No - config only |
| Custom output | `make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/scratch/myproject` | `/raw/combat.h5ad` | `/results` | No - config only |
| Debug shell | `make shell` | - | - | Yes - full `/app` access |

Note: `make shell` intentionally mounts the full `/app` directory for your own debugging. Do not share `make shell` access with external users.

---

#### Docker Troubleshooting

**"range of CPUs is from 0.01 to 4.00"**
Your machine has only 4 CPUs available. Use:
```bash
make run CPUS=4
```

**"Cannot allocate memory"**
Reduce CPU allocation or memory:
```bash
make run CPUS=2
```
Or edit Makefile to use `--memory=16g` instead of `32g`.

**"Cannot connect to Docker daemon"**
Make sure Docker is running:
```bash
docker ps
```

---

## Advanced Options

### Resume from Checkpoint
If the pipeline stops, resume from a specific step:
```bash
Rscript orchestrator.R --config my_config.json --resume cell_type_annotation
```

Available resume points:
- `load_data`
- `add_metadata`
- `process_seurat`
- `cell_type_annotation`
- `aggregate_ml`

### Run Individual Steps
```bash
# Modify config to run only specific steps
nano my_config.json
# Change: "steps": ["load_data", "process_seurat"]
```

### Adjust Performance
```json
{
  "execution": {
    "n_cores": 2,
    "mode": "chained"
  },
  "processing": {
    "features": {"n_features": 2000}
  }
}
```

### Rename Output Files
```json
{
  "output_files": {
    "cell_type_annotation": "annotated_seurat.rds",
    "aggregate_ml": "ml_ready_aggregate.xlsx"
  }
}
```
A timestamp is automatically appended in every mode: `annotated_seurat_20260423_180001.rds`. Note that `aggregate_ml` always produces an XLSX workbook - any other extension in its configured name is normalised to `.xlsx`.

---

## Troubleshooting

### Issue: "Cannot find file"
```bash
# Check file exists
ls /absolute/path/to/file.h5ad

# Get correct absolute path
realpath data/raw/rna.h5ad

# Update config
nano my_config.json
```

### Issue: "Out of memory"
```json
{
  "processing": {
    "features": {"n_features": 2000}
  }
}
```
Or use a machine with more RAM, or process a smaller dataset.

### Issue: "Cannot find script files"
```bash
# Verify structure
ls script/
# Should show 8 R files

ls python/
# Should show run_stcat.py

# If missing, re-run setup
bash setup.sh
```

### Issue: "Docker daemon not running"
```bash
# Start Docker
sudo systemctl start docker # Linux
open /Applications/Docker.app # Mac

# Or run without Docker
conda activate tcell-env
Rscript orchestrator.R --config my_config.json
```

### Issue: Python/STCAT packages not found
```bash
# Ensure environment is created and activated
conda env create -f environment.yml
conda activate tcell-env

# Verify
python -c "import scanpy; import STCAT; print('OK')"
```

### Issue: "libglpk.so.40: cannot open shared object"
This error occurs when running on a different machine where the Docker image needs to be rebuilt:

```
ERROR: unable to load shared object '/usr/local/lib/R/site-library/igraph/libs/igraph.so':
  libglpk.so.40: cannot open shared object file: No such file or directory
```

**Solution: Rebuild the Docker image on your machine**
```bash
# This compiles the R packages with your system's libraries
make build
```

### Check Logs for Errors

Log files are named per project (`tcell_workflow_<project>_<timestamp>.log`) and now capture **all** console output - QC counts, ADT cleaning, batch detection, STCAT output, etc. - not just milestone lines, so they're a complete trace of a run.

```bash
# View log file (project-named)
cat logs/tcell_workflow_*.log

# Find errors
grep -i "error\|fail" logs/tcell_workflow_*.log

# View last 50 lines
tail -50 logs/tcell_workflow_*.log
```

---

## File Structure

```
alternate_workflow/
+-- orchestrator.R * main entry point
+-- comprehensive_guide.md this file
+-- Dockerfile container definition
+-- Makefile Docker automation
+-- entrypoint.sh Docker entrypoint
+-- environment.yml Python environment spec
+-- setup.sh run once to set up
+-- verify.sh 10-point verification
+-- config_standard.json default config template
+-- config_minimal.json minimal config template
+-- config_advanced.json advanced config template
+-- config_schema.json config reference (documentation only - NOT enforced)
+-- my_config.json your config (not in git - edit this)
|
+-- python/
| +-- run_stcat.py STCAT cell type annotation script
|
+-- script/
| +-- logging_module.R
| +-- validation_module.R
| +-- load_data_updated.R
| +-- add_metadata.R
| +-- process_seurat_pipeline.R
| +-- cell_type_annotation.R
| +-- aggregate_obj_ML_processing.R
| +-- export_loupe.R Loupe (.cloupe) export step
|
+-- data/ not in git
| +-- raw/
| +-- rna_counts.h5ad your RNA data
| +-- adt_counts.h5ad your ADT data (optional)
| +-- metadata.xlsx sample metadata (optional)
|
+-- results/ not in git
| | step-by-step mode:
| +-- load_data_TIMESTAMP.rds
| +-- add_metadata_TIMESTAMP.rds
| +-- process_seurat_TIMESTAMP.rds
| +-- cell_type_annotation_TIMESTAMP.rds
| +-- <aggregate_ml>_TIMESTAMP.xlsx <- aggregate_ml writes ONLY this xlsx (no rds); name from output_files.aggregate_ml
| chained mode:
| +-- annotated_seurat_YYYYMMDD_HHMMSS.rds
| +-- <aggregate_ml>_TIMESTAMP.xlsx
|
+-- logs/ not in git
    +-- tcell_workflow_PROJECT_TIMESTAMP.log <- project-named; full console capture + config snapshot
                                            (+ <project>_TIMESTAMP_RNA.cloupe / _ADT.cloupe in results/ if export_loupe is enabled)
```

---

## Expected Runtimes

| Data Size | Time | RAM | CPU |
|-----------|------|-----|-----|
| 10k cells | 30 min | 8 GB | 4 cores |
| 50k cells | 1-2 hours | 16 GB | 8 cores |
| 100k cells | 2-4 hours | 32 GB | 16 cores |
| 500k cells | 5+ hours | 64 GB | 32 cores |

---

## Quick Commands Reference

```bash
# Setup
bash setup.sh && bash verify.sh

# Test configuration
Rscript orchestrator.R --config my_config.json --dry-run

# Run pipeline (local, no Docker)
conda activate tcell-env
Rscript orchestrator.R --config my_config.json

# Resume from checkpoint
Rscript orchestrator.R --config my_config.json --resume cell_type_annotation

# Docker: Build
make build

# Docker: Validate config only
make validate CONFIG=my_config.json

# Docker: Dry-run (simulate without executing)
make dry-run CONFIG=my_config.json

# Docker: Run on local machine (data in ./data/, results in ./results/)
make run CONFIG=my_config.json

# Docker: Run on research machine
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq

# Docker: Run with separate output location
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/scratch/myproject

# Docker: Run with custom CPU count
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq OUTPUT_DIR=/mnt/citeseq CPUS=4

# Docker: Interactive shell (for debugging)
make shell

# Docker: Cleanup stopped containers
make clean

# View logs
tail -f logs/tcell_workflow_*.log
```

---

## Next Steps After Pipeline Completes

### Load Results in R
```r
library(Seurat)
library(readxl)

# step-by-step mode - outputs are always timestamped; glob for the latest
f <- Sys.glob("results/cell_type_annotation_*.rds")
obj <- readRDS(f[which.max(file.mtime(f))])

# chained mode - same approach with the chained name
f <- Sys.glob("results/annotated_seurat_*.rds")
obj <- readRDS(f[which.max(file.mtime(f))])

# Check dimensions
dim(obj)

# View cell types
table(obj@meta.data$SingleR_main)
table(obj@meta.data$Prediction) # STCAT annotations

# Plot
DimPlot(obj, group.by = "SingleR_main", label = TRUE)
DimPlot(obj, group.by = "Prediction", label = TRUE)

# Load ML-ready data
ml_data <- read_excel(Sys.glob("results/*_aggregate_ml_*.xlsx")[1])
```

### Save Publication Figures
```r
ggsave("umap.pdf", p, width = 10, height = 8)
ggsave("celltype_composition.png", p, width = 8, height = 6, dpi = 300)
```

---

## Configuration Reference

This section provides a complete guide to all configuration parameters in `my_config.json`. Each parameter controls a specific aspect of the pipeline.

### Output Files Section

Controls the filename of each step's saved output. All keys are optional - omit any key to use the default name: `load_data.rds`, `add_metadata.rds`, `process_seurat.rds`, `cell_type_annotation.rds`, `aggregate_ml.xlsx`. Resume mode searches for these same base names, so a resume finds prior outputs whether or not you configure renames.

```json
"output_files": {
  "load_data": "MyProject_load_data.rds",
  "add_metadata": "MyProject_add_metadata.rds",
  "process_seurat": "MyProject_process_seurat.rds",
  "cell_type_annotation": "MyProject_cell_type_annotation.rds",
  "aggregate_ml": "MyProject_aggregate_ml.xlsx"
}
```

The run timestamp is appended before the extension in **every mode**. In **chained mode**, only the keys `cell_type_annotation` (default `annotated_seurat.rds`) and `aggregate_ml` are used:

```json
"output_files": {
  "cell_type_annotation": "annotated_seurat.rds",
  "aggregate_ml": "ml_ready_aggregate.xlsx"
}
```
Produces: `annotated_seurat_20260423_180001.rds` and `ml_ready_aggregate_20260423_180001.xlsx`. The `aggregate_ml` step always writes an XLSX workbook - any other extension in its configured name is normalised to `.xlsx`, and no RDS is saved for that step.

The timestamp prevents previous runs from being overwritten when you re-run with the same config.

---

### Metadata Section

Controls project-level information:

```json
"metadata": {
  "project_name": "MyProject",
  "description": "Optional description of your analysis"
}
```

- **project_name** - Name of your project; used in the run log filename (`tcell_workflow_MyProject_<timestamp>.log`) and as the default base name for Loupe exports
- **description** - Optional description for documentation purposes

### Paths Section

Specifies where data is located and where results should be saved:

```json
"paths": {
  "input_data": {
    "rna_path": "/absolute/path/to/rna_counts.h5ad",
    "adt_path": "/absolute/path/to/adt_counts.h5ad",
    "metadata_file": "/absolute/path/to/metadata.xlsx"
  },
  "output_dir": "/absolute/path/to/results",
  "logs_dir": "/absolute/path/to/logs"
}
```

- **rna_path** - Path to RNA count matrix (required). Supports: H5AD, H5, RDS, CSV, H5MU
- **adt_path** - Path to ADT/protein count matrix (optional). Supports same formats as RNA
- **metadata_file** - Path to sample metadata (optional). Excel (.xlsx), CSV, TSV, or RDS format
- **output_dir** - Container-internal path for results (use `/results`). Maps to `OUTPUT_DIR` on your host via Docker volume mount
- **logs_dir** - Container-internal path for logs (use `/logs`). Maps to the same `OUTPUT_DIR` host directory

**Important - Docker vs. local paths:**

When running via Docker (`make run`), paths in the config must be *container-internal* paths, not your host machine's paths. Docker mounts your data directory as `/raw` inside the container, so all input paths must start with `/raw/`:

```json
"rna_path": "/raw/combat.h5ad" [OK] container-internal - always correct
"rna_path": "/mnt/citeseq/data/raw/..." [X] host path - invisible inside container
"rna_path": "data/raw/rna.h5ad" [X] relative path - never works
```

Use the `DATA_DIR` Makefile variable to tell Docker *where* your data directory lives on the host. See the [Docker Path Resolution](#docker-path-resolution-and-data_dir) section for full details and examples.

### Execution Section

Controls how the pipeline runs:

```json
"execution": {
  "mode": "step-by-step",
  "steps": ["load_data", "process_seurat", "cell_type_annotation"],
  "verbose": true,
  "n_cores": 8
}
```

- **mode** - Controls what gets saved to disk:
  - `"step-by-step"` - saves a timestamped RDS after every Seurat-object step that runs (enables resuming from any checkpoint). The aggregate_ml and export_loupe steps write their own deliverables (XLSX / .cloupe), never an RDS.
  - `"chained"` - saves an RDS only after cell type annotation (`annotated_seurat_<timestamp>.rds`); the ML aggregation XLSX is written in every mode. Intermediates are held in memory only.
- **steps** - Array of steps to execute. Options:
  - `"load_data"` - Load RNA/ADT counts and optionally merge metadata
  - `"add_metadata"` - Add additional metadata columns to cells
  - `"process_seurat"` - Normalize, filter, find variable genes, reduce dimensions, cluster
  - `"cell_type_annotation"` - Annotate cell types using SingleR and/or STCAT
  - `"aggregate_ml"` - Aggregate data for machine learning (cell counts by type per sample)
  - `"export_loupe"` - *(optional)* Export `.cloupe` files (one per modality) from the annotated single-cell object for 10x Loupe Browser. Only runs when listed here, and reads the `cell_type_annotation` output (not the pseudobulk aggregate). See the **Export Loupe Section** below.

  Run all steps: `["load_data", "add_metadata", "process_seurat", "cell_type_annotation", "aggregate_ml"]`

  Add Loupe export at the end: `[..., "aggregate_ml", "export_loupe"]`

  Run specific steps only: `["load_data", "process_seurat"]`

- **verbose** - `true` for detailed progress messages, `false` for minimal output
- **n_cores** - Number of CPU cores used for parallel processing in SingleR annotation. Set to your machine's CPU count for best performance. Default: 8

### Load Data Section

Controls initial data loading:

```json
"load_data": {
  "project_name": "TestProject",
  "species_name": "human"
}
```

- **project_name** - Name for this dataset (shown in plots and outputs)
- **species_name** - `"human"` or `"mouse"` (used for mitochondrial and ribosomal gene detection)

### Processing Section

Controls all data processing parameters. This is the largest section:

#### Normalization

```json
"processing": {
  "normalization": {
    "rna_method": "LogNormalize",
    "adt_method": "CLR"
  },
  ...
}
```

- **rna_method** - Normalization for RNA counts. Options:
  - `"LogNormalize"` - Log normalization (default, most common)
  - `"SCT"` - Seurat's SCTransform (alternative, more sophisticated)

- **adt_method** - Normalization for ADT counts. Options:
  - `"CLR"` - Centered log ratio (default, recommended for ADT)
  - `"LogNormalize"` - Log normalization (alternative)

  *Note:* if the loaded ADT counts are already normalized (non-integer values, e.g. from a CELLxGENE-schema h5ad), CLR is **skipped automatically** and the provided values are used as-is, to avoid double-normalization.

#### Quality Control Filters

```json
"qc_filters": {
  "min_features": 200,
  "max_features": 6000,
  "max_percent_mt": 20,
  "max_percent_rp": 40
}
```

- **min_features** - Minimum genes detected per cell. Cells with fewer genes are removed (default: 200)
- **max_features** - Maximum genes detected per cell. Cells with more genes are potential doublets and removed (default: 6000)
- **max_percent_mt** - Maximum % mitochondrial reads. Cells exceeding this are removed (default: 20%)
- **max_percent_rp** - Maximum % ribosomal reads. Cells exceeding this are removed (default: 40%)
- **min_cells** - Minimum number of cells a gene must be detected in to be kept. Genes below this are removed *before* normalization, shrinking the matrix and lowering peak RAM (default: 3)
- **min_cells_per_sample** - Minimum cells a batch (`batch.column`, e.g. `orig.ident`) must have *after QC* to be kept. Batches below this are dropped entirely, since tiny samples give degenerate per-sample ADT densities that break ADTnorm and yield unreliable batch metrics (default: 500)

##### ADT panel cleaning (`processing.adt_filters`)

Non-marker features are removed from the ADT assay **before** normalization, so ADT CLR, PCA, clustering, WNN, and ADTnorm all operate on genuine surface markers.

**Sample-multiplexing tags are detected and removed automatically - no configuration required.** Hashtag/CellPlex/MULTI-seq oligos share the antibody chemistry, so they arrive inside the antibody matrix and are distinguishable only by name; left there they are scored as if they were surface markers. The built-in pattern covers `HTO*`, `Hashtag*`, `HT<n>`, `BC<n>`, `TotalSeq*Hashtag*`, `Set<n>-Sample<n>`, `Sample<n>`, `CMO*`, `LMO*` and `MULTI-*`. It is deliberately narrow: validated against three real datasets totalling ~90,000 features - including a 161-marker panel that carried 24 `Set1-Sample*` tags and a 192-marker `AB_`-prefixed panel - it caught every tag present and produced **zero** false positives on genuine markers such as `Ab.CD3-UCHT1`, `Ab.CD14-M5E2`, `Ab.CD141`, `Ab.CD303`, `Ab.Isotype-RTK2071` and `AB_CD56_NCAM`. Removal is skipped and logged if it would leave fewer than two features. This is panel hygiene only - the pipeline does not demultiplex, and these features are simply discarded.

`keep_pattern` and `drop_pattern` below are **additive** to that built-in pass. Matching is case-insensitive regex.

- **keep_pattern** - If set, keep only ADT features matching this regex (e.g. `"^Ab\\."`); `null` keeps all (default: `null`)
- **drop_pattern** - Additionally remove ADT features matching this regex. Default `"isotype"` drops isotype controls. Multiplexing tags are already handled automatically, so this is only needed for panel-specific oddities (default: `"isotype"`)

#### Feature Selection

```json
"features": {
  "n_features": 5000
}
```

- **n_features** - Number of variable genes to use for downstream analysis (default: 5000). Higher values include more genes but increase computation time. Common range: 2000-5000

#### Dimensionality Reduction

```json
"reduction": {
  "pca_dims_rna": 30,
  "pca_dims_adt": 10
}
```

- **pca_dims_rna** - Number of PCA dimensions to compute for RNA data (default: 30). Used for UMAP and clustering
- **pca_dims_adt** - Number of PCA dimensions to compute for ADT data (default: 10). **Important:** Must be less than the number of antibodies in your dataset (e.g., if you have 13 ADTs, max is 13)

#### Clustering

```json
"clustering": {
  "rna_resolution": 1.5,
  "wnn_resolution": 1.5,
  "k_param": 30,
  "algorithm": 1
}
```

- **rna_resolution** - Clustering resolution for RNA modality (default: 1.5). Higher values = more clusters (range: 0.1-3.0)
- **wnn_resolution** - Clustering resolution for weighted nearest neighbors combining RNA+ADT (default: 1.5). Higher values = more clusters
- **k_param** - Number of neighbors for KNN graph (default: 30). Higher values create more connected neighborhoods
- **algorithm** - Clustering algorithm to use (default: 1):
  - `1` - Louvain algorithm (faster, default)
  - `4` - Leiden algorithm (more refined). Runs via the memory-safe `igraph` method (avoids densifying the neighbor graph, which would OOM on large datasets) and requires the Python `leidenalg` package in the conda env; if it is unavailable the pipeline automatically falls back to Louvain.

### Cell Type Annotation Section

Controls automated cell type assignment:

```json
"cell_type_annotation": {
  "method": "both",
  "reference": "monaco",
  "cluster_column": "wsnn_clusters",
  "assay": "RNA",
  "singler_label_type": ["main", "fine"],
  "singler_workers": 8,
  "singler_fine_chunk_size": 20000,
  "stcat_python_script": null,
  "stcat_batch_column": "orig.ident",
  "prefix": null,
  "fallback_column": "SingleR_fine",
  "cleanup_temp": true
}
```

- **method** - Which annotation method(s) to use:
  - `"singler"` - Use only SingleR
  - `"stcat"` - Use only STCAT (requires Python script)
  - `"both"` - Use both SingleR and STCAT

- **reference** - Reference dataset for SingleR (case-insensitive). Options:
  - `"monaco"` - Monaco immune reference (recommended for immune cells)
  - `"hpca"` - Human Primary Cell Atlas reference
  - `"blueprint_encode"` - Blueprint/ENCODE reference
  - `"dice"` - DICE immune reference
  - `"immgen"` - Immunological Genome Project reference (mouse)

- **cluster_column** - Metadata column name containing cluster assignments. Default: `"RNA_clusters"`. Must exist in the object's metadata after the processing step (e.g., `"RNA_clusters"`, `"wsnn_clusters"`)

- **assay** - Which modality to use for annotation:
  - `"RNA"` - Use RNA counts (default)
  - `"ADT"` - Use protein counts (if available)

- **singler_label_type** - Which label types to return from SingleR:
  - `["main"]` - Broad cell types only
  - `["main", "fine"]` - Both broad and specific cell types
  - `["fine"]` - Specific types only

- **singler_workers** - Number of parallel workers for SingleR (default: 8). Should match `execution.n_cores`

- **singler_fine_chunk_size** - Cells per chunk in the per-cell fine annotation pass (default: 20000). Memory scales with the chunk, not the dataset; lower it on a small machine

- **stcat_python_script** - Path to custom STCAT script (usually null, uses default)

- **stcat_batch_column** - Metadata column identifying batches for the T-subset re-embedding sent to STCAT (defaults to `batch.column`, then `"orig.ident"`). When it has more than one level, that embedding is Harmony-corrected so STCAT's neighbour graph resolves T substructure instead of separating by batch

- **prefix** - Prefix for column names (e.g., `"stcat"` creates `stcat_Prediction`, `stcat_SingleR_main`, ...). Usually null

- **fallback_column** - Column that fills `Prediction` for cells NOT sent to STCAT (SingleR did not call them T cells) or that STCAT could not classify (default: `"SingleR_fine"`). Those cells get `Uncertainty_score = NA`

- **cleanup_temp** - `true` to remove temporary files after annotation (default: true)

**Notes:**
- With `method: "both"`, annotation runs as a sequence: SingleR (cluster-level main + fine) -> per-cell fine refinement (chunked) -> vocabulary rescue (Platelets/Erythrocytes from canonical markers, recorded in `vocab_rescue`) -> T cell gate (only cells whose `SingleR_main` OR `SingleR_fine` references T cells are sent to STCAT, on a T-subset re-embedding) -> SingleR fallback for everything else -> label harmonization to the STCAT vocabulary.
- Final metadata columns: `SingleR_main`, `SingleR_main_score`, `SingleR_fine`, `SingleR_fine_score` (per-cell, reference-faithful Monaco names, never harmonized), `vocab_rescue`, `Prediction`, and `Uncertainty_score` - which doubles as a provenance flag: non-`NA` means the label came from STCAT, `NA` means it came from SingleR. Working columns (`SingleR_*_delta`, `SingleR_fine_cluster*`, `Cluster`, `vocab_rescue_score`, `seurat_clusters`) are dropped before the object is written; the configured `cluster_column` is always kept.

### Add Metadata Section

Optional: Merge additional metadata with cells:

```json
"add_metadata": {
  "enabled": false,
  "merge_by": "MULTI_ID",
  "merge_by_metadata": "Specimen Barcode"
}
```

- **enabled** - `true` to merge external metadata, `false` to skip
- **merge_by** - Column name in your Seurat object to match on
- **merge_by_metadata** - Column name in your metadata file to match on

### Aggregation Section

Controls aggregation for downstream analysis:

```json
"aggregation": {
  "participant_col": "Participant ID",
  "timepoint_col": "Timepoint",
  "celltype_cols": ["Prediction"],
  "preserve_columns": ["Sample", "Treatment"]
}
```

- **participant_col** - Metadata column identifying unique participants/samples
- **timepoint_col** - Metadata column for timepoint/condition (optional)
- **celltype_cols** - Which cell type column(s) to aggregate by. Default: `["Prediction"]` (from STCAT). Could also be `["SingleR_main"]` for SingleR results
- **preserve_columns** - Additional metadata columns to preserve in aggregated output
**Output:** this step writes **only** an Excel file - no RDS and no CSV. Its name is taken from `output_files.aggregate_ml` (any extension is normalised to `.xlsx`), with the run timestamp inserted before the extension, e.g. `citeseq_rerun_aggregate_ml_YYYYMMDD_HHMMSS.xlsx`. Set the name you want in the `output_files` block.

---

## Export Loupe Section *(optional)*

Enables the `export_loupe` step (which must also be added to `execution.steps`). It loads the annotated single-cell object and writes one `.cloupe` per modality for 10x Loupe Browser.

```json
"export_loupe": {
  "source_step": "cell_type_annotation",
  "assays": ["RNA", "ADT"],
  "output_name": null,
  "executable_path": null
}
```

- **source_step** - Which step's output to load. Default `"cell_type_annotation"` (the single-cell annotated object). Do **not** point at `aggregate_ml` - that output is pseudobulk, not single cells.
- **assays** - Modalities to export, one `.cloupe` each. Default `["RNA", "ADT"]` -> `<name>_RNA.cloupe` and `<name>_ADT.cloupe`, both sharing the same UMAPs and cluster/annotation columns.
- **output_name** - Base name for the files; defaults to `metadata.project_name`.
- **executable_path** - Optional path to the `louper` binary; defaults to loupeR's cached install.

**Requirements:** the `loupeR` package and the `louper` executable (installed in the Docker image via `loupeR::setup()`). The step materializes counts into memory (Loupe cannot read BPCells on-disk matrices), so run it on a machine that still has the `bpcells_cache` and enough RAM (~19 GB for ~800k cells). Outputs land in `paths.output_dir`.

---

## Configuration Scenarios

### Scenario 1: Quick Test with Small Dataset

```json
{
  "execution": {
    "mode": "step-by-step",
    "n_cores": 4
  },
  "processing": {
    "features": {"n_features": 2000},
    "reduction": {"pca_dims_rna": 20, "pca_dims_adt": 8}
  },
  "cell_type_annotation": {
    "singler_workers": 4
  }
}
```

Use this to quickly test your pipeline with fewer genes and dimensions.

### Scenario 2: Production Run on Large Dataset

```json
{
  "execution": {
    "mode": "step-by-step",
    "n_cores": 16
  },
  "processing": {
    "features": {"n_features": 5000},
    "reduction": {"pca_dims_rna": 50, "pca_dims_adt": 15}
  },
  "cell_type_annotation": {
    "singler_workers": 16
  }
}
```

Use this for comprehensive analysis with maximum detail.

### Scenario 3: Constrained Memory Machine

```json
{
  "execution": {
    "mode": "chained",
    "n_cores": 4
  },
  "processing": {
    "features": {"n_features": 2000},
    "reduction": {"pca_dims_rna": 15, "pca_dims_adt": 5}
  },
  "cell_type_annotation": {
    "singler_workers": 4
  }
}
```

---

## Common Configuration Mistakes

**Mistake 1:** Using host machine paths in configs when running with Docker
```json
"rna_path": "/mnt/citeseq/raw/combat.h5ad" // [X] WRONG - host path, invisible inside container
"rna_path": "/raw/combat.h5ad" // [OK] CORRECT - container-internal path
```
Tell Docker where your data lives using `DATA_DIR`, not by editing the config:
```bash
make run CONFIG=my_config.json DATA_DIR=/mnt/citeseq
```

**Mistake 2:** Using relative paths
```json
"rna_path": "data/raw/rna.h5ad" // [X] WRONG
"rna_path": "/raw/rna.h5ad" // [OK] CORRECT (when using Docker)
```

**Mistake 3:** `pca_dims_adt` too high
```json
"pca_dims_adt": 20 // [X] WRONG if you only have 13 antibodies
"pca_dims_adt": 10 // [OK] CORRECT (less than number of ADTs)
```

**Mistake 4:** cluster_column doesn't exist
```json
"cluster_column": "my_clusters" // [X] Will fail if this column doesn't exist
// First check what columns exist:
table(obj@meta.data) // In R, check available columns
```

**Mistake 5:** Wrong reference for your cell type
```json
"reference": "immgen" // [X] Wrong for non-immune cells
"reference": "monaco" // [OK] For immune cells (T cells, B cells, etc.)
```

**Mistake 6:** n_cores higher than available CPUs
```json
"n_cores": 32 // [X] WRONG if your machine has only 8 CPUs
"n_cores": 8 // [OK] Match your actual CPU count
```

---

## Doublet Detection and Removal

Doublets are two cells captured in one droplet and sequenced as though they were a single cell. In CITE-seq they produce impossible phenotypes - a T-B doublet reads CD3+CD19+, a T-monocyte doublet CD3+CD14+ - which no classifier or hierarchy can place correctly. The step runs inside `process_seurat` as QC sub-step 2b.

### Configuration

```json
{
  "processing": {
    "doublets": {
      "detect": true,
      "remove": true,
      "sample_column": null,
      "use_adt": true,
      "adt_quantile": 0.90,
      "min_cells": 50,
      "seed": 42
    }
  }
}
```

| Key | Default | Meaning |
|---|---|---|
| `detect` | `true` | Run detection. Set `false` to skip entirely. |
| `remove` | `true` | Remove detected doublets, as for other QC filters. `false` = flag and retain. |
| `sample_column` | `null` | Column identifying the capture/library. Falls back to `batch.column`, then `orig.ident`. |
| `use_adt` | `true` | Use the ADT cross-lineage signature as confirmation (CITE-seq only). |
| `adt_quantile` | `0.90` | Per-marker upper-tail quantile defining ADT positivity. Raise to be stricter. |
| `min_cells` | `50` | Captures smaller than this are skipped and logged. |
| `seed` | `42` | Random seed for `scDblFinder`, for reproducible calls. |

### The two signals

| Signal | Method | Catches | Blind to |
|---|---|---|---|
| **RNA** | `scDblFinder`, per capture | Heterotypic doublets - two cells of different types in one droplet | Homotypic doublets (two cells of the same type), an information limit shared by all expression-based methods |
| **ADT** (confirmatory) | Cross-lineage co-expression | Multiplets that are unambiguous in surface phenotype (CD3+CD19+, CD3+CD14+) | Anything without a mutually exclusive marker pair in the panel |

They are complementary rather than redundant, which is why both run when ADT is present. `doublet_method` records which fired for each cell (`RNA`, `ADT`, or `RNA+ADT`).

### Why per capture

Doublets form within a single droplet emulsion. Running a detector across a merged multi-sample object produces systematically wrong calls, so each capture is processed separately. On disk-backed (BPCells) data each capture's counts are materialised one at a time and released immediately, so peak memory is bounded by the **largest single capture**, not the whole dataset.

### Outputs

Cell-level metadata: `doublet_score` (scDblFinder score), `doublet_class_rna`, `adt_crosslineage` (logical), `doublet_class` (consensus `doublet`/`singlet`), `doublet_method` (which signals fired).

Files and summaries: `<output_prefix>_doublet_per_sample.csv` with per-capture cells, doublet count, percentage and whether the capture was skipped; `obj@misc$doublet_qc` with the same table plus totals and whether removal was applied. The doublet rate per sample is also written to the run log - an unusually high rate points to overloading at capture.

### Notes and limits

The ADT check is deliberately conservative: a cell must sit in the upper tail of **two** mutually exclusive lineage groups (T/B/NK/myeloid/platelet) to be flagged, so it confirms obvious multiplets rather than competing with the transcriptome detector.

Marker names are matched as **whole tokens** after punctuation is normalised, so vendor decoration resolves automatically (`CD3`, `Ab.CD3`, `Hu-CD3`, `CD3-TotalSeqC`, `anti-human CD3` all match `CD3`) while digit-extended markers stay distinct (`CD3` never matches `CD34`/`CD38`, `CD14` never matches `CD141`, `CD20` never matches `CD200`).

Whole-token matching rather than substring matching is deliberate, and was validated against a real dataset: across the 32,367 gene symbols of a CELLxGENE h5ad, a substring rule for `CD3` also matches the unrelated genes `ABCD3`, `HACD3`, `NUDCD3`, `PLCD3` and `PTCD3`. For the same reason the lineage groups list only CD/antibody nomenclature - gene-symbol aliases such as `NCAM1` and `MS4A1` are excluded, because no antibody panel uses them and `NCAM1` collides with the antisense lncRNA `NCAM1-AS1`. One consequence to be aware of: a marker glued to a prefix with no separator (`AbCD3`) will not match, because it cannot be distinguished from the gene `ABCD3` - use a separator or the bare marker name.

If `scDblFinder` is unavailable the RNA signal is skipped with a logged warning and the pipeline continues. If detection would leave fewer than two cells, removal is refused and logged rather than emptying the object.

## Batch Effect Detection and Correction

The `process_seurat` step includes automatic batch effect detection and optional correction. This is controlled by the `batch` section in your config file.

### Configuration

```json
{
  "batch": {
    "column": "orig.ident",
    "detect": true,
    "correct": "auto",
    "correct_rna": "auto",
    "correct_adt": "auto",
    "harmony_theta": 2,
    "harmony_niter": 10,
    "lisi_threshold": 0.7,
    "silhouette_threshold": 0.15,
    "silhouette_n_cells": 5000,
    "wasserstein_threshold": 0.3,
    "min_cells_per_batch": 50,
    "adt_correction_method": "harmony",
    "adt_kbet_threshold": 0.6,
    "adt_guardrail_min_improvement": 0.05
  }
}
```

### How It Works

**Detection** (runs automatically when `detect: true`):

- **RNA batch detection** - computed on the PCA embedding:
  - **Silhouette width** is the RNA batch-effect metric (see below). (iLISI/`lisi` has been removed - the pipeline no longer depends on the `lisi` package; the `lisi_threshold` config key is retained but unused.)
  - **Silhouette width**: measures how well cells cluster by batch vs biology. High values = batch-driven clustering. Subsampled to `silhouette_n_cells` cells (O(n2) memory).

- **ADT batch detection** - two complementary views:
  - **Wasserstein distance** (or KS statistic fallback), computed per-protein on CLR-normalized data: measures distribution shift between batches for each protein. Flagged proteins are logged, and the full per-protein table is written to `<output>_adt_per_protein_distance.csv` (plus `_adt_pairwise_distance.csv`) in `paths.output_dir`.
  - **kBET-style rejection rate** on a provisional ADT PCA: measures how poorly batches mix in the joint embedding - the view that predicts whether WNN clustering will split by batch.
  - Correction is recommended when >10% of proteins exceed `wasserstein_threshold` **or** the kBET rejection rate exceeds `adt_kbet_threshold` (default 0.6).

**Correction** (controlled by `correct`, or per-modality by `correct_rna` / `correct_adt`):

| Setting | Behavior |
|---|---|
| `"auto"` | Correct only when detection flags a batch effect |
| `true` | Always correct, regardless of detection results |
| `false` | Never correct |

RNA and ADT correction are **independent**. `correct` sets both; `correct_rna` and
`correct_adt` override it per modality (each takes `"auto"` / `true` / `false`). So
you can correct **RNA only** (`"correct_rna": true, "correct_adt": false`), **ADT
only** (`"correct_rna": false, "correct_adt": true`), both, or neither. When a
per-modality key is omitted it falls back to `correct`. (String booleans like
`"true"`/`"false"` are also accepted, so a quoted value no longer silently
disables correction.)

- **RNA**: Harmony integration on the PCA embedding -> stored as the `"harmony"` reduction (full mode). In **sketch mode**, Harmony corrects the sketch PCA and the corrected embedding is projected to all cells as `"harmony.full"`.
- **ADT**: `adt_correction_method` chooses HOW correction is applied:
  - `"harmony"` (default) - Harmony corrects the `apca` embedding into `"harmony.adt"`, which WNN/UMAP then consume. The ADT marker VALUES (CLR `"data"` layer) are never touched.
  - `"adtnorm"` - degenerate markers (constant within a batch, all-zero, or low-complexity) are removed, then ADTnorm runs on the raw counts - whole-panel first, then **per-marker** on failure so a single unmodellable marker is dropped rather than aborting correction. Corrected values replace the `"data"` layer.
  - `"none"` - ADT correction is skipped even when detection recommends it.
  - Either way, a **guardrail** re-measures batch mixing after correction and only ADOPTS it if the kBET rejection rate drops by at least `adt_guardrail_min_improvement` (default 0.05); otherwise the correction is reverted and the object keeps the uncorrected ADT. The outcome is recorded in `misc$batch_qc` (`adt_correction_method`, `adt_guardrail`, `adt_correction_applied`).

**Post-correction validation**: After Harmony correction, detection metrics are re-run on the corrected embedding to verify improvement. Results stored in `batch_qc$rna_post_correction`.

### Accessing Results

```r
obj <- readRDS("results/process_seurat_YYYYMMDD_HHMMSS.rds")

# Full batch QC record
obj@misc$batch_qc

# Pre-correction metrics
obj@misc$batch_qc$rna$mean_silhouette # RNA batch metric (LISI decoupled; median_lisi is NA)
obj@misc$batch_qc$adt$flagged_proteins

# Post-correction metrics (if correction was applied)
obj@misc$batch_qc$rna_post_correction$mean_silhouette
```

---

## Orchestrator Architecture (Dispatch-Loop)

The orchestrator uses a **dispatch-loop architecture** rather than hardcoded sequential `if` blocks. Each pipeline step is a self-contained function registered in a dispatch table:

```r
step_dispatch <- list(
  load_data = run_load_data,
  add_metadata = run_add_metadata,
  process_seurat = run_process_seurat,
  cell_type_annotation = run_cell_type_annotation,
  aggregate_ml = run_aggregate_ml
)
```

The main loop iterates over `config$execution$steps` **in the order you specify**, so `add_metadata` can appear at any position after `load_data`:

```json
// Traditional order
"steps": ["load_data", "add_metadata", "process_seurat", "cell_type_annotation", "aggregate_ml"]

// Metadata after annotation
"steps": ["load_data", "process_seurat", "cell_type_annotation", "add_metadata", "aggregate_ml"]

// Skip metadata entirely
"steps": ["load_data", "process_seurat", "cell_type_annotation", "aggregate_ml"]
```

When resuming with `--resume <step>`, the orchestrator dynamically loads the most recent prior step's output by searching backwards through the actual step order - not a hardcoded list.

## Summary

```
1. git clone <repo> && cd cite-seq/alternate_workflow
2. bash setup.sh && bash verify.sh
3. Edit my_config.json (set your absolute paths)
4. conda activate tcell-env
5. Rscript orchestrator.R --config my_config.json
6. Check results/ and logs/ folders
```

**Total time:** 5-10 minutes setup + pipeline runtime (30 min to 5+ hours)
