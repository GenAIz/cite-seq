#!/bin/bash
#
# Workflow setup: creates the directory layout and organizes source files.
# Run from the alternate_workflow root (the directory containing orchestrator.R).
#
# Usage: bash setup.sh

set -e

echo "Workflow setup"
echo "=============="
echo ""

if [ ! -f "orchestrator.R" ]; then
  echo "ERROR: orchestrator.R not found in the current directory."
  echo "Run this script from the alternate_workflow root."
  exit 1
fi

echo "Creating directories..."
for d in script python data/raw results logs; do
  if [ -d "$d" ]; then
    echo "  exists: $d/"
  else
    mkdir -p "$d"
    echo "  created: $d/"
  fi
done
echo ""

echo "Organizing Python scripts..."
if [ -f "run_stcat.py" ] && [ ! -f "python/run_stcat.py" ]; then
  mv "run_stcat.py" "python/"
  echo "  moved: run_stcat.py -> python/"
fi

echo "Organizing R scripts..."
SCRIPTS=(
  "logging_module.R"
  "validation_module.R"
  "load_data_updated.R"
  "add_metadata.R"
  "process_seurat_pipeline.R"
  "cell_type_annotation.R"
  "aggregate_obj_ML_processing.R"
)
moved_count=0
for script in "${SCRIPTS[@]}"; do
  if [ -f "$script" ] && [ ! -f "script/$script" ]; then
    mv "$script" "script/"
    echo "  moved: $script -> script/"
    ((moved_count++))
  fi
done
echo ""

echo "Preparing configuration..."
if [ -f "config_standard.json" ] && [ ! -f "my_config.json" ]; then
  cp config_standard.json my_config.json
  echo "  created my_config.json from config_standard.json (edit before running)"
elif [ -f "my_config.json" ]; then
  echo "  my_config.json already present"
else
  echo "  config_standard.json not found; provide a config before running"
fi
echo ""

echo "Environment:"
if command -v R &> /dev/null; then
  echo "  R:      $(R --version | head -n 1)"
else
  echo "  R:      not found in PATH"
fi
if command -v conda &> /dev/null; then
  echo "  conda:  $(conda --version)"
else
  echo "  conda:  not found in PATH"
fi
if command -v docker &> /dev/null; then
  echo "  docker: $(docker --version)"
else
  echo "  docker: not found in PATH (optional)"
fi
echo ""

echo "Setup complete."
echo "Scripts moved this run: $moved_count"
echo ""
echo "Next steps:"
echo "  1. Edit my_config.json to set your data paths and parameters."
echo "  2. Place input files in data/raw/."
echo "  3. Create the Python environment: conda env create -f environment.yml"
echo "  4. Validate: Rscript orchestrator.R --config my_config.json --dry-run"
echo "  5. Run: Rscript orchestrator.R --config my_config.json"

exit 0
