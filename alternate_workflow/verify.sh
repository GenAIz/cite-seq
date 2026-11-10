#!/bin/bash
#
# Workflow verification: checks that the project is set up correctly.
# Run from the alternate_workflow root (the directory containing orchestrator.R).
#
# Usage: bash verify.sh

echo "Workflow verification"
echo "====================="
echo ""

checks_passed=0
checks_failed=0

# 1. orchestrator.R
if [ -f "orchestrator.R" ]; then
  echo "[PASS] orchestrator.R found"
  ((checks_passed++))
else
  echo "[FAIL] orchestrator.R not found"
  ((checks_failed++))
fi

# 2. script/ folder
if [ -d "script" ]; then
  echo "[PASS] script/ folder exists"
  ((checks_passed++))
else
  echo "[FAIL] script/ folder not found"
  ((checks_failed++))
fi

# 3. R scripts in script/
required_scripts=(
  "logging_module.R"
  "validation_module.R"
  "load_data_updated.R"
  "add_metadata.R"
  "process_seurat_pipeline.R"
  "cell_type_annotation.R"
  "aggregate_obj_ML_processing.R"
  "export_loupe.R"
)
scripts_missing=0
for script in "${required_scripts[@]}"; do
  [ -f "script/$script" ] || { echo "  missing: script/$script"; ((scripts_missing++)); }
done
if [ "$scripts_missing" -eq 0 ]; then
  echo "[PASS] all ${#required_scripts[@]} R scripts present in script/"
  ((checks_passed++))
else
  echo "[FAIL] $scripts_missing R script(s) missing from script/"
  ((checks_failed++))
fi

# 4. python/ and run_stcat.py
if [ -f "python/run_stcat.py" ]; then
  echo "[PASS] python/run_stcat.py found"
  ((checks_passed++))
else
  echo "[FAIL] python/run_stcat.py not found"
  ((checks_failed++))
fi

# 5. configuration files
config_files_found=0
[ -f "config_standard.json" ] && ((config_files_found++))
[ -f "my_config.json" ] && ((config_files_found++))
if [ "$config_files_found" -gt 0 ]; then
  echo "[PASS] configuration file(s) found"
  ((checks_passed++))
else
  echo "[FAIL] no configuration files found"
  ((checks_failed++))
fi

# 6. environment files
env_ok=0
[ -f "environment.yml" ] && ((env_ok++))
[ -f "Dockerfile" ] && ((env_ok++))
if [ "$env_ok" -gt 0 ]; then
  echo "[PASS] environment files found"
  ((checks_passed++))
else
  echo "[WARN] environment files missing (environment.yml / Dockerfile)"
fi

# 7. directory structure
dirs_ok=0
[ -d "data/raw" ] && ((dirs_ok++))
[ -d "results" ] && ((dirs_ok++))
if [ "$dirs_ok" -gt 0 ]; then
  echo "[PASS] project directories present"
  ((checks_passed++))
else
  echo "[INFO] run setup.sh to create data/raw and results"
fi

# 8. R installation
if command -v R &> /dev/null; then
  echo "[PASS] R installed ($(R --version | head -n 1))"
  ((checks_passed++))
else
  echo "[FAIL] R not found in PATH (https://www.r-project.org/)"
  ((checks_failed++))
fi

# 9. JSON validity
if command -v python3 &> /dev/null && [ -f "config_standard.json" ]; then
  if python3 -m json.tool config_standard.json > /dev/null 2>&1; then
    echo "[PASS] config_standard.json is valid JSON"
    ((checks_passed++))
  else
    echo "[FAIL] config_standard.json is not valid JSON"
    ((checks_failed++))
  fi
else
  echo "[INFO] skipping JSON validation (python3 or config_standard.json absent)"
fi

# 10. orchestrator references script/
script_refs=$(grep -c "script/" orchestrator.R 2>/dev/null || echo "0")
if [ "$script_refs" -ge 7 ]; then
  echo "[PASS] orchestrator.R references script/ ($script_refs)"
  ((checks_passed++))
else
  echo "[FAIL] orchestrator.R has only $script_refs script/ references (expected 7+)"
  ((checks_failed++))
fi

echo ""
echo "Results: $checks_passed passed, $checks_failed failed"
echo ""

if [ "$checks_failed" -eq 0 ]; then
  echo "Setup is complete."
  echo "  Rscript orchestrator.R --config my_config.json --dry-run"
  echo "  Rscript orchestrator.R --config my_config.json"
  exit 0
elif [ "$checks_failed" -le 2 ]; then
  echo "Setup is mostly complete. Remaining items:"
  [ -d "data/raw" ] || echo "  mkdir -p data/raw"
  [ -d "results" ] || echo "  mkdir -p results"
  { [ ! -f "my_config.json" ] && [ -f "config_standard.json" ]; } && echo "  cp config_standard.json my_config.json"
  exit 0
else
  echo "Setup needs attention. Run: bash setup.sh"
  exit 1
fi
