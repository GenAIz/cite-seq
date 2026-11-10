#!/bin/bash
set -e

CONFIG_FILE="${FSP_CONFIG:-config.json}"
CONFIG_PATH="${SF_INPUT_PATH:-run/in}/configs/${CONFIG_FILE}"

echo "============================================================"
echo "  CITE-seq Functional State Profile"
echo "============================================================"
echo "Config : ${CONFIG_PATH}"
echo "Input  : ${SF_INPUT_PATH:-run/in}"
echo "Output : ${SF_OUTPUT_PATH:-run/out}"
echo ""

exec Rscript /app/main.R --config "${CONFIG_PATH}"
