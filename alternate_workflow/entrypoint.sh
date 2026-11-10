#!/usr/bin/env bash
set -e

if [ -n "$HOST_UID" ] && [ -n "$HOST_GID" ]; then
  if [ "$(id -u)" = "0" ]; then
    groupadd -g "$HOST_GID" runner 2>/dev/null || true
    useradd  -u "$HOST_UID" -g "$HOST_GID" -m -s /bin/bash runner 2>/dev/null || true
    exec gosu runner "$0" "$@"
  fi
fi

COMMAND=$1
CONFIG=$2
STEP=$3

case "$COMMAND" in

  run)
    echo "Running full pipeline..."
    Rscript /app/orchestrator.R --config /app/$CONFIG
    ;;

  resume)
    echo "Resuming pipeline from step: $STEP"
    Rscript /app/orchestrator.R --config /app/$CONFIG --resume $STEP
    ;;

  validate)
    echo "Validating configuration..."
    Rscript /app/orchestrator.R --config /app/$CONFIG --validate
    ;;

  dry-run)
    echo "Checking pipeline without execution..."
    Rscript /app/orchestrator.R --config /app/$CONFIG --dry-run
    ;;

  shell)
    echo "Opening shell..."
    exec /bin/bash
    ;;

  *)
    echo "Usage:"
    echo "  run <config>"
    echo "  resume <config> <step>"
    echo "  validate <config>"
    echo "  dry-run <config>"
    echo "  shell"
    exit 1
    ;;

esac
