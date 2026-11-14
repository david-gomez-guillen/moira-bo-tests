#!/usr/bin/env bash
# Runs the standalone MOIRA-BO on a problem of this repository, creating the
# conda environment of environment.yml the first time.
#
#   ./run_moira_bo.sh                           1-bimodal, to 56 evaluations
#   ./run_moira_bo.sh 7-valley --seed 3         another problem and seed
#
# The problem is a folder of the repository, with or without its number, and
# every option after it goes to moira-bo/run.R (see moira-bo/README.md).
set -euo pipefail
cd "$(dirname "$0")"

ENV_NAME=moira-bo
CONDA=${CONDA_EXE:-$(command -v conda || true)}
if [ -z "$CONDA" ]; then
  echo "conda was not found: install Miniforge or Miniconda first" >&2
  exit 1
fi

if ! "$CONDA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Creating the conda environment '$ENV_NAME', which takes a few minutes the first time"
  "$CONDA" env create -f environment.yml
fi

if [ $# -eq 0 ] || [[ $1 == -* ]]; then
  set -- 1-bimodal "$@"
fi
exec "$CONDA" run -n "$ENV_NAME" --no-capture-output Rscript moira-bo/run.R "$@"
