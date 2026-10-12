#!/usr/bin/env bash
# One-time setup (macOS / Linux): Python environment, R packages, then an install check.
#   ./setup.sh
# R itself must already be installed (https://cran.r-project.org).
set -euo pipefail
cd "$(dirname "$0")"

VENV_DIR="${VENV_DIR:-venv}"

# Newest Python 3.11+ available
PYTHON=""
for candidate in python3.14 python3.13 python3.12 python3.11 python3; do
  if command -v "$candidate" >/dev/null 2>&1 &&
     "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
    PYTHON="$candidate"
    break
  fi
done
if [ -z "$PYTHON" ]; then
  echo "Python 3.11 or newer is required. Install it (e.g. from https://www.python.org) and rerun ./setup.sh"
  exit 1
fi

if [ ! -x "$VENV_DIR/bin/python" ]; then
  echo "Creating Python environment in $VENV_DIR/ with $("$PYTHON" --version)..."
  "$PYTHON" -m venv "$VENV_DIR"
fi

echo "Installing Python packages..."
"$VENV_DIR/bin/python" -m pip install --quiet --upgrade pip
"$VENV_DIR/bin/python" -m pip install --quiet -r requirements.txt

if command -v Rscript >/dev/null 2>&1; then
  echo "Installing R packages (only missing ones; DESeq2 can take a while the first time)..."
  Rscript R/install_packages.R
else
  echo "Rscript not found: install R from https://cran.r-project.org, then rerun ./setup.sh"
fi

echo
"$VENV_DIR/bin/python" run_pipeline.py --check
echo
echo "Activate the environment with:  source $VENV_DIR/bin/activate"
