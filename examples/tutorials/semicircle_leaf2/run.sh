#!/bin/bash
# Usage:  bash run.sh <results_folder_name>
# Example: bash run.sh pinpd_t112_6x20_single_source
#
# Regenerates input files, runs the simulator, then copies the vtk output and
# the generator script into results/<name>/ so every run is self-documented.

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if [ -z "$1" ]; then
    echo "Usage: bash run.sh <results_folder_name>"
    exit 1
fi

RESULTS="results/$1"

echo "=== Regenerating input files ==="
python3 generate_tissue_files2.py

echo "=== Running simulator ==="
../../../bin/simulator semicircle_leaf2.model semicircle_leaf2.init semicircle_leaf2.solver

echo "=== Saving results to $RESULTS ==="
mkdir -p "$RESULTS"
cp vtk/*.vtu vtk/tissue.pvd "$RESULTS/"
cp generate_tissue_files2.py "$RESULTS/generate_tissue_files2.py.snapshot"
cp semicircle_leaf2.model "$RESULTS/"
cp semicircle_leaf2.solver "$RESULTS/"

echo "=== Done. Results saved to $RESULTS ==="
echo "    Open $RESULTS/tissue.pvd in ParaView."
