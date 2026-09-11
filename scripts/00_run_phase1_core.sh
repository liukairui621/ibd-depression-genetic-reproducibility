#!/usr/bin/env bash
set -euo pipefail
ROOT=/root/IBD/20_Reproducibility_Ladder
LOG="$ROOT/logs/00_run_phase1_core.log"
mkdir -p "$ROOT/logs"
exec > >(tee "$LOG") 2>&1

python3 "$ROOT/scripts/01_inventory_inputs.py"
bash "$ROOT/scripts/02_prepare_core_inputs.sh"
bash "$ROOT/scripts/04_munge_core.sh"
bash "$ROOT/scripts/05_run_ldsc_core.sh"
python3 "$ROOT/scripts/06_parse_ldsc_core.py"
python3 "$ROOT/scripts/07_summarize_layer1_core.py"
bash "$ROOT/scripts/08_capture_versions.sh"

date --iso-8601=seconds > "$ROOT/provenance/phase1_core_completed_at.txt"
