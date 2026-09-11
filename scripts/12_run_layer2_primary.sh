#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
LOG_DIR="$ROOT/logs"
mkdir -p "$LOG_DIR" "$ROOT/results/local"

export LAVA_WORKERS=6

Rscript "$ROOT/scripts/10_run_lava_primary_pair.R" \
  HowardMDD__FinnGenIBD \
  > "$LOG_DIR/10_lava_HowardMDD__FinnGenIBD.log" 2>&1 &
pid1=$!

Rscript "$ROOT/scripts/10_run_lava_primary_pair.R" \
  FinnGenDEP__deLangeIBD \
  > "$LOG_DIR/10_lava_FinnGenDEP__deLangeIBD.log" 2>&1 &
pid2=$!

printf "HowardMDD__FinnGenIBD PID=%s\n" "$pid1"
printf "FinnGenDEP__deLangeIBD PID=%s\n" "$pid2"

status=0
wait "$pid1" || status=1
wait "$pid2" || status=1
if [[ "$status" -ne 0 ]]; then
  printf "At least one LAVA pair failed. Inspect logs under %s\n" "$LOG_DIR" >&2
  exit 1
fi

python3 "$ROOT/scripts/11_summarize_lava_reproducibility.py"
python3 "$ROOT/scripts/12_local_power_diagnostic.py"
Rscript "$ROOT/scripts/13_verify_corrected_local_results.R"

sha256sum \
  "$ROOT/scripts/10_run_lava_primary_pair.R" \
  "$ROOT/scripts/11_summarize_lava_reproducibility.py" \
  "$ROOT/scripts/12_local_power_diagnostic.py" \
  "$ROOT/scripts/13_verify_corrected_local_results.R" \
  "$ROOT/scripts/12_run_layer2_primary.sh" \
  >> "$ROOT/provenance/code_and_protocol_sha256.txt"

date --iso-8601=seconds > "$ROOT/provenance/layer2_primary_completed_at.txt"
