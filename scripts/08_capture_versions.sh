#!/usr/bin/env bash
set -euo pipefail
ROOT=/root/IBD/20_Reproducibility_Ladder
OUT="$ROOT/provenance/software_versions.txt"
{
  date --iso-8601=seconds
  hostname
  uname -a
  /root/anaconda3/envs/ldsc/bin/python --version
  /root/anaconda3/envs/ldsc/bin/python - <<'PY'
import numpy, pandas, scipy
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
PY
  python3 --version
  R --version | head -n 1
  md5sum /root/ldsc/ldsc.py /root/ldsc/munge_sumstats.py
} > "$OUT" 2>&1

find "$ROOT/scripts" "$ROOT/provenance" -type f -print0 \
  | sort -z | xargs -0 sha256sum > "$ROOT/provenance/code_and_protocol_sha256.txt"
