#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
RAW=/root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz
OUT="$ROOT/data/munged_core/MDD_Howard2019"
PY=/root/anaconda3/envs/ldsc/bin/python
PY3=/root/anaconda3/bin/python
MUNGE=/root/ldsc/munge_sumstats.py
ALLELES=/root/IBD/00_RawData/Reference/LD_EUR/w_hm3_alleles.snplist
LOG="$ROOT/logs/01_munge_howard_no23andme.log"

mkdir -p "$ROOT/data/munged_core" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

if [[ -e "$OUT.sumstats.gz" || -e "$OUT.log" ]]; then
  echo "Refusing to overwrite an existing corrected munged file" >&2
  exit 1
fi

echo "Started: $(date -Iseconds)"
sha256sum "$RAW" "$ALLELES"

"$PY" "$MUNGE" \
  --sumstats "$RAW" \
  --snp MarkerName \
  --a1 A1 \
  --a2 A2 \
  --frq Freq \
  --p P \
  --signed-sumstats LogOR,0 \
  --N-cas 170756 \
  --N-con 329443 \
  --merge-alleles "$ALLELES" \
  --chunksize 500000 \
  --out "$OUT"

"$PY3" - "$OUT.sumstats.gz" <<'PY'
import gzip
import math
import sys

path = sys.argv[1]
rows = 0
n_values = set()
with gzip.open(path, "rt") as handle:
    header = handle.readline().split()
    n_idx = header.index("N")
    for line in handle:
        fields = line.split()
        if not fields:
            continue
        rows += 1
        n_values.add(float(fields[n_idx]))

expected_n = 500199.0
if n_values != {expected_n}:
    raise SystemExit(f"Unexpected N values: {sorted(n_values)[:10]}")
if rows == 0:
    raise SystemExit("Corrected munged file contains no variants")
print(f"Verified variants: {rows}")
print(f"Verified constant N: {expected_n:.0f}")
PY

sha256sum "$OUT.sumstats.gz" "$OUT.log" > "$ROOT/provenance/HOWARD_NO23ANDME_MUNGED_SHA256SUMS.txt"
echo "Completed: $(date -Iseconds)"
