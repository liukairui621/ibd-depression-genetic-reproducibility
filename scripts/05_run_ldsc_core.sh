#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
M="$ROOT/data/munged_core"
OUT="$ROOT/results/global_core"
LOG="$ROOT/logs/05_run_ldsc_core.log"
PY=/root/anaconda3/envs/ldsc/bin/python
LDSC=/root/ldsc/ldsc.py
LD=/root/IBD/00_RawData/Reference/LD_EUR/LDscore/LDscore.
W=/root/IBD/00_RawData/Reference/LD_EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
mkdir -p "$OUT/h2" "$OUT/rg" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

DEP_TRAITS=(MDD_Howard2019 DEP_FinnGen_R12)
IBD_TRAITS=(IBD_deLange2017 IBD_FinnGen_R12 CD_FinnGen_R12 UC_FinnGen_R12)

run_h2() {
  local t="$1"
  "$PY" "$LDSC" --h2 "$M/$t.sumstats.gz" \
    --ref-ld-chr "$LD" --w-ld-chr "$W" --out "$OUT/h2/$t"
}
export PY LDSC M LD W OUT
export -f run_h2
printf '%s\n' "${DEP_TRAITS[@]}" "${IBD_TRAITS[@]}" | sort -u \
  | xargs -n1 -P4 bash -c 'run_h2 "$0"'

PAIR_FILE="$ROOT/provenance/layer1_core_pair_manifest.tsv"
printf 'pair_class\ttrait1\ttrait2\n' > "$PAIR_FILE"
for d in "${DEP_TRAITS[@]}"; do
  for i in "${IBD_TRAITS[@]}"; do
    printf 'cross_disease\t%s\t%s\n' "$d" "$i" >> "$PAIR_FILE"
  done
done
printf 'within_depression\tMDD_Howard2019\tDEP_FinnGen_R12\n' >> "$PAIR_FILE"
for ((a=0; a<${#IBD_TRAITS[@]}; a++)); do
  for ((b=a+1; b<${#IBD_TRAITS[@]}; b++)); do
    printf 'within_IBD\t%s\t%s\n' "${IBD_TRAITS[$a]}" "${IBD_TRAITS[$b]}" >> "$PAIR_FILE"
  done
done

run_rg() {
  local cls="$1"
  local t1="$2"
  local t2="$3"
  "$PY" "$LDSC" --rg "$M/$t1.sumstats.gz,$M/$t2.sumstats.gz" \
    --ref-ld-chr "$LD" --w-ld-chr "$W" \
    --out "$OUT/rg/${cls}__${t1}__${t2}"
}
export -f run_rg
tail -n +2 "$PAIR_FILE" | tr '\t' ' ' \
  | xargs -n3 -P4 bash -c 'run_rg "$0" "$1" "$2"'

date --iso-8601=seconds > "$ROOT/provenance/ldsc_core_completed_at.txt"
