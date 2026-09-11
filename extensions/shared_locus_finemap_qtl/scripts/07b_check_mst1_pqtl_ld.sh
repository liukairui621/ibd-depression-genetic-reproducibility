#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl}"
PLINK=/root/anaconda3/bin/plink
REF=/root/IBD/00_RawData/Reference/LD_EUR/1000G_EUR_Phase3_plink/1000G.EUR.QC.3
OUT="$ROOT/results/opentargets/ld_checks"
mkdir -p "$OUT" "$ROOT/logs"

"$PLINK" --bfile "$REF" \
  --ld-snp rs11130213 \
  --ld-window 99999 \
  --ld-window-kb 1000 \
  --ld-window-r2 0 \
  --r2 \
  --out "$OUT/MST1_pqtl_leads" \
  > "$ROOT/logs/07b_MST1_pqtl_ld.log" 2>&1

grep -E 'rs111812549|rs113465792|rs2029591|rs6774202|rs34762726|rs1131095|rs9835157' \
  "$OUT/MST1_pqtl_leads.ld" > "$OUT/MST1_selected_pairwise_ld.tsv"

