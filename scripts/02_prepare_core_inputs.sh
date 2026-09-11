#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
RAW="$ROOT/data/raw"
mkdir -p "$RAW" "$ROOT/logs" "$ROOT/provenance"

ln -sfn /root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz "$RAW/MDD_Howard2019.txt.gz"
ln -sfn /root/IBD/19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz "$RAW/DEP_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/deLange2017/deLange2017_harmonised.tsv.gz "$RAW/IBD_deLange2017.h.tsv.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/IBD_STRICT_FinnGen_R12.gz "$RAW/IBD_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/CD_STRICT2_FinnGen_R12.gz "$RAW/CD_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/UC_STRICT2_FinnGen_R12.gz "$RAW/UC_FinnGen_R12.gz"

find "$RAW" -maxdepth 1 -type l -print0 | sort -z | xargs -0 -r sha256sum \
  > "$ROOT/provenance/core_raw_input_sha256.txt"
