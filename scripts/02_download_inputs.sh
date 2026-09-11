#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
RAW="$ROOT/data/raw"
LOG="$ROOT/logs/02_download_inputs.log"
mkdir -p "$RAW" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

download_plain_and_gzip() {
  local url="$1"
  local out="$2"
  local expected_md5="$3"
  local tmp="${out%.gz}"

  if [[ -s "$out" ]]; then
    echo "SKIP existing $out"
    return
  fi

  curl -fL --retry 5 --retry-delay 5 --continue-at - "$url" -o "$tmp"
  echo "$expected_md5  $tmp" | md5sum -c -
  gzip -n -f "$tmp"
}

download_gzip() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "SKIP existing $out"
    return
  fi
  curl -fL --retry 5 --retry-delay 5 --continue-at - "$url" -o "$out"
  gzip -t "$out"
}

download_plain_and_gzip \
  https://ndownloader.figshare.com/files/21356610 \
  "$RAW/UKB_LifetimeMDD.covararraypcs.assoc.logistic.gz" \
  ccb58b7a3d3b2a32d0ac651fa4282bbb &
download_plain_and_gzip \
  https://ndownloader.figshare.com/files/21356598 \
  "$RAW/UKB_MDDRecur.covararraypcs.assoc.logistic.gz" \
  2f6730c7e95583691b8e796ec0df2830 &
download_plain_and_gzip \
  https://ndownloader.figshare.com/files/21356688 \
  "$RAW/UKB_GPpsy.covararraypcs.assoc.logistic.gz" \
  d1677647fa42f4f59f27b2aa296a141e &
download_plain_and_gzip \
  https://ndownloader.figshare.com/files/21356586 \
  "$RAW/UKB_ICD10Dep.covararraypcs.assoc.logistic.gz" \
  1337d06b801074403d2fc1b74ed53a7f &

download_gzip \
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/harmonised/28067908-GCST004132-EFO_0000384.h.tsv.gz \
  "$RAW/CD_deLange2017.h.tsv.gz" &
download_gzip \
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz \
  "$RAW/UC_deLange2017.h.tsv.gz" &

wait

ln -sfn /root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz "$RAW/MDD_Howard2019.txt.gz"
ln -sfn /root/IBD/19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz "$RAW/DEP_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/deLange2017/deLange2017_harmonised.tsv.gz "$RAW/IBD_deLange2017.h.tsv.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/IBD_STRICT_FinnGen_R12.gz "$RAW/IBD_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/CD_STRICT2_FinnGen_R12.gz "$RAW/CD_FinnGen_R12.gz"
ln -sfn /root/IBD/00_RawData/GWAS/IBD/UC_STRICT2_FinnGen_R12.gz "$RAW/UC_FinnGen_R12.gz"

find "$RAW" -maxdepth 1 -type f -o -type l | sort | xargs -r sha256sum > "$ROOT/provenance/raw_input_sha256.txt"
date --iso-8601=seconds > "$ROOT/provenance/download_completed_at.txt"
