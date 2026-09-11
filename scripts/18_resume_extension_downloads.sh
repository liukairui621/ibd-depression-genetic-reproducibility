#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
PART="$ROOT/data/partial_downloads"
RAW="$ROOT/data/raw"
LOG="$ROOT/logs/18_resume_extension_downloads.log"
mkdir -p "$PART" "$RAW" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

download_plain() {
  local file="$1"
  local file_id="$2"
  local remote_file="$3"
  local expected_size="$4"
  local expected_md5="$5"
  local url="https://pfigshare-u-files.s3-eu-west-1.amazonaws.com/$file_id/$remote_file"

  if [[ ! -f "$PART/$file" ]] ||
     [[ "$(stat -c %s "$PART/$file")" != "$expected_size" ]] ||
     ! echo "$expected_md5  $PART/$file" | md5sum -c -; then
    aria2c \
      --continue=true \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --content-disposition=false \
      --file-allocation=none \
      --max-connection-per-server=16 \
      --split=16 \
      --min-split-size=1M \
      --retry-wait=5 \
      --max-tries=0 \
      --connect-timeout=30 \
      --timeout=60 \
      --dir="$PART" \
      --out="$file" \
      "$url"
  fi

  [[ "$(stat -c %s "$PART/$file")" == "$expected_size" ]]
  echo "$expected_md5  $PART/$file" | md5sum -c -
  gzip -n -c "$PART/$file" > "$RAW/$file.gz"
  gzip -t "$RAW/$file.gz"
}

download_gzip() {
  local file="$1"
  local url="$2"
  local expected_size="$3"

  if [[ ! -f "$PART/$file" ]] ||
     [[ "$(stat -c %s "$PART/$file")" != "$expected_size" ]] ||
     ! gzip -t "$PART/$file"; then
    aria2c \
      --continue=true \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --content-disposition=false \
      --file-allocation=none \
      --max-connection-per-server=8 \
      --split=8 \
      --min-split-size=1M \
      --retry-wait=5 \
      --max-tries=0 \
      --connect-timeout=30 \
      --timeout=60 \
      --dir="$PART" \
      --out="$file" \
      "$url"
  fi

  [[ "$(stat -c %s "$PART/$file")" == "$expected_size" ]]
  gzip -t "$PART/$file"
  cp -f "$PART/$file" "$RAW/$file"
}

download_plain \
  UKB_LifetimeMDD.covararraypcs.assoc.logistic \
  21356610 \
  LifetimeMDD.covararraypcs.assoc.logistic \
  763812382 \
  ccb58b7a3d3b2a32d0ac651fa4282bbb &
pid1=$!
download_plain \
  UKB_MDDRecur.covararraypcs.assoc.logistic \
  21356598 \
  MDDRecur.covararraypcs.assoc.logistic \
  763547170 \
  2f6730c7e95583691b8e796ec0df2830 &
pid2=$!
download_plain \
  UKB_GPpsy.covararraypcs.assoc.logistic \
  21356688 \
  GPpsy.covararraypcs.assoc.logistic \
  777948024 \
  d1677647fa42f4f59f27b2aa296a141e &
pid3=$!
download_plain \
  UKB_ICD10Dep.covararraypcs.assoc.logistic \
  21356586 \
  ICD10Dep.covararraypcs.assoc.logistic \
  772599146 \
  1337d06b801074403d2fc1b74ed53a7f &
pid4=$!
download_gzip \
  CD_deLange2017.h.tsv.gz \
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/harmonised/28067908-GCST004132-EFO_0000384.h.tsv.gz \
  296840782 &
pid5=$!
download_gzip \
  UC_deLange2017.h.tsv.gz \
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz \
  299230653 &
pid6=$!

status=0
for pid in "$pid1" "$pid2" "$pid3" "$pid4" "$pid5" "$pid6"; do
  wait "$pid" || status=1
done
[[ "$status" == 0 ]]

{
  printf "file\tbytes\tmd5\tsha256\n"
  for file in \
    "$RAW/UKB_LifetimeMDD.covararraypcs.assoc.logistic.gz" \
    "$RAW/UKB_MDDRecur.covararraypcs.assoc.logistic.gz" \
    "$RAW/UKB_GPpsy.covararraypcs.assoc.logistic.gz" \
    "$RAW/UKB_ICD10Dep.covararraypcs.assoc.logistic.gz" \
    "$RAW/CD_deLange2017.h.tsv.gz" \
    "$RAW/UC_deLange2017.h.tsv.gz"; do
    printf "%s\t%s\t%s\t%s\n" \
      "$(basename "$file")" \
      "$(stat -c %s "$file")" \
      "$(md5sum "$file" | cut -d' ' -f1)" \
      "$(sha256sum "$file" | cut -d' ' -f1)"
  done
} > "$ROOT/provenance/extension_download_checksums.tsv"

date --iso-8601=seconds > "$ROOT/provenance/extension_download_completed_at.txt"
