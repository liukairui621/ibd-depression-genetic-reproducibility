#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
RAW="$ROOT/data/raw"
MUNGED="$ROOT/data/munged"
LOG="$ROOT/logs/19_prepare_extension.log"
mkdir -p "$RAW" "$MUNGED" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

required=(
  UKB_LifetimeMDD.covararraypcs.assoc.logistic.gz
  UKB_MDDRecur.covararraypcs.assoc.logistic.gz
  UKB_GPpsy.covararraypcs.assoc.logistic.gz
  UKB_ICD10Dep.covararraypcs.assoc.logistic.gz
  CD_deLange2017.h.tsv.gz
  UC_deLange2017.h.tsv.gz
)
for file in "${required[@]}"; do
  [[ -s "$RAW/$file" ]]
  gzip -t "$RAW/$file"
done

ln -sfn \
  /root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz \
  "$RAW/MDD_Howard2019.txt.gz"
ln -sfn \
  /root/IBD/19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz \
  "$RAW/DEP_FinnGen_R12.gz"
ln -sfn \
  /root/IBD/00_RawData/GWAS/IBD/deLange2017/deLange2017_harmonised.tsv.gz \
  "$RAW/IBD_deLange2017.h.tsv.gz"
ln -sfn \
  /root/IBD/00_RawData/GWAS/IBD/IBD_STRICT_FinnGen_R12.gz \
  "$RAW/IBD_FinnGen_R12.gz"
ln -sfn \
  /root/IBD/00_RawData/GWAS/IBD/CD_STRICT2_FinnGen_R12.gz \
  "$RAW/CD_FinnGen_R12.gz"
ln -sfn \
  /root/IBD/00_RawData/GWAS/IBD/UC_STRICT2_FinnGen_R12.gz \
  "$RAW/UC_FinnGen_R12.gz"

traits=(
  MDD_Howard2019
  DEP_FinnGen_R12
  UKB_LifetimeMDD
  UKB_MDDRecur
  UKB_GPpsy
  UKB_ICD10Dep
  IBD_deLange2017
  CD_deLange2017
  UC_deLange2017
  IBD_FinnGen_R12
  CD_FinnGen_R12
  UC_FinnGen_R12
)

need_munge=0
for trait in "${traits[@]}"; do
  file="$MUNGED/$trait.sumstats.gz"
  if [[ ! -s "$file" ]] || ! /root/anaconda3/bin/python -c \
    'import gzip,sys
header=gzip.open(sys.argv[1], "rt").readline().rstrip("\n").split("\t")
raise SystemExit(0 if set(header)=={"SNP","A1","A2","Z","N"} else 1)' \
    "$file"; then
    need_munge=1
    break
  fi
done

if [[ "$need_munge" == 1 ]]; then
  /root/anaconda3/bin/python \
    "$ROOT/scripts/03_normalize_delange.py"
  bash "$ROOT/scripts/04_munge_sumstats.sh"
else
  echo "All 12 munged files passed the required-column check; reusing them."
fi

manifest="$ROOT/provenance/extension_munged_manifest.tsv"
printf "trait\trows_including_header\tbytes\tsha256\n" > "$manifest"
for trait in "${traits[@]}"; do
  file="$MUNGED/$trait.sumstats.gz"
  [[ -s "$file" ]]
  gzip -t "$file"
  header=$(
    /root/anaconda3/bin/python -c \
      'import gzip,sys; print(gzip.open(sys.argv[1], "rt").readline().rstrip("\n"))' \
      "$file"
  )
  HEADER="$header" /root/anaconda3/bin/python -c \
    'import os
header=os.environ["HEADER"].split("\t")
raise SystemExit(0 if set(header)=={"SNP","A1","A2","Z","N"} else 1)'
  rows=$(gzip -cd "$file" | wc -l)
  printf "%s\t%s\t%s\t%s\n" \
    "$trait" \
    "$rows" \
    "$(stat -Lc %s "$file")" \
    "$(sha256sum "$file" | cut -d' ' -f1)" \
    >> "$manifest"
done

date --iso-8601=seconds > "$ROOT/provenance/extension_munging_completed_at.txt"
