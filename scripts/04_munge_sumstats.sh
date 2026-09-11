#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
RAW="$ROOT/data/raw"
PREP="$ROOT/data/prepared"
OUT="$ROOT/data/munged"
LOG="$ROOT/logs/04_munge_sumstats.log"
PY=/root/anaconda3/envs/ldsc/bin/python
MUNGE=/root/ldsc/munge_sumstats.py
HM3=/root/IBD/00_RawData/Reference/LD_EUR/w_hm3_alleles.snplist
mkdir -p "$OUT" "$ROOT/logs"
exec > >(tee "$LOG") 2>&1

if [[ ! -s "$ROOT/data/munged_core/MDD_Howard2019.sumstats.gz" ]]; then
  bash "$ROOT/scripts/01_munge_howard_no23andme.sh"
fi
ln -sfn "$ROOT/data/munged_core/MDD_Howard2019.sumstats.gz" "$OUT/MDD_Howard2019.sumstats.gz"
python3 "$ROOT/scripts/01b_validate_howard_no23andme.py"
ln -sfn /root/IBD/16_QualityUpgrade/01_LDSC_deLange_MDD/munged/deLange2017_IBD.sumstats.gz "$OUT/IBD_deLange2017.sumstats.gz"
ln -sfn /root/IBD/03_LDSC/munged/IBD_STRICT_FinnGen.sumstats.gz "$OUT/IBD_FinnGen_R12.sumstats.gz"
ln -sfn /root/IBD/03_LDSC/munged/CD_STRICT2_FinnGen.sumstats.gz "$OUT/CD_FinnGen_R12.sumstats.gz"
ln -sfn /root/IBD/03_LDSC/munged/UC_STRICT2_FinnGen.sumstats.gz "$OUT/UC_FinnGen_R12.sumstats.gz"

munge_delange() {
  local trait="$1"
  local n="$2"
  "$PY" "$MUNGE" \
    --sumstats "$PREP/${trait}.normalized.tsv.gz" \
    --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 \
    --N "$n" --merge-alleles "$HM3" --chunksize 500000 \
    --out "$OUT/$trait"
}

munge_ukb() {
  local trait="$1"
  "$PY" "$MUNGE" \
    --sumstats "$RAW/${trait}.covararraypcs.assoc.logistic.gz" \
    --snp SNP --a1 A1 --a2 A0 --p P --signed-sumstats OR,1 \
    --N-col NMISS --frq A1FREQ --merge-alleles "$HM3" --chunksize 500000 \
    --out "$OUT/$trait"
}

"$PY" "$MUNGE" \
  --sumstats "$RAW/DEP_FinnGen_R12.gz" \
  --snp rsids --a1 alt --a2 ref --p pval --signed-sumstats beta,0 \
  --frq af_alt --N-cas 59333 --N-con 434831 \
  --merge-alleles "$HM3" --chunksize 500000 \
  --out "$OUT/DEP_FinnGen_R12"

munge_delange CD_deLange2017 40266
munge_delange UC_deLange2017 45975

munge_ukb UKB_LifetimeMDD
munge_ukb UKB_MDDRecur
munge_ukb UKB_GPpsy
munge_ukb UKB_ICD10Dep

find "$OUT" -maxdepth 1 -name '*.sumstats.gz' -printf '%f\t%s\n' | sort
