#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
M="$ROOT/data/munged"
OUT="$ROOT/results/global_extension"
LOG="$ROOT/logs/20_run_ldsc_extension.log"
PY=/root/anaconda3/envs/ldsc/bin/python
LDSC=/root/ldsc/ldsc.py
LD=/root/IBD/00_RawData/Reference/LD_EUR/LDscore/LDscore.
W=/root/IBD/00_RawData/Reference/LD_EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
PAIR_FILE="$ROOT/provenance/extension_pair_manifest.tsv"
mkdir -p "$OUT/h2" "$OUT/rg" "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$LOG") 2>&1

DEP_TRAITS=(
  MDD_Howard2019
  DEP_FinnGen_R12
  UKB_LifetimeMDD
  UKB_MDDRecur
  UKB_GPpsy
  UKB_ICD10Dep
)
IBD_TRAITS=(
  IBD_deLange2017
  CD_deLange2017
  UC_deLange2017
  IBD_FinnGen_R12
  CD_FinnGen_R12
  UC_FinnGen_R12
)

declare -A DEP_COHORT=(
  [MDD_Howard2019]=Howard_2019
  [DEP_FinnGen_R12]=FinnGen_R12
  [UKB_LifetimeMDD]=UK_Biobank_Cai2020
  [UKB_MDDRecur]=UK_Biobank_Cai2020
  [UKB_GPpsy]=UK_Biobank_Cai2020
  [UKB_ICD10Dep]=UK_Biobank_Cai2020
)
declare -A IBD_COHORT=(
  [IBD_deLange2017]=IIBDGC_deLange2017
  [CD_deLange2017]=IIBDGC_deLange2017
  [UC_deLange2017]=IIBDGC_deLange2017
  [IBD_FinnGen_R12]=FinnGen_R12
  [CD_FinnGen_R12]=FinnGen_R12
  [UC_FinnGen_R12]=FinnGen_R12
)
declare -A IBD_SUBTYPE=(
  [IBD_deLange2017]=IBD
  [CD_deLange2017]=CD
  [UC_deLange2017]=UC
  [IBD_FinnGen_R12]=IBD
  [CD_FinnGen_R12]=CD
  [UC_FinnGen_R12]=UC
)

run_h2() {
  local trait="$1"
  "$PY" "$LDSC" \
    --h2 "$M/$trait.sumstats.gz" \
    --ref-ld-chr "$LD" \
    --w-ld-chr "$W" \
    --out "$OUT/h2/$trait"
}
export PY LDSC M LD W OUT
export -f run_h2
printf "%s\n" "${DEP_TRAITS[@]}" "${IBD_TRAITS[@]}" \
  | sort -u \
  | xargs -n1 -P4 bash -c 'run_h2 "$0"'

printf "pair_class\ttrait1\ttrait2\tdepression_cohort\tibd_cohort\tibd_subtype\toverlap_class\tprimary_extension\n" \
  > "$PAIR_FILE"

for dep in "${DEP_TRAITS[@]}"; do
  for ibd in "${IBD_TRAITS[@]}"; do
    dep_cohort="${DEP_COHORT[$dep]}"
    ibd_cohort="${IBD_COHORT[$ibd]}"
    overlap=no_known_overlap
    if [[ "$dep_cohort" == FinnGen_R12 && "$ibd_cohort" == FinnGen_R12 ]]; then
      overlap=expected_participant_and_control_overlap
    elif [[ "$dep_cohort" == UK_Biobank_Cai2020 && "$ibd_cohort" == IIBDGC_deLange2017 ]]; then
      overlap=shared_UK_sampling_cannot_be_excluded
    elif [[ "$dep_cohort" == Howard_2019 && "$ibd_cohort" == IIBDGC_deLange2017 ]]; then
      overlap=meta_components_not_individually_auditable
    fi
    primary=0
    if [[ "$dep_cohort" == UK_Biobank_Cai2020 ]]; then
      primary=1
    fi
    printf "cross_disease\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$dep" "$ibd" "$dep_cohort" "$ibd_cohort" \
      "${IBD_SUBTYPE[$ibd]}" "$overlap" "$primary" >> "$PAIR_FILE"
  done
done

UKB_TRAITS=(
  UKB_LifetimeMDD
  UKB_MDDRecur
  UKB_GPpsy
  UKB_ICD10Dep
)
for ((i=0; i<${#UKB_TRAITS[@]}; i++)); do
  for ((j=i+1; j<${#UKB_TRAITS[@]}; j++)); do
    printf "within_depression_ukb\t%s\t%s\tUK_Biobank_Cai2020\tNA\tNA\texpected_overlap\t0\n" \
      "${UKB_TRAITS[$i]}" "${UKB_TRAITS[$j]}" >> "$PAIR_FILE"
  done
done

for cohort in deLange2017 FinnGen_R12; do
  if [[ "$cohort" == deLange2017 ]]; then
    group=(IBD_deLange2017 CD_deLange2017 UC_deLange2017)
    cohort_name=IIBDGC_deLange2017
  else
    group=(IBD_FinnGen_R12 CD_FinnGen_R12 UC_FinnGen_R12)
    cohort_name=FinnGen_R12
  fi
  for ((i=0; i<${#group[@]}; i++)); do
    for ((j=i+1; j<${#group[@]}; j++)); do
      printf "within_ibd_cohort\t%s\t%s\tNA\t%s\tNA\texpected_overlap\t0\n" \
        "${group[$i]}" "${group[$j]}" "$cohort_name" >> "$PAIR_FILE"
    done
  done
done

printf "within_ibd_cross_cohort\tIBD_deLange2017\tIBD_FinnGen_R12\tNA\tmixed\tIBD\tno_known_overlap\t0\n" >> "$PAIR_FILE"
printf "within_ibd_cross_cohort\tCD_deLange2017\tCD_FinnGen_R12\tNA\tmixed\tCD\tno_known_overlap\t0\n" >> "$PAIR_FILE"
printf "within_ibd_cross_cohort\tUC_deLange2017\tUC_FinnGen_R12\tNA\tmixed\tUC\tno_known_overlap\t0\n" >> "$PAIR_FILE"

run_rg() {
  local pair_class="$1"
  local trait1="$2"
  local trait2="$3"
  "$PY" "$LDSC" \
    --rg "$M/$trait1.sumstats.gz,$M/$trait2.sumstats.gz" \
    --ref-ld-chr "$LD" \
    --w-ld-chr "$W" \
    --print-cov \
    --print-delete-vals \
    --out "$OUT/rg/${pair_class}__${trait1}__${trait2}"
}
export -f run_rg
tail -n +2 "$PAIR_FILE" \
  | cut -f1-3 \
  | tr "\t" " " \
  | xargs -n3 -P4 bash -c 'run_rg "$0" "$1" "$2"'

date --iso-8601=seconds > "$ROOT/provenance/extension_pairwise_ldsc_completed_at.txt"
