#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov
PROV="$ROOT/provenance"
CONFIG="$ROOT/config"
SCRIPTS="$ROOT/scripts"
REF=/root/IBD/00_RawData/Reference

mkdir -p "$PROV" "$ROOT/results" "$ROOT/reports" "$ROOT/figures" "$ROOT/logs"

if find "$ROOT/results" -mindepth 1 -print -quit | grep -q .; then
  echo "Refusing to freeze: results directory is not empty" >&2
  exit 20
fi

HOWARD=/root/IBD/20_Reproducibility_Ladder/corrections/20260806_howard_no23andme/work/munged/MDD_Howard2019_no23andMe.sumstats.gz
FG_DEP=/root/IBD/20_Reproducibility_Ladder/data/munged_core/DEP_FinnGen_R12.sumstats.gz
DELANGE_IBD=/root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_deLange2017.sumstats.gz
FG_IBD=/root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_FinnGen_R12.sumstats.gz

printf '%s  %s\n' 41f9475df6d0f74f0b3eb69a9167f34aae9e9160e87a977b439c6f2b87ff4024 "$HOWARD" | sha256sum -c -
printf '%s  %s\n' 878aabd8f1cf51654a553a82c9664fafb3c66d9fc6776ebb26c20f3086b6277b "$FG_DEP" | sha256sum -c -
printf '%s  %s\n' 3edc165150571dfe1395a147ed0f7a4e392e2900608e8d29f510b6e5177b3dbc "$DELANGE_IBD" | sha256sum -c -
printf '%s  %s\n' 9059c0b7662bf63815e7007e0b7177b85398a22715707c2594859bf3fcdb592a "$FG_IBD" | sha256sum -c -

for chr in $(seq 1 22); do
  test -s "$REF/SLDSC/baselineLD.${chr}.annot.gz"
  test -s "$REF/SLDSC/baselineLD.${chr}.l2.ldscore.gz"
  test -s "$REF/SLDSC/baselineLD.${chr}.l2.M_5_50"
  test -s "$REF/LD_EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.${chr}.l2.ldscore.gz"
  test -s "$REF/LD_EUR/1000G_Phase3_frq/1000G.EUR.QC.${chr}.frq"
done

zcat "$REF/SLDSC/baselineLD.1.annot.gz" | sed -n '1p' > "$PROV/baselineLD_header.tsv"
while IFS=$'\t' read -r order annotation role; do
  if [ "$order" = "order" ] || [ -z "$annotation" ]; then
    continue
  fi
  grep -qw "$annotation" "$PROV/baselineLD_header.tsv"
done < "$CONFIG/primary_annotations.tsv"

cat > "$PROV/input_manifest.tsv" <<EOF
trait_id\trole\tphenotype_definition\tn_cases\tn_controls\tn_total\tpath\tsha256
IBD_deLange2017\tdiscovery_IBD\tIBD European meta-analysis\t25042\t34915\t59957\t$DELANGE_IBD\t3edc165150571dfe1395a147ed0f7a4e392e2900608e8d29f510b6e5177b3dbc
DEP_FinnGen_R12\tdiscovery_depression\tFinnGen R12 F5_DEPRESSIO\t59333\t434831\t494164\t$FG_DEP\t878aabd8f1cf51654a553a82c9664fafb3c66d9fc6776ebb26c20f3086b6277b
IBD_FinnGen_R12\treplication_IBD\tFinnGen R12 IBD_STRICT\t10960\t489388\t500348\t$FG_IBD\t9059c0b7662bf63815e7007e0b7177b85398a22715707c2594859bf3fcdb592a
MDD_Howard2019_no23andMe\treplication_depression\tPublic PGC-UK Biobank depression meta-analysis excluding 23andMe\t170756\t329443\t500199\t$HOWARD\t41f9475df6d0f74f0b3eb69a9167f34aae9e9160e87a977b439c6f2b87ff4024
EOF

find "$REF/SLDSC" -maxdepth 1 -type f \( -name 'baselineLD.*.annot.gz' -o -name 'baselineLD.*.l2.ldscore.gz' -o -name 'baselineLD.*.l2.M_5_50' \) -print0 | sort -z | xargs -0 sha256sum > "$PROV/baselineLD_v2.2_manifest.sha256"
find "$REF/LD_EUR/1000G_Phase3_weights_hm3_no_MHC" -maxdepth 1 -type f -name 'weights.hm3_noMHC.*.l2.ldscore.gz' -print0 | sort -z | xargs -0 sha256sum > "$PROV/weights_noMHC_manifest.sha256"
find "$REF/LD_EUR/1000G_Phase3_frq" -maxdepth 1 -type f -name '1000G.EUR.QC.*.frq' -print0 | sort -z | xargs -0 sha256sum > "$PROV/frequency_manifest.sha256"

Rscript -e '.libPaths(c("/root/IBD/20_Reproducibility_Ladder/software/R_lib",.libPaths())); d<-packageDescription("GenomicSEM"); cat("R=",as.character(getRversion()),"\nGenomicSEM=",d$Version,"\nGenomicSEM_RemoteSha=",d$RemoteSha,"\n",sep="")' > "$PROV/software_versions.txt"
printf 'n_blocks=200\nexclude_cont=TRUE\nscale=observed\nrandom_seed=not_applicable_deterministic_jackknife\nweights=HapMap3_noMHC\n' >> "$PROV/software_versions.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$PROV/freeze_timestamp_utc.txt"

sha256sum \
  "$ROOT/PREANALYSIS_FREEZE.md" \
  "$CONFIG/primary_annotations.tsv" \
  "$SCRIPTS/00_validate_and_freeze.sh" \
  "$SCRIPTS/01_run_sldsc_pair.R" \
  "$SCRIPTS/02_summarize_results.R" \
  "$SCRIPTS/run_pipeline.sh" \
  > "$PROV/frozen_analysis_manifest.sha256"

sha256sum /root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl/reports/EXECUTION_REPORT.md > "$PROV/local_adjudication_link.sha256"
echo "Freeze complete: $(cat "$PROV/freeze_timestamp_utc.txt")"
