#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov
PROV="$ROOT/provenance"

test -s "$ROOT/results/PRIMARY_ANNOTATION_DECISION.tsv"
test -s "$ROOT/reports/ANNOTATION_STRATIFIED_GCOV_REPORT.md"
test -s "$ROOT/figures/annotation_covariance_forest.pdf"
test -s "$ROOT/IMPLEMENTATION_AMENDMENT_001.md"

grep -v 'scripts/01_run_sldsc_pair.R' "$PROV/frozen_analysis_manifest.sha256" > "$PROV/frozen_unchanged_manifest.sha256"
sha256sum -c "$PROV/frozen_unchanged_manifest.sha256"
sha256sum -c "$PROV/amendment_001_manifest.sha256"

find "$ROOT" -type f \
  ! -path "$ROOT/software/*" \
  ! -path "$ROOT/reference/downloads/*" \
  ! -path "$PROV/output_manifest.sha256" \
  ! -path "$PROV/file_inventory.tsv" \
  ! -path "$PROV/completion_timestamp_utc.txt" \
  -print0 | sort -z | xargs -0 sha256sum > "$PROV/output_manifest.sha256"
find "$ROOT" -type f -printf '%P\t%s\n' | sort > "$PROV/file_inventory.tsv"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$PROV/completion_timestamp_utc.txt"

echo "Provenance finalised: $(cat "$PROV/completion_timestamp_utc.txt")"
