#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov

test -s "$ROOT/results/PRIMARY_ANNOTATION_DECISION.tsv"
test -s "$ROOT/reports/ANNOTATION_STRATIFIED_GCOV_REPORT.md"
test -s "$ROOT/figures/annotation_covariance_forest.pdf"

find "$ROOT/results" "$ROOT/reports" "$ROOT/figures" "$ROOT/logs" -type f -print0 | sort -z | xargs -0 sha256sum > "$ROOT/provenance/output_manifest.sha256"
find "$ROOT" -type f -printf '%P\t%s\n' | sort > "$ROOT/provenance/file_inventory.tsv"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$ROOT/provenance/completion_timestamp_utc.txt"

sha256sum -c "$ROOT/provenance/frozen_analysis_manifest.sha256"
echo "Provenance finalised: $(cat "$ROOT/provenance/completion_timestamp_utc.txt")"
