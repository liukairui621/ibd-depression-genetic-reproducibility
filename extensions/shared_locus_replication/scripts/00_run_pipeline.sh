#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication
mkdir -p "$ROOT/logs" "$ROOT/provenance"
exec > >(tee "$ROOT/logs/00_run_pipeline.log") 2>&1

date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ" | tee "$ROOT/provenance/started_at.txt"

python3 "$ROOT/scripts/01_harmonize_pairs.py"
Rscript "$ROOT/scripts/02_compute_placo_cutoffs.R"
python3 "$ROOT/scripts/03_filter_candidates.py"
Rscript "$ROOT/scripts/04_run_placo.R"
python3 "$ROOT/scripts/05_prepare_replication.py"
Rscript "$ROOT/scripts/06_score_replication_tests.R"
python3 "$ROOT/scripts/07_finalize.py"
Rscript "$ROOT/scripts/08_validate_official_parameters.R"
Rscript "$ROOT/scripts/09_external_placo_plus_sensitivity.R"
python3 "$ROOT/scripts/10_direction_audit.py"
python3 "$ROOT/scripts/11_write_postrun_audit.py"

{
  python3 --version
  Rscript --version
  bgzip --version | head -n 1
  tabix --version | head -n 1
  Rscript -e 'cat("data.table=", as.character(packageVersion("data.table")), "\n", sep="")'
  cd "$ROOT/software/PLACO"
  git rev-parse HEAD
  sha256sum PLACO_v0.2.0.R
} > "$ROOT/provenance/SOFTWARE_VERSIONS.txt" 2>&1

date -u +"completed_utc=%Y-%m-%dT%H:%M:%SZ" | tee "$ROOT/provenance/completed_at.txt"

find "$ROOT" -type f ! -path '*/.git/*' ! -name 'SHA256SUMS.tsv' ! -name '12_checksum_verification.log' -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$ROOT/provenance/SHA256SUMS.tsv"
sha256sum -c "$ROOT/provenance/SHA256SUMS.tsv" > "$ROOT/logs/12_checksum_verification.log"
echo "PIPELINE_COMPLETE"
