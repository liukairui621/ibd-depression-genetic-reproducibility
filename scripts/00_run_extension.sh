#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
LOG="$ROOT/logs/00_run_extension.log"
mkdir -p "$ROOT/logs"
exec > >(tee "$LOG") 2>&1

bash "$ROOT/scripts/17_write_extension_source_manifest.sh"
if [[ ! -s "$ROOT/provenance/extension_download_completed_at.txt" ]]; then
  bash "$ROOT/scripts/18_resume_extension_downloads.sh"
fi
bash "$ROOT/scripts/19_prepare_extension.sh"
bash "$ROOT/scripts/20_run_ldsc_extension.sh"
/root/anaconda3/bin/python "$ROOT/scripts/21_parse_ldsc_extension.py"
Rscript "$ROOT/scripts/22_run_genomicsem_extension.R"
Rscript "$ROOT/scripts/09_primary_pair_dependence_diagnostic.R"
Rscript "$ROOT/scripts/23_attribute_instability.R"
Rscript "$ROOT/scripts/24_sensitivity_exclude_gcov_flags.R"
Rscript "$ROOT/scripts/24_plot_extension.R"
/root/anaconda3/bin/python "$ROOT/scripts/25_write_extension_report.py"
/root/anaconda3/bin/python "$ROOT/scripts/26_finalize_extension.py"
