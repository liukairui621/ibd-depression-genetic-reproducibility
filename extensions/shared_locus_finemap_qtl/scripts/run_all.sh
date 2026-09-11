#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl}"
cd "$ROOT"

bash scripts/00_run_gwas_stage.sh
python3 scripts/04_map_candidate_genes.py --root "$ROOT"
software/qtl_env/bin/python scripts/05_retrieve_prepare_qtl.py --root "$ROOT"
Rscript scripts/06_run_qtl_coloc.R "$ROOT"
software/qtl_env/bin/python scripts/07_query_opentargets_pqtl.py --root "$ROOT"
bash scripts/07b_check_mst1_pqtl_ld.sh "$ROOT"
software/qtl_env/bin/python scripts/08_summarize_and_plot.py --root "$ROOT"
software/qtl_env/bin/python scripts/09_validate_results.py --root "$ROOT"
bash scripts/10_finalize_archive.sh
