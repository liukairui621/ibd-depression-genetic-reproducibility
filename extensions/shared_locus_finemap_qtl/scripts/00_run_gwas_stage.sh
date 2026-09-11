#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl
mkdir -p "${ROOT}"/{data/{prepared,ld},results/gwas,logs,provenance,scripts,reports}

python3 "${ROOT}/scripts/01_prepare_gwas_inputs.py" --root "${ROOT}" 2>&1 | tee "${ROOT}/logs/01_prepare_gwas_inputs.log"
bash "${ROOT}/scripts/02_build_ld.sh" "${ROOT}" 2>&1 | tee "${ROOT}/logs/02_build_ld.log"
Rscript "${ROOT}/scripts/03_run_gwas_susie_coloc.R" "${ROOT}" 2>&1 | tee "${ROOT}/logs/03_run_gwas_susie_coloc.log"
