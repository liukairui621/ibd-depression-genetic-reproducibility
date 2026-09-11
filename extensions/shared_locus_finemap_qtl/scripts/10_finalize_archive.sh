#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl}"
cd "$ROOT"

mkdir -p provenance reports

{
  printf 'archive_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host\t%s\n' "$(hostname)"
  printf 'plink\t'
  /root/anaconda3/bin/plink --version 2>&1 | head -n 1
  Rscript -e 'cat("R\t", R.version.string, "\n", sep=""); for (p in c("data.table","susieR","coloc")) cat(p, "\t", as.character(packageVersion(p)), "\n", sep="")'
  software/qtl_env/bin/python - <<'PY'
import platform
from importlib.metadata import version
print(f"python\t{platform.python_version()}")
for package in ("pandas", "numpy", "matplotlib", "pysam", "requests", "certifi"):
    print(f"{package}\t{version(package)}")
PY
  printf 'ld_reference\t1000 Genomes Phase 3 EUR\n'
  printf 'gene_annotation\tGENCODE v19\n'
  printf 'random_seed\t20260817\n'
} > provenance/software_versions.txt

printf 'path\tsize_bytes\tmodified_utc\n' > reports/FILE_INVENTORY.tsv
find data figures logs provenance reports results scripts -type f \
  ! -path 'provenance/CHECKSUMS.sha256' \
  ! -path 'reports/CHECKSUM_VERIFICATION.log' \
  -printf '%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TSZ\n' \
  | LC_ALL=C sort >> reports/FILE_INVENTORY.tsv

find data figures logs provenance reports results scripts -type f \
  ! -path 'provenance/CHECKSUMS.sha256' \
  ! -path 'reports/CHECKSUM_VERIFICATION.log' \
  -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > provenance/CHECKSUMS.sha256

sha256sum -c provenance/CHECKSUMS.sha256 > reports/CHECKSUM_VERIFICATION.log
printf 'Checked files: %s\n' "$(wc -l < provenance/CHECKSUMS.sha256)" >> reports/CHECKSUM_VERIFICATION.log
