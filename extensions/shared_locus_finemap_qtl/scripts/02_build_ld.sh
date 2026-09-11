#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?project root required}
PLINK=/root/anaconda3/bin/plink
REF=/root/IBD/00_RawData/Reference/LD_EUR/1000G_EUR_Phase3_plink
OUT=${ROOT}/data/ld
mkdir -p "${OUT}" "${ROOT}/results/gwas"

while read -r locus chrom; do
  "${PLINK}" --bfile "${REF}/1000G.EUR.QC.${chrom}" \
    --extract "${ROOT}/data/prepared/${locus}_common_snps.txt" \
    --maf 0.01 --geno 0.05 --keep-allele-order --allow-no-sex \
    --make-bed --out "${OUT}/${locus}_1000G_EUR"
  "${PLINK}" --bfile "${OUT}/${locus}_1000G_EUR" \
    --keep-allele-order --allow-no-sex --freq \
    --out "${OUT}/${locus}_1000G_EUR"
  "${PLINK}" --bfile "${OUT}/${locus}_1000G_EUR" \
    --keep-allele-order --allow-no-sex --r square bin \
    --out "${OUT}/${locus}_LD"
done <<'EOF'
block154 1
block464 3
block1671 11
EOF

printf 'locus\tn_variants\tn_samples\tld_bytes\n' > "${ROOT}/results/gwas/ld_build_qc.tsv"
for locus in block154 block464 block1671; do
  printf '%s\t%s\t%s\t%s\n' "${locus}" \
    "$(wc -l < "${OUT}/${locus}_1000G_EUR.bim")" \
    "$(wc -l < "${OUT}/${locus}_1000G_EUR.fam")" \
    "$(stat -c '%s' "${OUT}/${locus}_LD.ld.bin")" \
    >> "${ROOT}/results/gwas/ld_build_qc.tsv"
done
