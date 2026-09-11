#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder
OUT="$ROOT/provenance/extension_source_manifest.tsv"
mkdir -p "$ROOT/provenance"

printf "trait\tcategory\tcohort\tphenotype\tcases\tcontrols\ttotal_n\tfile_id\texpected_bytes\texpected_md5\tsource_url\tpmid\tdoi\tnote\n" > "$OUT"
printf "UKB_LifetimeMDD\tdepression\tUK_Biobank_Cai2020\tLifetime MDD using DSM symptom and impairment criteria\t16301\t50870\t67171\t21356610\t763812382\tccb58b7a3d3b2a32d0ac651fa4282bbb\thttps://figshare.com/ndownloader/files/21356610\t32231276\t10.1038/s41588-020-0594-5\tPer-variant NMISS is used by LDSC\n" >> "$OUT"
printf "UKB_MDDRecur\tdepression\tUK_Biobank_Cai2020\tRecurrent MDD\t10302\t49083\t59385\t21356598\t763547170\t2f6730c7e95583691b8e796ec0df2830\thttps://figshare.com/ndownloader/files/21356598\t32231276\t10.1038/s41588-020-0594-5\tPer-variant NMISS is used by LDSC\n" >> "$OUT"
printf "UKB_GPpsy\tdepression\tUK_Biobank_Cai2020\tGP consultation for nerves anxiety tension or depression\t113262\t219360\t332622\t21356688\t777948024\td1677647fa42f4f59f27b2aa296a141e\thttps://figshare.com/ndownloader/files/21356688\t32231276\t10.1038/s41588-020-0594-5\tBroad help-seeking phenotype; per-variant NMISS is used\n" >> "$OUT"
printf "UKB_ICD10Dep\tdepression\tUK_Biobank_Cai2020\tLinked-record ICD-10 depression\t9176\t203235\t212411\t21356586\t772599146\t1337d06b801074403d2fc1b74ed53a7f\thttps://figshare.com/ndownloader/files/21356586\t32231276\t10.1038/s41588-020-0594-5\tPer-variant NMISS is used by LDSC\n" >> "$OUT"
printf "CD_deLange2017\tIBD\tIIBDGC_deLange2017\tCrohn disease\t12194\t28072\t40266\tGCST004132\t296840782\tNA\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/harmonised/28067908-GCST004132-EFO_0000384.h.tsv.gz\t28067908\t10.1038/ng.3760\tFixed total N is used because the harmonized file lacks per-variant N\n" >> "$OUT"
printf "UC_deLange2017\tIBD\tIIBDGC_deLange2017\tUlcerative colitis\t12366\t33609\t45975\tGCST004133\t299230653\tNA\thttps://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz\t28067908\t10.1038/ng.3760\tFixed total N is used because the harmonized file lacks per-variant N\n" >> "$OUT"

sha256sum "$OUT" > "$ROOT/provenance/extension_source_manifest.sha256"
