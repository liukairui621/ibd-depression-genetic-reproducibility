# Howard 2019 public-file identity audit

## Determination

The project input `/root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz` is the public `PGC_UKB_depression_genome-wide.txt` release, not the full meta-analysis including 23andMe.

Evidence:

- The compressed input SHA256 is `c422cbcdd4296958273fb38e4a29e98a498022acd21f3e4b13d93232511f17e3`.
- The decompressed file contains 8,483,301 data rows and the fields `MarkerName`, `A1`, `A2`, `Freq`, `LogOR`, `StdErrLogOR`, and `P`.
- The University of Edinburgh DataShare ReadMe states that the public genome-wide file excludes 23andMe and contains 8,483,301 variants from 170,756 cases and 329,443 controls (N = 500,199).
- The same ReadMe states that the full 807,553-person meta-analysis is publicly represented only by a 10,000-variant file because 23andMe restricts public release. The project does not hold a 23andMe data-transfer agreement or the full restricted file.

Authoritative source: University of Edinburgh DataShare, dataset DOI `10.7488/ds/2458`, item `https://datashare.ed.ac.uk/handle/10283/3203`, ReadMe `https://datashare.ed.ac.uk/server/api/core/bitstreams/71f87d2b-5137-46f8-8da3-09adaebca549/content`.

## Correction

The previously assigned counts of 246,363 cases, 561,190 controls, and N = 807,553 were incorrect for this file. All corrected metadata and binary-trait analyses use 170,756 cases, 329,443 controls, and N = 500,199. The public summary statistics are described as a PGC-UK Biobank depression meta-analysis excluding 23andMe, not as a strictly ascertained MDD cohort.
