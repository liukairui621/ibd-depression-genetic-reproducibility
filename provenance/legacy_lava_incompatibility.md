# Legacy LAVA incompatibility record

Recorded during Layer-2 execution on 2026-07-26, before the new local results
were available.

The earlier analysis under `/root/IBD/05_LAVA` is not reused as a replication
arm for this project for two reasons:

1. Its input table encoded Howard 2019 depression as `cases=807553, controls=NA`.
   LAVA 0.1.5 therefore treated depression as continuous rather than as a
   binary case-control phenotype. The corrected analysis uses the public
   genome-wide release excluding 23andMe: 170,756 cases and 329,443 controls.
2. The earlier analysis selected local univariate signals using BH FDR. The
   frozen reproducibility protocol uses Bonferroni correction across 2,495
   loci and four primary traits (`P < 0.05/(2495*4)`), followed by BH FDR
   across all 110 eligible bivariate tests from both reciprocal pairs.

Consequently, the previously reported NOVA1-, PANK2-, and FHIT-region findings
are historical exploratory results only. They are not carried forward unless
they independently satisfy the new reciprocal Layer-2 criteria.
