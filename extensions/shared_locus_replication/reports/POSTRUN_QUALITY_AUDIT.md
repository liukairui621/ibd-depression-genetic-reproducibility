# Post-run quality audit

Generated: 2026-08-17T17:08:09.900697+00:00

Overall status: **PASS**.

## Full-data parameter re-estimation

The official PLACO `var.placo()` and `cor.pearson()` functions were rerun on all harmonized SNPs for both primary pairs. Both comparisons passed the absolute tolerance of 1e-10. Exact values are in `OFFICIAL_PARAMETER_VALIDATION.tsv`.

## External-pair PLACO+ sensitivity

The external Howard-de Lange pair was originally assigned regular PLACO under the frozen no-known-overlap assumption. Because its null-Z correlation was 0.0338, all integration candidates were re-scored with PLACO+ as a post-run robustness analysis. The external pair had 27 genome-wide-significant blocks under PLACO+. The frozen cross-cohort endpoint retained 3 replicated non-MHC blocks after replacing every external-pair discovery or validation P value with PLACO+.

## Trait-specific direction audit

The effect alleles were re-aligned across the two primary pairs for every frozen replicated index SNP. 4/4 directional replication rows had concordant depression effects and concordant IBD effects separately. This is stricter than comparing only the sign of the product statistic.

## Implementation correction

The first run's positional allele parser was corrected to use semantic LDSC `A1`/`A2` column names. The complete pipeline was rerun, and the pre-correction snapshot remains archived. See `provenance/IMPLEMENTATION_CORRECTION_A1_SEMANTICS.md`.

No gene, mechanism, or clinical claim is made by this audit.
