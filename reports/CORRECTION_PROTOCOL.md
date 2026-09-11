# Howard 2019 no-23andMe correction protocol

## Status and scope

This protocol was frozen before inspecting any corrected LDSC or LAVA results. It applies only to the Howard et al. depression input and downstream analyses that depend on it. Unaffected FinnGen depression x de Lange IBD raw results will not be recomputed, although pooled multiple-testing correction will be recalculated if the number or ordering of eligible local tests changes.

## Input correction

- Source summary statistics: `/root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz`.
- Public release identity: no-23andMe genome-wide file with 8,483,301 variants.
- Correct sample counts: 170,756 cases, 329,443 controls, total N = 500,199.
- Superseded counts: 246,363 cases, 561,190 controls, total N = 807,553. These counts describe the full meta-analysis including 23andMe and must not be assigned to the public no-23andMe file.
- The original raw file is immutable. A corrected munged file will be generated in the correction workspace, with logs and checksums retained.

## Prespecified analyses

1. Re-run single-trait LDSC heritability QC for Howard depression using the corrected sample size.
2. Re-run all global LDSC comparisons involving Howard depression. Report corrected estimates and compare them with the frozen pre-correction values.
3. Re-run Howard depression x FinnGen strict IBD LAVA using the corrected case-control counts.
4. Combine the corrected Howard x FinnGen eligible local tests with the unchanged FinnGen depression x de Lange IBD eligible tests and reapply BH-FDR across the pooled bivariate family.
5. Re-evaluate local estimability and cross-pair reproducibility using the original criteria: a block must be estimable in both primary pairings, have concordant direction, and meet pooled q < 0.05 in both pairings. Distinguish non-estimable, estimable but non-significant, and directionally discordant blocks.
6. Remove the two-estimate fixed-effect meta-analysis, Cochran Q, and I-squared from the manuscript and figures irrespective of corrected results. The primary correlations share FinnGen participants and are not statistically independent.
7. Repeat the global phenotype-extension summary after excluding the six cells with prespecified genetic-covariance-intercept flags. Treat this as a sensitivity analysis, not a replacement primary analysis.

## Interpretation boundaries

- Do not claim that smaller N must reduce the number of LAVA-eligible blocks; eligibility and local estimates will be described from the rerun.
- The global conclusion requires positive direction in both primary pairings, not a pooled effect estimate.
- A failure to identify a reproduced local component means only that none was detected under the available power and prespecified estimability/FDR criteria.
- No causal, mechanistic, individual-prediction, screening, treatment-selection, or candidate-gene claim will be inferred from this correction.
- No new molecular follow-up analysis will be added during this correction cycle.

## Required downstream updates

After numerical verification, update the manuscript, tables, figures, supplementary files, numeric trace, public-code snapshot, provenance amendment, checksums, and submission archive. Preserve the frozen pre-correction archive and clearly label all superseded artifacts.
