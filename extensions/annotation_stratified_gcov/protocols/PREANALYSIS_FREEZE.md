# Pre-analysis freeze: annotation-stratified genetic covariance

Freeze status: rules must be hashed before any `s_ldsc` result is generated.

## Scientific question

Does the small positive genome-wide genetic covariance between inflammatory bowel disease (IBD) and depression reproducibly localize to broad functional genomic annotations when both the IBD cohort and the depression cohort are changed?

This is a bounded global-to-functional follow-up. It does not reopen locus discovery. The existing three-locus fine-mapping project remains the local adjudication result and will not be altered or reinterpreted to nominate a shared causal gene.

## Dataset roles

Discovery pair:

- IBD: de Lange et al. 2017 IBD meta-analysis, 25,042 cases and 34,915 controls.
- Depression: FinnGen R12 F5_DEPRESSIO, 59,333 cases and 434,831 controls.

Reciprocal replication pair:

- IBD: FinnGen R12 IBD_STRICT, 10,960 cases and 489,388 controls.
- Depression: Howard et al. 2019 public PGC-UK Biobank depression meta-analysis excluding 23andMe, 170,756 cases and 329,443 controls.

The two pairs use four distinct GWAS datasets. The within-FinnGen depression-IBD pair is excluded because it does not provide cohort-independent replication and has a flagged cross-trait intercept in the preceding global analysis.

## Method

- Software: `GenomicSEM::s_ldsc`, installed version and git SHA recorded at execution.
- Reference: 1000 Genomes Phase 3 European baseline-LD v2.2.
- Regression weights: HapMap3 weights with the extended MHC excluded.
- Scale: observed scale (`sample.prev` and `population.prev` are `NA`).
- Jackknife: 200 deterministic blocks; no random seed is used by this analysis.
- Estimand: zero-order genetic covariance within each annotation. Genetic correlation is descriptive and reported only when both annotation-specific heritabilities are positive.
- All 97 baseline-LD v2.2 annotations are included jointly in the regression. Formal inference is restricted to the ten annotations frozen below.

Because the regression SNP set is intersected with `weights.hm3_noMHC`, the primary analysis is MHC-excluded by design. A second post-hoc HLA removal analysis will not be performed.

## Frozen annotation family

1. `Coding_UCSC`
2. `Conserved_LindbladToh`
3. `DHS_Trynka`
4. `Enhancer_Andersson`
5. `H3K27ac_Hnisz`
6. `H3K4me1_Trynka`
7. `Promoter_UCSC`
8. `SuperEnhancer_Hnisz`
9. `TFBS_ENCODE`
10. `Repressed_Hoffman`

The family captures coding, conserved, accessible, enhancer, promoter, transcription-factor-binding and repressed chromatin. `Repressed_Hoffman` is retained as a negative/comparator annotation. Flanking and continuous annotations are not members of the primary family.

## Decision rules

For each frozen annotation:

1. Discovery uses a two-sided covariance P value, with Benjamini-Hochberg correction across the ten annotations.
2. Replication must have the same covariance direction and a directional one-sided P value below 0.05.
3. Both trait-specific annotation heritabilities must be positive in both pairs; otherwise the annotation is labelled not estimable for replicated interpretation.
4. A **replicated annotation-level covariance** requires all three conditions above and discovery FDR below 0.05.
5. A **replicated concentration** additionally requires covariance enrichment above 1 in both pairs. Enrichment is the annotation share of genome-wide covariance divided by the annotation SNP proportion. It is secondary because the genome-wide covariance is small and the ratio can be unstable.
6. Discovery-only, replication-only, opposite-direction, or non-estimable results are not replicated findings.

## Stopping rule and claim boundary

- If no annotation meets the replicated annotation-level covariance rule, the result is reported as no reproducible functional concentration detected under the frozen model. No additional annotation family, cell type, pathway, locus or gene will be searched in this extension.
- If an annotation meets the rule, it supports localization of covariance to a broad genomic context, not a shared causal pathway, cell type, gene, clinical biomarker or treatment target.
- The prior fine-mapping conclusion remains: broad local-correlation blocks decomposed into cohort- or disease-specific molecular components and did not identify a cross-disease shared causal gene.

## Non-primary outputs

The full 97-annotation output, tau estimates and cross-trait intercepts are retained for audit. They are not eligible to change the prespecified decision.

