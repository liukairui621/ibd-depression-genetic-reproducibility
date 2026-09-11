# Prespecified protocol: depression-definition by IBD-subtype extension

Frozen before inspecting any UKB-definition x IBD-subtype LDSC result.

## Objective

Determine whether variation in the observed cross-trait genetic correlation is
more closely associated with depression phenotype definition or with IBD
subtype. This is an attribution of variability in observed LDSC genetic
correlations, not a causal decomposition of biological mechanisms.

## Trait panel

Depression definitions:

1. Howard 2019 major-depression meta-analysis.
2. FinnGen R12 F5_DEPRESSIO.
3. UK Biobank lifetime MDD.
4. UK Biobank recurrent MDD.
5. UK Biobank GP consultation for nerves, anxiety, tension, or depression.
6. UK Biobank linked-record ICD-10 depression.

IBD phenotypes:

1. IBD, Crohn disease, and ulcerative colitis from de Lange 2017.
2. Strict IBD, strict Crohn disease, and strict ulcerative colitis from
   FinnGen R12.

UK Biobank depression definitions share participants and are phenotype
sensitivity analyses, not independent replication cohorts. FinnGen
depression x FinnGen IBD-family estimates are within-cohort sensitivity
analyses because participant and control overlap is expected.

## Primary comparison structure

The extension is evaluated on two orthogonal axes.

1. Depression-definition axis: compare the four UK Biobank depression
   definitions against each fixed IBD phenotype. This isolates phenotype
   definition within one depression cohort as closely as the public summary
   data permit.
2. IBD-subtype axis: compare IBD, Crohn disease, and ulcerative colitis within
   each fixed IBD cohort and depression definition.

Cross-cohort reciprocal anchors from the completed core analysis remain:

1. Howard 2019 MDD x FinnGen R12 IBD family.
2. FinnGen R12 depression x de Lange 2017 IBD family.

## Analysis

All traits are munged against the same HapMap3 allele list and analyzed using
the same European LD-score reference and weights. For every pair, report rg,
SE, z, P, BH-FDR across the full extension matrix, both trait heritabilities
and intercepts, genetic covariance, and genetic-covariance intercept.

The pairwise matrix is summarized without treating overlapping estimates as
independent:

1. Within each fixed IBD phenotype, calculate the range, standard deviation,
   and median absolute deviation of rg across the four UK Biobank depression
   definitions.
2. Within each fixed depression definition and IBD cohort, calculate the same
   quantities across IBD, Crohn disease, and ulcerative colitis.
3. Report the distributions of these within-stratum dispersion statistics and
   their ratio (depression-definition dispersion / IBD-subtype dispersion).
4. Use leave-one-definition and leave-one-subtype sensitivity summaries.
5. Pairwise Cochran-style contrasts are used only where sampling covariance is
   available from a joint multivariable LDSC model. Otherwise, comparisons are
   descriptive and no exact percentage attribution is claimed.

## Interpretation rules

1. Depression-definition variability is considered larger only if its
   dispersion is consistently larger across both de Lange and FinnGen IBD
   anchors and is not driven solely by the broad GP consultation phenotype.
2. IBD-subtype variability is considered larger only if CD-versus-UC
   separation is directionally consistent across multiple fixed depression
   definitions and both IBD cohorts.
3. If dispersion ratios vary materially by IBD cohort or leave-one-out
   analysis, the result is classified as mixed or indeterminate.
4. Statistical significance alone is not used to classify stability because
   power and sample size differ markedly between definitions and subtypes.
5. No local LAVA, SuSiE, colocalization, or molecular-QTL analysis is added in
   this extension.

## Reproducibility

The project stores scripts, commands, logs, source URLs, checksums, software
versions, intermediate files, raw LDSC outputs, parsed tables, figures, and
reports on the analysis server only.
