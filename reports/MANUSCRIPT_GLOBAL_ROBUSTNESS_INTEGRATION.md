# Manuscript integration: global robustness analysis

## Editorial decision

The UK Biobank depression-definition and CD/UC extension should be presented
as a global robustness analysis. Its role is to test whether the direction and
magnitude of the genome-wide IBD-depression correlation depend on one
depression definition or one IBD phenotype. It does not support a main
conclusion that depression phenotype definition or IBD subtype drives the
observed magnitude variability.

## Abstract replacement sentence

Across four UK Biobank depression definitions and IBD, Crohn disease, and
ulcerative colitis phenotypes from two IBD cohorts, all 24 genetic-correlation
estimates were positive, although variation in magnitude could not be
attributed primarily to depression definition or IBD subtype.

## Methods addendum

We evaluated global robustness by crossing four UK Biobank depression
definitions (lifetime major depressive disorder, recurrent major depressive
disorder, general-practitioner consultation for nerves, anxiety, tension, or
depression, and linked-record ICD-10 depression) with IBD, Crohn disease, and
ulcerative colitis GWAS from de Lange et al. and FinnGen R12. Summary statistics
were harmonized to the same HapMap3 allele list and analyzed using European
LD-score references. Pairwise LDSC estimates included genetic correlation,
standard error, P value, FDR, heritability, trait intercepts, genetic
covariance, and genetic-covariance intercept. Because the UK Biobank
definitions share participants, inferential comparisons used the sampling
covariance matrix from a 12-trait multivariable LDSC model implemented in
GenomicSEM. We compared adjusted dispersion across depression definitions and
IBD subtypes using 50,000 multivariate-normal draws (seed 42), together with
leave-one-definition, leave-one-subtype, and cross-cohort CD-versus-UC
direction checks.

## Results paragraph

The positive genome-wide association was robust to alternative depression
definitions and IBD phenotypes. All 24 UK Biobank depression-by-IBD estimates
were positive (`rg=0.045-0.230`), of which 15 passed BH-FDR < 0.05. Adjusted
marginal correlations ranged from 0.111 for the broad GP-consultation
definition to 0.163 for recurrent MDD and from 0.117 for ulcerative colitis to
0.153 for Crohn disease. The adjusted standard deviation across depression
definitions was 0.025, compared with 0.020 across IBD subtypes, giving a
definition-to-subtype dispersion ratio of 1.26 (95% CI 0.38-6.21). This
tendency was not robust to exclusion of the GP-consultation definition
(ratio=0.985). Although the IBD-subtype omnibus test was significant
(`P=0.005`), CD had a higher correlation than UC for all four definitions in
de Lange but a lower correlation for three of four definitions in FinnGen.
Neither pairwise differences among depression definitions nor CD-versus-UC
contrasts survived FDR correction. The prespecified classification was
therefore mixed or indeterminate.

## Discussion replacement paragraph

This extension distinguishes robustness of direction from stability of
magnitude. The uniformly positive estimates indicate that the genome-wide
IBD-depression correlation is not an artifact of a single depression
definition or a single IBD phenotype. However, the wide interval around the
dispersion ratio, sensitivity to removal of the broad GP-consultation
definition, and opposite CD-versus-UC ordering across IBD cohorts prevent a
credible attribution of magnitude variability to either axis. The results
therefore strengthen the conclusion of a small positive genome-wide overlap
while narrowing its interpretation: they do not support depression-definition
specificity, CD- or UC-specific shared architecture, or a claim that either
phenotype axis explains the absence of replicated local genetic components.

## Conclusion replacement sentence

IBD and depression show a small positive genome-wide genetic overlap whose
direction is robust across the tested phenotype definitions and IBD subtypes,
but whose magnitude variability cannot be assigned primarily to either axis.

## Claims removed from the main narrative

The main text should not describe the overlap as
"cohort/definition-dependent," state that broad depression phenotyping drives
the observed instability, or infer stronger shared architecture for CD or UC.
The absence of replicated local components should instead be interpreted as
compatible with diffuse polygenicity or limited local power, without assigning
that result to phenotype definition or IBD subtype.
