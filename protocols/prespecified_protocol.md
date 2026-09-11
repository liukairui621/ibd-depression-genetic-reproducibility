# Prespecified protocol

Frozen before reading any new UK Biobank, FinnGen-depression, de Lange CD/UC,
local-correlation, or molecular-QTL result from this project.

Freeze date: 2026-07-26

## Primary question

Does the positive IBD-depression genetic correlation reproduce in reciprocal
cross-cohort comparisons, and is any reproducible genome-wide signal carried by
stable local and molecular components?

## Cohort roles

### Depression

- Howard 2019 major depression: high-power meta-analysis discovery phenotype.
- FinnGen R12 F5_DEPRESSIO: EHR-derived independent replication phenotype.
- UK Biobank Cai 2020 definitions: within-cohort phenotype-definition panel,
  not independent replication of Howard 2019.
- BIONIC DSM-5 lifetime MDD: clinically harmonized future anchor, pending a
  machine-downloadable summary-statistics accession.

### IBD

- de Lange 2017 European IBD/CD/UC: consortium cohort.
- FinnGen R12 strict IBD/CD/UC: EHR-biobank cohort.
- FinnGen + de Lange meta-analysis: excluded from replication claims because
  it contains both component cohorts.

## Primary reciprocal replication pairs

1. Howard 2019 MDD x FinnGen strict IBD.
2. FinnGen depression x de Lange 2017 IBD.

Subtype extensions replace IBD with CD or UC. The two primary IBD pairs use
different depression cohorts and different IBD cohorts. Within-FinnGen
depression-IBD comparisons are sensitivity analyses because participants and
controls may overlap.

## Layer 0: QC

- SNP heritability z-score must be at least 4 to enter local analysis.
- LDSC intercept, ratio, SNP count, and mean chi-square are reported without
  selectively constraining the intercept.
- A genetic-covariance intercept flag is recorded when
  `abs(gcov_intercept / gcov_intercept_se) >= 2`.
- Alleles are aligned to HapMap3 merge alleles and the European LD-score
  reference used in the original project.

## Layer 1: genome-wide reproducibility

Primary endpoints are the two reciprocal IBD pairs. The family-wise threshold
is 0.025 for these two tests. BH FDR is additionally reported across all
cross-disease tests.

Classification:

- **Fully reproduced**: both reciprocal primary estimates are positive and
  pass P < 0.025.
- **Partially reproduced**: both estimates are positive; one passes P < 0.025
  and the other passes nominal P < 0.05, with no significant Cochran-Q
  heterogeneity.
- **Inconclusive**: directions agree but the replication confidence interval
  is wide and includes both zero and the discovery estimate.
- **Definition/cohort dependent**: directions disagree, heterogeneity is
  significant, or the replication confidence interval excludes the discovery
  estimate while including zero.

Failure to reject zero is not treated as proof of no shared liability.

## Layer 2: local reproducibility

Only positive Layer-1 pairs whose two traits pass Layer-0 QC enter LAVA.

- Local univariate tests use Bonferroni correction across loci and traits.
- Bivariate local tests are performed only where both traits pass the local
  univariate threshold.
- BH FDR is applied across eligible bivariate tests.
- A local signal is replicated only when the same LD block, or genomic
  intervals with at least 50% reciprocal overlap, has the same direction and
  q < 0.05 in both reciprocal cross-cohort comparisons.

## Layer 3: fine-mapping and colocalization

Only Layer-2 replicated loci enter multi-signal analysis.

- SuSiE is fitted with the same European LD reference and deterministic seed.
- All retained credible sets must pass convergence and purity checks.
- Primary `coloc.susie` support: PP.H4 >= 0.80.
- Sensitivity support: PP.H4 >= 0.50 under the stricter prior
  `p12 = 1e-6`, with the primary prior `p12 = 1e-5`.
- A locus-level claim requires a shared component to reproduce across the two
  reciprocal disease-cohort comparisons.

## Layer 4: molecular reproducibility

Molecular datasets are queried only for Layer-3 components. Evidence is kept in
separate resolution tracks:

- cell-type eQTL;
- bulk-tissue eQTL;
- context-specific and bulk sQTL;
- pQTL.

Results are classified as:

- same-resolution replicated;
- cross-resolution concordant;
- context-specific;
- unsupported;
- technically inconclusive.

A null bulk eQTL is not treated as direct refutation of a cell-type eQTL.
DENND1B is retained only as a case study if it lies in a replicated component;
it is not a prespecified centerpiece.

## Fixed software settings

- LDSC 1.0.1 in the existing `/root/anaconda3/envs/ldsc` environment.
- LAVA 0.1.5 and R 4.3.3.
- Random seed: 42 for stochastic downstream methods.
- Genome build: GRCh37 for all GWAS/local analyses.
- LD reference: 1000 Genomes European Phase 3.
