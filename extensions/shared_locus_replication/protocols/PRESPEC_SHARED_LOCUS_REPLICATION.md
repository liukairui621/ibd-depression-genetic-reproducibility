# Prespecified cross-cohort shared-locus analysis

Frozen: 2026-08-17, before inspection of any PLACO/PLACO+ result from this extension.

## Question

Do variant-level pleiotropic associations between depression and inflammatory bowel disease (IBD) reproduce when both phenotype cohorts are replaced?

This is an exploratory extension of the reproducibility project, but the rules below are prospective for this analysis. A new list of shared loci is not itself treated as novel because shared-locus analyses have already been reported for IBD and psychiatric phenotypes. The intended contribution is cross-cohort reproducibility.

## Frozen inputs

All inputs are European-ancestry GWAS summary statistics already quality controlled with LDSC `munge_sumstats.py`, restricted to HapMap3-compatible, non-ambiguous common SNPs, and aligned to the same HapMap3 allele list.

1. Howard et al. public depression meta-analysis excluding 23andMe: 170,756 cases and 329,443 controls (N=500,199).
2. de Lange et al. 2017 IBD meta-analysis: 25,042 cases and 34,915 controls (N=59,957).
3. FinnGen R12 depression endpoint: 59,333 cases and 434,831 controls (N=494,164).
4. FinnGen R12 strict IBD endpoint: 10,960 cases and 489,388 controls (N=500,348).

The exact paths, real paths, SHA256 checksums, schemas, and retained SNP counts will be written to `manifests/INPUT_MANIFEST.tsv`.

## Analysis pairs

### Primary independent cohort split

- `External_Howard_deLange`: Howard depression x de Lange IBD. Original PLACO is used because the two contributing studies are treated as non-overlapping.
- `FinnGen_DEP_IBD`: FinnGen depression x FinnGen strict IBD. PLACO+ is used because both GWAS arise from FinnGen and may share participants and controls. The Z-score correlation is estimated with the official `cor.pearson()` rule: exclude a SNP if either marginal P is below 1e-4.

These two pair analyses use different cohort sources and form the primary cross-cohort replication comparison.

### Secondary diagonal sensitivity

- `Diagonal_Howard_FinnGenIBD`: Howard depression x FinnGen strict IBD.
- `Diagonal_FinnGenDEP_deLange`: FinnGen depression x de Lange IBD.

Original PLACO is used within each diagonal pair because the two studies in a pair are non-overlapping. These analyses are sensitivity checks only: because the two diagonal analyses reuse FinnGen across different phenotypes, agreement between them is not labelled independent replication.

## Harmonization and test implementation

1. Retain SNPs with numeric Z scores in both traits of a pair.
2. Harmonize the second trait to the first trait's displayed effect allele; reverse its Z score when alleles are swapped. Exclude allele mismatches.
3. Estimate trait-specific Z-score variances using the official PLACO `var.placo()` exclusion rule: remove variants for which both marginal P values are below 1e-4.
4. For PLACO+, estimate Z-score correlation using the official `cor.pearson()` exclusion rule described above.
5. Use the unmodified official PLACO v0.2.0 R source. Numerical integration is performed by `placo()` or `placo.plus()` as specified for each pair.
6. To avoid millions of unnecessary integrations, derive separate positive- and negative-product cutoffs corresponding to PLACO P=1e-4. Only variants beyond those cutoffs are integrated. Monotonicity is checked numerically, and exact PLACO P values are retained only after integration. This screen cannot omit a genome-wide significant result.
7. Use two-sided PLACO/PLACO+ P values. Genome-wide significance is P<5e-8.

## Locus definition

Genome-wide significant variants are assigned to the frozen 2,495 LAVA/1000G EUR LD blocks on GRCh37. Within each block and pair, the SNP with the smallest PLACO P is the index SNP. The extended MHC (chr6:25-34 Mb) is reported separately and cannot satisfy the primary replication endpoint.

## Primary replication endpoint

The two primary pair analyses are evaluated in both directions to avoid selecting the more favourable discovery cohort after seeing results. Family-wise alpha is split across directions.

A non-MHC locus is called cross-cohort replicated when:

1. its index SNP reaches P<5e-8 in the discovery pair;
2. the exact index SNP is present in the other primary pair;
3. its validation PLACO/PLACO+ P is below 0.025 divided by the number of non-MHC discovery blocks in that direction;
4. both validation marginal trait P values are below 0.05; and
5. the sign of `Z_depression x Z_IBD` is the same in discovery and validation.

Replicated loci from the two directions are deduplicated by LD block. No LD proxy substitution is allowed because all inputs use the same HapMap3 SNP backbone.

## Secondary summaries

- P<1e-6 loci are reported as exploratory only.
- Same-block P<1e-6 signals in both primary pairs are described as exploratory convergence, not replication.
- The two diagonal pairs are queried at primary lead SNPs and summarized as sensitivity evidence. They cannot rescue a failed primary replication endpoint.
- MHC findings are tabulated separately.

## Go/no-go rule

- `GO`: at least one non-MHC locus satisfies the frozen primary replication endpoint. Only those loci may enter later gene mapping, TWAS, colocalization, or functional follow-up.
- `NO-GO`: no non-MHC locus satisfies the endpoint. In that case no candidate gene is promoted from this analysis, even if a pair-specific or exploratory signal is present.

No manuscript conclusion, figure, or table is changed automatically by this run.
