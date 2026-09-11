# Prespecified rules for replicated-locus fine-mapping and molecular colocalisation

Frozen: 2026-08-17, before inspecting any SuSiE or QTL-colocalisation output.

## Scope

Only three non-MHC LAVA blocks that passed reciprocal PLACO replication are eligible:

| locus | build-37 interval | directional PLACO leads |
|---|---|---|
| block154 | chr1:200134006-201067952 | rs169850; rs3861929 |
| block464 | chr3:47588462-50387742 | rs9862080 |
| block1671 | chr11:60515106-61717117 | rs174581 |

No locus may be added after examining the results. The analysis is a restricted follow-up of replicated variant-level pleiotropy, not a new genome-wide scan.

## GWAS fine-mapping

The four frozen inputs are Howard 2019 public depression excluding 23andMe (170756 cases, 329443 controls), FinnGen R12 F5_DEPRESSIO (59333/434831), de Lange 2017 IBD (25042/34915), and FinnGen R12 strict IBD (10960/489388). Variants are restricted to biallelic rsIDs present in the 1000 Genomes Phase 3 European reference, MAF >= 0.01, genotype missingness <= 0.05, non-palindromic alleles, and the four-GWAS common intersection in each frozen block. Alleles are aligned to the PLINK BIM allele order.

Primary fine-mapping uses SuSiE-RSS through coloc::runsusie with L=10, 95% credible sets, maxit=2000, external 1000G EUR LD, estimate_residual_variance=FALSE, and seed 20260817. A ridge-adjusted LD matrix, R*0.9999 + I*0.0001, is a numerical sensitivity analysis. Results are interpretable only when SuSiE converges and LD/summary-statistic mismatch diagnostics are reported.

## GWAS-GWAS colocalisation and replication

The six frozen pair types are Howard-de Lange, FinnGen depression-FinnGen IBD, the two crossed-cohort disease pairs, Howard-FinnGen depression, and de Lange-FinnGen IBD. Multi-signal coloc.susie is primary with p1=p2=1e-4 and p12=1e-5; p12=1e-6 is a prior sensitivity. Single-causal coloc.abf is secondary only.

PP.H4 >= 0.80 is strong shared-component support; 0.50-0.80 is suggestive; <0.50 is no support. Cross-cohort credible-set replication requires an exact shared variant or top variants in 1000G EUR LD r2 >= 0.80. A PLACO-replicated block with no such result is classified as pleiotropic but fine-mapping-unresolved, not as a shared causal mechanism.

## Candidate-gene mapping

Candidate genes are mapped only after GWAS fine-mapping. Mapping evidence is kept in separate columns: protein-coding gene overlap, distance within 100 kb of a 95% credible-set variant, Open Targets variant-to-gene/L2G evidence when available, and molecular-QTL colocalisation. Nearest-gene assignment alone cannot nominate a gene.

## Molecular layers and tissues

Eligible molecular layers are cis-eQTL, cis-sQTL, and plasma pQTL. The primary tissue set is whole blood, monocyte/immune-cell datasets, colon sigmoid, colon transverse, terminal ileum, and brain amygdala. Other brain tissues and plasma are sensitivity contexts. Local QTL-SuSiE colocalisation uses the same H4 thresholds and both p12 priors. Precomputed Open Targets colocalisation is reported separately from locally recomputed coloc and is not treated as an independent replication of the same source data.

Implementation detail frozen before QTL retrieval: a molecular trait must retain at least 200 allele-aligned SNPs from the locus reference and have at least one nominal cis-QTL association at P < 1e-5 to enter local SuSiE. For leafcutter sQTL data, at most one junction per candidate gene and dataset is carried forward, selected by the smallest minimum cis-QTL P value; the number of screened junctions is retained and this layer is labelled exploratory. Expression-QTL analyses use the gene-level molecular trait directly. Plasma pQTL evidence is taken from the Open Targets precomputed colocalisation index and kept separate from local eQTL/sQTL results.

To control computation without suppressing negative results, coloc.abf is run for every eligible local QTL-GWAS pair. Multi-signal SuSiE is then run only when the primary-prior ABF result has PP.H4 >= 0.50 or PP.H3 >= 0.50, indicating evidence for shared or distinct association signals that warrants multi-signal resolution. All ABF rows and all pairs not entering SuSiE are retained. This is a computational gate, not a gene-evidence threshold.

## Evidence tiers

- Tier 1: cross-cohort GWAS credible-set replication plus H4 >= 0.80 molecular-QTL colocalisation for the same mapped gene in an eligible context.
- Tier 2: either cross-cohort GWAS credible-set replication or H4 >= 0.80 molecular-QTL colocalisation, with independent positional or functional mapping support.
- Tier 3: positional/nearest-gene evidence only, suggestive H4, or a single-cohort signal.

All negative and non-estimable analyses remain in the exported tables. No candidate is promoted from Tier 3 to preserve a narrative.
