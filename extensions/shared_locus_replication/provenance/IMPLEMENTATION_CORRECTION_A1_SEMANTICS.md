# Implementation correction: LDSC A1 semantics

Detected: 2026-08-17 during the post-run allele-direction audit, before scientific interpretation or candidate promotion.

## Issue

The first pipeline run selected the first two non-statistic columns positionally and treated the first as the effect allele. The Howard and de Lange LDSC files are ordered `SNP A1 A2 Z N`, but the FinnGen files are ordered `SNP A2 A1 Z N`. LDSC Z is defined with respect to the column named `A1`, regardless of physical column order.

## Consequence

- The effect/other allele labels were wrong for FinnGen rows.
- Within the two primary pair analyses, PLACO/PLACO+ P values were unaffected because both traits in each pair used the same allele convention and the product statistic is unchanged when both Z scores are expressed for the opposite allele.
- In the diagonal sensitivity pairs, product signs and estimated null correlations could be wrong because the two files have different physical allele-column orders. PLACO P values based on the absolute product were not expected to change materially, but all downstream direction summaries required correction.

## Correction

`01_harmonize_pairs.py` now identifies effect and other alleles by the semantic column names `A1` and `A2`. The complete pipeline is rerun from harmonization, not patched at the result-table level.

The first-run project snapshot is preserved at:

`/root/IBD/20_Reproducibility_Ladder/extensions/archive/20260817_shared_locus_pre_A1_fix`

This document records an implementation correction, not a change to the frozen scientific thresholds or replication endpoint.
