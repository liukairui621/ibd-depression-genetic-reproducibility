# Execution report

## Completed stages

1. Four-GWAS harmonisation and 1000 Genomes EUR LD construction.
2. SuSiE-RSS fine-mapping with ridge-free primary and 1e-4 ridge sensitivity analyses.
3. Cross-cohort credible-set overlap and GWAS-GWAS multi-signal colocalisation.
4. GENCODE v19 positional candidate mapping.
5. Local eQTL/sQTL retrieval, eligibility filtering, ABF screening and SuSiE colocalisation.
6. Open Targets pQTL query and LD separation of MST1 protein signals.
7. Evidence-tier tables and publication-ready summary figure.

## Counts

- Positional candidate rows: 78
- Eligible QTL traits: 78
- ABF result rows: 624
- Primary SuSiE gene-trait H4>=0.8 rows: 6
- Open Targets region-matched pQTL-GWAS rows: 614
- CEDAR technical retrieval failures: 9

## Software

- Python: 3.12.3
- pandas: 2.2.3
- numpy: 2.1.3
- matplotlib: 3.9.2
- R/SuSiE/coloc versions are recorded separately in `provenance/software_versions.txt`.
