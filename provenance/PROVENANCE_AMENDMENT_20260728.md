# Provenance Amendment — 2026-07-28

## Purpose

This amendment resolves the record-hygiene issues identified during the read-only Claude Code audit. It does not alter GWAS inputs, analysis scripts, statistical outputs, figures, or scientific conclusions.

## Frozen core manifest

`provenance/final_core_file_sha256.tsv` remains unchanged as the historical checksum record of the frozen core analysis. Two reports were subsequently extended when the phenotype-definition analysis was integrated:

| File | Core-freeze SHA256 | Core-freeze bytes | Current SHA256 | Current bytes |
|---|---:|---:|---:|---:|
| `reports/EXECUTION_REPORT.md` | `19079a881ae50131a485f61b60233504325519ff20bc9c92a662ca7ca0f637cf` | 3829 | `c1e5a5b7237e98a9c6a1383c84c80ba0a3b76e2262275f838c60e8d88689816e` | 4990 |
| `reports/layer2_evidence_report.md` | `b81758bf333d34a5548cd49670865cea47bf7c97e3497a78184bc607dd0ad1b5` | 3188 | `ea466a757d686b5091a89909d1de2a22a7a6f9c1a7a930cb8f82e9ba730c262d` | 3397 |

The pre-integration versions remain preserved under `reports/archive/pre_global_robustness_integration_20260728/`. The current report versions and both archived versions are recorded in `provenance/extension_final_sha256.tsv`. The core manifest is therefore retained as a historical freeze rather than regenerated after the fact.

## Project status correction

`results/project_status.tsv` was corrected after the extension had completed:

| Item | Previous record | Current record |
|---|---|---|
| phenotype extension | `pending_source_transfer` | `completed_global_robustness` |
| reason | UKB definition panel and de Lange CD/UC downloads incomplete | 24/24 estimates were positive; 15/24 passed BH-FDR; instability attribution was mixed_or_indeterminate |
| SHA256 | `865aa0cb9d2c3cb8d91e21e14d7add76daf8f5265671a37a8dca7aa5c0bd3491` | `35d03c72c5c7597732fe5fea1fa617dc954447360300b16e8520475915492e9b` |

This is a metadata correction only. It records results already frozen in `reports/EXTENSION_UKB_DEPRESSION_CD_UC_REPORT.md`, `reports/MANUSCRIPT_GLOBAL_ROBUSTNESS_INTEGRATION.md`, and `reports/EXTENSION_FINAL_AUDIT.md`.

## Correlation-estimate provenance

The manuscript and manuscript-review handoff distinguish the two correlation sources:

1. Pairwise LDSC estimates describe the 24-cell genetic-correlation matrix, nominal and FDR significance, ranges, and cross-trait intercept quality-control flags.
2. Joint multivariable LDSC implemented through GenomicSEM supplies the sampling covariance and standardized correlation estimates used for marginal means, dispersion estimates, omnibus tests, and contrasts.

These are complementary outputs and are not interchangeable.

## Additional reporting decisions

- The identity sample-overlap matrix used in LAVA is reported as a limitation because exact cohort-overlap fractions were unavailable.
- The UK Biobank ICD-10 depression definition remains in the phenotype-sensitivity panel because it passed the prespecified heritability gate, but its lower mean chi-square and wider uncertainty are acknowledged.
- Existing CRLF checksum files are retained byte-for-byte. New audit-package manifests are written with LF line endings.
- Source-data symlinks are retained for computational continuity. The manuscript audit relies on source manifests, checksums, and retrieval records rather than assuming that symlinks are portable.

