# Verification status

Revision date: 14 September 2026.

| Check | Result | Scope |
|---|---|---|
| Isolated workspace staging | Passed | Existing raw GWAS/reference assets and installed software; no cached statistical outputs copied |
| Six raw-to-munged GWAS reconstructions | Passed | All decompressed SHA256 hashes matched frozen corrected inputs |
| Core LDSC rerun | Passed | Six heritability rows and 15 genetic-correlation rows matched exported values |
| Portable downstream preflight | Passed | Stage commands, configured interpreters, extension output directories, archived-path substitution, and shell syntax on the original Linux server; analyses were not started |
| Local-analysis BH correction | Passed | 110 eligible tests; four q<0.05; three jointly estimable blocks; zero same-direction dual-significant replications |
| Phenotype-extension BH correction | Passed | 24 positive estimates, 15 q<0.05; unflagged subset 18 positive, nine q<0.05 |
| Annotation covariance BH correction | Passed | Ten frozen tests, two discovery q<0.05; model outputs contain 84 annotations |
| Full downstream end-to-end rerun | Not performed in this revision | LAVA, PLACO, molecular-QTL retrieval, multivariable extension, and S-LDSC retain archived execution results |
| Fresh software installation on an independent machine | Not performed | Environment versions and setup inputs are documented; no clean-room reproduction claim |
| GitHub publication | Public tag `v2.2.1-plos-one` | Submission-matched code and derived outputs |

The raw reconstruction and LDSC rerun were executed in a newly staged directory
on the original Linux server. The reference files and installed tools were
reused; source GWAS preprocessing was redone. This tests removal of private
preprocessed-input dependencies, not independence of the computing environment.
The downstream preflight validates only that the staged workflow can resolve its
entry points and required directories; it is not a downstream statistical rerun.

Run `python3 reproducibility/verify_published_results.py` to repeat the
export-level tests. `FROZEN_CORE_INPUT_HASHES.json` contains decompressed input
hashes for comparisons after raw reconstruction. Gzip container hashes may
differ because of compression timestamps even when content is identical.

The original full-branch execution protocols, versions, and results remain
available in their respective directories. Protocol-description errors and
interpretation corrections are enumerated in
`provenance/PROVENANCE_AMENDMENT_20260905.md`; old protocol wording should not
be read as a description of the corrected manuscript.
