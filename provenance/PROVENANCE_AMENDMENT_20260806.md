# Provenance amendment: Howard public release correction

Date: 2026-08-06

The pre-correction project and submission package were frozen before any file
was replaced. Historical checksum files are retained with the
`historical_*_pre_correction` prefix and must not be used to validate the
current tree.

The Howard source file was identified as the public genome-wide release
excluding 23andMe by its repository documentation, schema, row count, and
checksum. It was remunged with 170,756 cases and 329,443 controls (N=500,199).
All Howard-dependent LDSC and LAVA calculations were rerun. The unaffected
FinnGen depression x de Lange IBD LAVA output was retained, and BH-FDR was
recomputed across the complete corrected family of 110 tests.

The prior fixed-effect synthesis of the two reciprocal genome-wide estimates
was removed for a separate methodological reason: the pairings share FinnGen
participation and have non-zero sampling dependence. Current files report the
two estimates separately.

Current-tree hashes are recorded in `public_repository_sha256.txt`. The dated
correction protocol, input-identity report, before-after comparisons, and
execution notes are available under `reports/`.
