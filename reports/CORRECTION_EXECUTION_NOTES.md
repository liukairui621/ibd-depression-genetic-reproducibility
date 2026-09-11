# Correction execution notes

## Frozen archive

The pre-correction scripts, results, manuscript, reports, provenance, figures, submission package, and existing release archives were copied to `archive_pre_correction/`. File-level checksums are stored in `provenance/PRE_CORRECTION_SHA256SUMS.txt`. The correction protocol was frozen before corrected LDSC or LAVA results were inspected.

## Munging compatibility event

`munge_sumstats.py` completed successfully and wrote the corrected summary-statistics file. The first post-munging validator then failed because it was invoked with the LDSC Python 2 interpreter while containing Python 3 f-string syntax. No scientific output was altered. The validator was changed to the project Python 3 runtime and rerun against the completed file.

Validation then confirmed:

- 9,997,231 reference-template rows;
- 6,465,253 nonmissing summary statistics;
- constant nonmissing N = 500,199.

The corrected script and separate validator are preserved in `scripts/`, and both attempts are retained in `logs/`.

## Analysis boundaries

Only analyses affected by Howard sample metadata were rerun. The unaffected FinnGen depression x de Lange IBD LAVA output remains frozen; its BH q values are recalculated jointly with the corrected Howard x FinnGen tests. The fixed-effect synthesis of the two primary global correlations is removed independently of the corrected numerical results because shared FinnGen participation violates independence.
