# Implementation amendment 001: annotation-name normalization

Amendment scope: output extraction only. Scientific inputs, model, annotations, multiple-testing family, decision thresholds and stopping rule are unchanged.

## Trigger

The discovery `s_ldsc` computation completed successfully and produced:

- Raw RDS SHA256: `26b14a443e16a02352a29c330ed14d33a39ef270cbe533e7696f025c1b5dbaa3`
- Completed log timestamp: `2026-08-24 09:47:26`

The extraction step then stopped because the installed GenomicSEM version appends `L2` to annotation names (`baseL2`, `Coding_UCSCL2`, etc.), whereas the frozen extractor expected `base`.

## Mechanical correction

`scripts/01_run_sldsc_pair.R` was changed to:

1. remove the terminal `L2` suffix when exporting annotation names; and
2. reuse an already completed raw RDS rather than rerunning the same S-LDSC computation.

Original script SHA256:

`7db44b0edb76bf628cb8df8c5d03f9659525ba60b901b0cadfbfd6977b45b6d2`

Amended script SHA256:

`e35305eafcab1b1aca2fc9e79561f5a4823a62d0dc2786c3b91e18237de7cd6e`

## Annotation-count clarification

The baseline-LD files contain 97 annotation columns. The frozen call used `exclude_cont=TRUE`; therefore the installed function retained 84 binary/non-continuous annotation matrices for estimation. This is the prespecified behavior, not a post-result exclusion. Formal inference remains restricted to the same ten frozen binary annotations.

