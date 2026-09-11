# Execution report: cross-cohort shared-locus analysis

Generated: 2026-08-17T17:03:33.318042+00:00

## Scope

This run tested variant-level depression-IBD pleiotropy with a cohort-separated primary comparison and two diagonal sensitivity pairs. It used the frozen protocol in `prespec/PRESPEC_SHARED_LOCUS_REPLICATION.md`. No manuscript, submitted figure, or submitted table was modified.

## Harmonization and model parameters

- `External_Howard_deLange`: 6,131,698 harmonized SNPs; VarZ=(1.5065, 1.4231); null-Z correlation=0.0338; method=placo.
- `FinnGen_DEP_IBD`: 6,607,911 harmonized SNPs; VarZ=(1.4675, 1.3448); null-Z correlation=0.0293; method=placo_plus.
- `Diagonal_Howard_FinnGenIBD`: 5,853,763 harmonized SNPs; VarZ=(1.5355, 1.3196); null-Z correlation=0.0134; method=placo.
- `Diagonal_FinnGenDEP_deLange`: 6,051,449 harmonized SNPs; VarZ=(1.4662, 1.4452); null-Z correlation=0.0195; method=placo.

All four inputs retained the common European-reference-matched variant set after numeric-Z and allele checks. Exact paths, schemas, counts, and SHA256 values are in `manifests/INPUT_MANIFEST.tsv`.

## Pair-level results

- `External_Howard_deLange`: minimum P=2.086e-17; 28 genome-wide-significant LD block(s), including 9 MHC block(s); 55 block(s) at P<1e-6.
- `FinnGen_DEP_IBD`: minimum P=1.087e-24; 23 genome-wide-significant LD block(s), including 12 MHC block(s); 42 block(s) at P<1e-6.
- `Diagonal_Howard_FinnGenIBD`: minimum P=3.497e-13; 21 genome-wide-significant LD block(s), including 11 MHC block(s); 45 block(s) at P<1e-6.
- `Diagonal_FinnGenDEP_deLange`: minimum P=2.097e-24; 38 genome-wide-significant LD block(s), including 10 MHC block(s); 60 block(s) at P<1e-6.

## Frozen primary replication endpoint

- External_to_FinnGen, block 908 (rs60689680): discovery P=2.086e-17, validation P=0.0427531501334467, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 1545 (rs10822050): discovery P=4.133e-16, validation P=0.20336152072983, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 151 (rs2477077): discovery P=1.291e-14, validation P=0.0801959206943639, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 952 (rs401763): discovery P=4.187e-14, validation P=2.68855286484466e-05, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 961 (rs148571474): discovery P=4.431e-14, validation P=0.857537851859491, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 2133 (rs2066844): discovery P=6.782e-14, validation P=0.141095311594814, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 960 (rs2523568): discovery P=1.430e-13, validation P=0.00264030954764962, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 967 (rs114595701): discovery P=1.703e-13, validation P=0.0497867058983603, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 57 (rs11209026): discovery P=3.671e-12, validation P=0.00111301151778551, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 2408 (rs6131010): discovery P=9.730e-12, validation P=0.00844795516023337, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 887 (rs12517950): discovery P=1.361e-11, validation P=0.0346936711894394, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 464 (rs148734725): discovery P=2.214e-11, validation P=NA, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 968 (rs57651384): discovery P=2.544e-11, validation P=0.16720667611248, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 2281 (rs7236656): discovery P=3.106e-11, validation P=0.894394312028222, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 2203 (rs8069176): discovery P=8.127e-11, validation P=0.0348672524283181, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 154 (rs169850): discovery P=2.066e-10, validation P=2.34690471350521e-07, threshold=0.00131578947368, same direction=TRUE, replicated=TRUE.
- External_to_FinnGen, block 951 (rs66462181): discovery P=5.405e-10, validation P=0.00475930692499271, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 268 (rs10191548): discovery P=1.359e-09, validation P=4.32911645747989e-06, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 1377 (rs36051895): discovery P=3.751e-09, validation P=0.00058372205945947, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 965 (rs9271549): discovery P=1.106e-08, validation P=4.53135155009745e-13, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 964 (rs2760994): discovery P=1.466e-08, validation P=4.03694308132342e-05, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 1875 (rs12427851): discovery P=1.507e-08, validation P=0.902828156132831, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 1192 (rs13234982): discovery P=1.559e-08, validation P=0.0304643829858558, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 62 (rs10889997): discovery P=1.958e-08, validation P=0.206157903001076, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 950 (rs34661691): discovery P=2.670e-08, validation P=0.114812585592055, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 998 (rs6455092): discovery P=2.773e-08, validation P=0.0249382485053917, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- External_to_FinnGen, block 159 (rs3024493): discovery P=2.862e-08, validation P=8.2674229045997e-13, threshold=0.00131578947368, same direction=FALSE, replicated=FALSE.
- External_to_FinnGen, block 248 (rs6715059): discovery P=4.951e-08, validation P=0.073973091115712, threshold=0.00131578947368, same direction=TRUE, replicated=FALSE.
- FinnGen_to_External, block 464 (rs9862080): discovery P=1.087e-24, validation P=9.18482525882017e-11, threshold=0.00227272727273, same direction=TRUE, replicated=TRUE.
- FinnGen_to_External, block 965 (rs9271580): discovery P=2.055e-14, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 465 (rs7613681): discovery P=3.189e-14, validation P=0.000642622537435382, threshold=0.00227272727273, same direction=TRUE, replicated=FALSE.
- FinnGen_to_External, block 962 (rs9268951): discovery P=5.271e-14, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 159 (rs3122605): discovery P=1.236e-13, validation P=3.96370766244866e-07, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 1100 (rs12672490): discovery P=1.086e-12, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 964 (rs28366327): discovery P=2.231e-12, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 956 (rs2516714): discovery P=2.697e-12, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 1579 (rs6584283): discovery P=3.319e-12, validation P=3.01703920806635e-07, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 963 (rs71549254): discovery P=4.596e-12, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 966 (rs3830060): discovery P=5.586e-12, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 955 (rs9260028): discovery P=4.615e-11, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 1671 (rs174581): discovery P=2.474e-10, validation P=1.81059597548045e-05, threshold=0.00227272727273, same direction=TRUE, replicated=TRUE.
- FinnGen_to_External, block 954 (rs2734975): discovery P=4.127e-10, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 154 (rs3861929): discovery P=8.324e-10, validation P=8.38194072992613e-08, threshold=0.00227272727273, same direction=TRUE, replicated=TRUE.
- FinnGen_to_External, block 967 (rs3104407): discovery P=1.510e-09, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 961 (rs492899): discovery P=1.843e-09, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 957 (rs2844637): discovery P=9.849e-09, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 1097 (rs798514): discovery P=1.067e-08, validation P=0.000732357548824927, threshold=0.00227272727273, same direction=TRUE, replicated=FALSE.
- FinnGen_to_External, block 960 (rs56044559): discovery P=1.709e-08, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 901 (rs34021572): discovery P=1.966e-08, validation P=0.409867688484451, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 887 (rs12719454): discovery P=2.614e-08, validation P=NA, threshold=0.00227272727273, same direction=FALSE, replicated=FALSE.
- FinnGen_to_External, block 2205 (rs12937642): discovery P=3.243e-08, validation P=0.11814229495631, threshold=0.00227272727273, same direction=TRUE, replicated=FALSE.

Primary replicated non-MHC blocks: **3**.

Decision: **GO**.

## Secondary observations

The two primary pair analyses shared 18 exploratory LD block(s) with lead PLACO P<1e-6 in both; 11 had the same sign of Z_depression x Z_IBD. These are descriptive convergence signals only and do not meet the prespecified replication definition. Results at all primary lead SNPs across the four pairings are retained in `results/loci/four_pair_lead_snp_support.tsv`.

## Interpretation boundary

The endpoint concerns statistical sharing at a SNP/LD-block level. It does not establish a causal gene, molecular mechanism, direction of causation, or clinical utility. Pair-specific or P<1e-6 signals are not promoted to candidate genes. The FinnGen within-cohort analysis used PLACO+ to account for correlated Z scores; the external pair used original PLACO under the no-known-overlap assumption. Exact-SNP replication is deliberately stringent and may miss replication expressed through a different proxy SNP in the same LD block.

## Reproducibility

- Official PLACO source commit and checksum are retained under `software/PLACO`.
- Full harmonized pair files, integration candidates, exact PLACO outputs, locus tables, and validation tests remain on the server under this project.
- Software versions and R session information are under `provenance/`.
- A project-wide SHA256 manifest is generated after pipeline completion.


## Post-run quality audit

The official full-data parameter re-estimation passed at absolute tolerance 1e-10. A PLACO+ sensitivity analysis for the external Howard-de Lange pair retained all **3** replicated non-MHC blocks (LAVA blocks 154, 464, and 1671). Allele-aligned trait-specific direction checks passed for all 4 directional replication rows. The implementation-level A1 column-order correction and complete rerun are documented in `provenance/IMPLEMENTATION_CORRECTION_A1_SEMANTICS.md`; full audit results are in `reports/POSTRUN_QUALITY_AUDIT.md`.
