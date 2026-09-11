# UKB depression-definition and CD/UC extension

## Question and design

This extension tests whether variability in IBD-depression genetic correlation is more closely associated with depression phenotype definition or IBD subtype. Four UK Biobank depression definitions were crossed with IBD, CD, and UC from de Lange 2017 and FinnGen R12. The UKB definitions share participants and are phenotype-sensitivity analyses, not independent replications.

## Layer-0 QC

| Trait | h2_obs | SE | z | Intercept | Mean chi2 | h2 z>=4 |
|---|---:|---:|---:|---:|---:|:---:|
| CD_FinnGen_R12 | 0.005 | 0.001 | 4.17 | 1.0499 | 1.114 | yes |
| CD_deLange2017 | 0.303 | 0.035 | 8.53 | 1.1474 | 1.447 | yes |
| DEP_FinnGen_R12 | 0.031 | 0.002 | 16.21 | 1.1737 | 1.518 | yes |
| IBD_FinnGen_R12 | 0.018 | 0.002 | 8.05 | 1.1510 | 1.365 | yes |
| IBD_deLange2017 | 0.217 | 0.023 | 9.62 | 1.1794 | 1.497 | yes |
| MDD_Howard2019 | 0.030 | 0.001 | 22.85 | 1.0594 | 1.589 | yes |
| UC_FinnGen_R12 | 0.014 | 0.002 | 7.61 | 1.1102 | 1.276 | yes |
| UC_deLange2017 | 0.166 | 0.020 | 8.23 | 1.1222 | 1.309 | yes |
| UKB_GPpsy | 0.045 | 0.003 | 17.19 | 1.0569 | 1.385 | yes |
| UKB_ICD10Dep | 0.015 | 0.002 | 7.84 | 0.9941 | 1.064 | yes |
| UKB_LifetimeMDD | 0.058 | 0.008 | 7.46 | 1.0249 | 1.115 | yes |
| UKB_MDDRecur | 0.060 | 0.008 | 7.46 | 1.0208 | 1.101 | yes |

## Crossed UKB definition by IBD phenotype matrix

Cells are rg (SE); a dagger marks |genetic-covariance intercept z| >= 2.

| Depression definition | IBD_deLange2017 | CD_deLange2017 | UC_deLange2017 | IBD_FinnGen_R12 | CD_FinnGen_R12 | UC_FinnGen_R12 |
|---|---:|---:|---:|---:|---:|---:|
| UKB_LifetimeMDD | 0.137 (0.060) | 0.163 (0.059) | 0.080 (0.069) | 0.230 (0.057)† | 0.154 (0.094) | 0.229 (0.067)† |
| UKB_MDDRecur | 0.147 (0.061) | 0.162 (0.060) | 0.105 (0.075) | 0.220 (0.058)† | 0.191 (0.099) | 0.209 (0.070)† |
| UKB_GPpsy | 0.124 (0.035) | 0.121 (0.035) | 0.107 (0.043) | 0.125 (0.037)† | 0.111 (0.058) | 0.099 (0.044)† |
| UKB_ICD10Dep | 0.167 (0.061) | 0.202 (0.059) | 0.087 (0.077) | 0.141 (0.067) | 0.045 (0.110) | 0.107 (0.067) |

All 24 estimates were positive, ranging from 0.045 (UKB_ICD10Dep x CD_FinnGen_R12) to 0.230 (UKB_LifetimeMDD x IBD_FinnGen_R12); 15/24 passed BH-FDR < 0.05. These counts describe signal detection and are not used to attribute instability because precision differs across cells.

## Attribution of instability

Classification: **mixed_or_indeterminate**.

The 95% interval for the adjusted dispersion ratio includes 1, or interaction/leave-one-out sensitivity prevents a stable attribution.

The adjusted between-definition SD was 0.025; the adjusted between-subtype SD was 0.020. Their ratio was 1.258 (simulation 95% CI 0.379-6.214; Pr[ratio>1]=0.753).

The point estimate therefore leaned toward greater variation across depression definitions, but the interval was wide and crossed 1. Excluding the broad GP-consultation definition moved the ratio to 0.985, so this tendency was not leave-one-definition robust.

IBD subtype showed an omnibus difference (P=0.005), whereas depression definition did not (P=0.470). However, CD exceeded UC for 4/4 UKB definitions in de Lange, while CD was lower than UC for 3/4 in FinnGen. This opposite cross-cohort direction and the absence of FDR-significant CD-versus-UC contrasts prevent a subtype-dominant interpretation.

Taken together, the extension supports a consistently positive IBD-depression genetic-correlation direction, but it does not support assigning the observed magnitude instability primarily to either depression definition or IBD subtype.

### Omnibus and interaction tests

| Test | Statistic | df | P |
|---|---:|---:|---:|
| depression_definition | 2.530 | 3 | 0.470 |
| ibd_subtype | 10.546 | 2 | 0.005 |
| ibd_cohort | 3.93e-04 | 1 | 0.984 |
| depression_definition_x_ibd_subtype | 6.069 | 6 | 0.416 |
| ibd_subtype_x_ibd_cohort | 1.721 | 2 | 0.423 |

### Leave-one-definition sensitivity

| Excluded UKB definition | Definition SD | Subtype SD | Ratio |
|---|---:|---:|---:|
| UKB_LifetimeMDD | 0.027 | 0.020 | 1.346 |
| UKB_MDDRecur | 0.024 | 0.020 | 1.205 |
| UKB_GPpsy | 0.020 | 0.020 | 0.985 |
| UKB_ICD10Dep | 0.029 | 0.020 | 1.448 |

### Leave-one-subtype sensitivity

| Excluded subtype | Definition SD | Retained-subtype SD | Ratio |
|---|---:|---:|---:|
| IBD | 0.025 | 0.026 | 0.978 |
| CD | 0.025 | 0.023 | 1.092 |
| UC | 0.025 | 0.003 | 9.320 |

### CD-versus-UC direction across UKB definitions

| IBD cohort | CD > UC | CD < UC | Prespecified direction rule |
|---|---:|---:|:---:|
| deLange2017 | 4 | 0 | pass |
| FinnGen_R12 | 1 | 3 | pass |

## Contrast summary

Sampling-covariance-aware contrasts identified 0 FDR-significant pairwise differences among UKB depression definitions and 0 FDR-significant CD-versus-UC differences. 6 of 24 primary extension cells had an absolute genetic-covariance intercept z-score of at least 2 and are flagged rather than removed.

## Interpretation boundaries

This analysis attributes variability in observed LDSC rg estimates. It does not establish which disease causes the other, and it does not partition biological mechanisms. The UKB definitions are correlated because they use overlapping participants. The multivariable LDSC sampling covariance matrix was therefore used for inferential contrasts. Significance counts alone were not used to classify stability.

## Reproducibility files

- Prespecified protocol: `provenance/extension_prespecified_protocol.md`
- Pairwise LDSC output: `results/global_extension/`
- GenomicSEM S and V matrices: `results/global_extension/genomicsem/`
- Attribution models and contrasts: `results/global_extension/attribution/`
- Figure panels: `figures/extension/`
- Final checksums: `provenance/extension_final_sha256.tsv`
