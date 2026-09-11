# Annotation-stratified genetic covariance report

## Frozen decision

- Replicated annotation-level covariance: 2 of 10 annotations.
- Replicated concentration above SNP proportion: 2 of 10 annotations.
- Primary regression used HapMap3 weights excluding the extended MHC.

## Genome-wide base annotation

- Discovery (deLange IBD x FinnGen depression): covariance 0.0112784 (SE 0.0033681, P 0.0008122).
- Replication (FinnGen IBD x Howard public no-23andMe depression): covariance 0.0046862 (SE 0.00114757, P 4.435e-05).

## Annotation decisions

- Coding_UCSC: discovery covariance 0.000390989 (P 0.8409, BH-FDR 0.9454); replication covariance 0.000467287 (directional P 0.2381); not_discovery_fdr.
- Conserved_LindbladToh: discovery covariance 0.00399249 (P 0.2348, BH-FDR 0.587); replication covariance -0.000710236 (directional P 0.7609); not_discovery_fdr.
- DHS_Trynka: discovery covariance -0.00316378 (P 0.6754, BH-FDR 0.9454); replication covariance 0.00294587 (directional P 0.8782); not_discovery_fdr.
- Enhancer_Andersson: discovery covariance 8.81758e-05 (P 0.954, BH-FDR 0.954); replication covariance -7.93414e-06 (directional P 0.5062); not_discovery_fdr.
- H3K27ac_Hnisz: discovery covariance 0.00603464 (P 0.04335, BH-FDR 0.1445); replication covariance 0.00336679 (directional P 0.001148); not_discovery_fdr.
- H3K4me1_Trynka: discovery covariance -0.00481413 (P 0.3797, BH-FDR 0.7595); replication covariance 0.00511387 (directional P 0.9992); not_discovery_fdr.
- Promoter_UCSC: discovery covariance 0.000524979 (P 0.8509, BH-FDR 0.9454); replication covariance 0.000723819 (directional P 0.2295); not_discovery_fdr.
- SuperEnhancer_Hnisz: discovery covariance 0.00470035 (P 0.007955, BH-FDR 0.03978); replication covariance 0.00240057 (directional P 8.279e-05); replicated_concentration.
- TFBS_ENCODE: discovery covariance -0.00423861 (P 0.5151, BH-FDR 0.8585); replication covariance 0.00414964 (directional P 0.9886); not_discovery_fdr.
- Repressed_Hoffman: discovery covariance 0.0155385 (P 0.00346, BH-FDR 0.0346); replication covariance 0.00332551 (directional P 0.04549); replicated_concentration.

## Interpretation boundary

At least one frozen annotation carried reproducible covariance across reciprocal cohort pairs. This localizes covariance to a broad genomic context but does not identify a shared causal pathway, cell type, gene or treatment target.

The prior three-locus fine-mapping result remains the local adjudication: no cross-disease shared causal gene was identified.
