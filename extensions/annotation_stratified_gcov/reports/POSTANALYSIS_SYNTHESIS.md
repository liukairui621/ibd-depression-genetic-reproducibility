# Post-analysis synthesis with local fine-mapping adjudication

## What changed

The new analysis adds a functional-annotation layer without altering the existing locus-level verdict.

- Genome-wide covariance was positive in both reciprocal cohort pairs: 0.01128 (SE 0.00337, P = 8.12e-4) for deLange IBD x FinnGen depression and 0.00469 (SE 0.00115, P = 4.44e-5) for FinnGen IBD x Howard public no-23andMe depression.
- `SuperEnhancer_Hnisz` met the frozen discovery-replication rule: discovery covariance 0.00470 (BH-FDR = 0.0398), replication covariance 0.00240 (directional P = 8.28e-5), with covariance-enrichment point estimates of 2.49 and 3.06.
- `Repressed_Hoffman`, the prespecified negative/comparator annotation, also met the literal rule: discovery covariance 0.01554 (BH-FDR = 0.0346) and replication covariance 0.00333 (directional P = 0.0455). This is a marginal result because the replication two-sided P value is 0.091 and the annotation is broad (46.0% of SNPs). It should not be converted into a mechanistic claim.
- `H3K27ac_Hnisz` was directionally consistent and nominal in discovery but did not pass the discovery multiple-testing threshold (BH-FDR = 0.145).
- The other seven frozen annotations did not replicate under the prespecified rule; several changed direction.

## Relationship to the three fine-mapped blocks

The local analysis remains unchanged:

- block154 supported an IBD component near GPR25, not a replicated depression component;
- block464 supported an IBD component near MST1, while IBD and depression protein signals were distinct low-LD signals;
- block1671 yielded FinnGen-depression molecular candidates near FADS1/TMEM258 without reciprocal depression or IBD support;
- no gene was supported as a shared causal gene across both IBD and both depression GWAS.

The global and local layers therefore answer different questions. The annotation analysis suggests that a portion of the small cross-trait covariance is reproducibly carried by broad regulatory sequence, most clearly super-enhancers. Fine-mapping shows that this broad covariance does not resolve into the same causal gene or molecular signal in the three previously highlighted local blocks.

## Manuscript-level conclusion

The appropriate combined conclusion is:

> IBD and depression share a small, reproducible genome-wide genetic covariance with evidence of concentration in super-enhancer annotations, but the examined local-correlation blocks decompose into disease- or cohort-specific molecular components rather than replicated shared causal genes.

This result weakens a purely nonspecific-null interpretation of the global correlation, but it does not rescue a shared-gene, shared-pathway, causal or clinical-target claim. The repressed-chromatin comparator and the negative replication intercept flag require the super-enhancer result to remain a bounded secondary finding.

## Scope boundary

These baseline-LD annotations are broad genomic contexts, not brain, gut or immune cell-type annotations. The analysis cannot identify an organ, cell population, pathway, biomarker, risk-stratification rule or treatment target. Under the frozen stopping rule, it does not trigger another post-hoc search for shared genes or pathways.

