# Reporting and reproducibility amendment 5 September 2026

No frozen statistical output or analysis threshold was changed in this revision.
The original server analysis directories and the 31 August submission package
remain the historical records. This public revision corrects executable report
templates and adds reproducibility infrastructure.

1. **Model specification.** The annotation freeze described all 97 baseline-LD
columns as jointly fitted. The executed GenomicSEM call used
`exclude_cont=TRUE`, retaining 84 binary annotations. This is a discrepancy
between the original protocol description and the implemented model, now
reported explicitly. The frozen ten-annotation testing family is unchanged.
The original protocol text is retained for chronology, not silently rewritten.

2. **Inference.** Annotation P values test covariance against zero. Enrichment
ratios are descriptive point estimates; no standard errors or tests against
enrichment of 1 were calculated. The legacy `replicated_concentration` flag
combines covariance gates with enrichment point estimates above 1; its label
must not be interpreted as significant enrichment.

3. **Molecular resolution.** No candidate had support across all four GWAS.
Absence of a SuSiE credible set does not establish disease specificity. Block
464 included an IBD/FinnGen-depression colocalisation but no Howard depression
credible set. Detailed positive and nonpositive molecular-QTL output tables
are now included in `extensions/shared_locus_finemap_qtl/results/full/`.

4. **SNP scope.** The archived `w_hm3_alleles.snplist` name was misleading.
It contains 9,997,231 full EUR reference variants, reconstructed from the
chromosome BIM files. Munging and PLACO used this broader set; LDSC regression
uses HapMap3 LD scores and weights. FinnGen matching used rsID and alleles
against the GRCh37 reference, not an undocumented coordinate-liftover step.

5. **Raw reconstruction.** Six fresh core munged inputs were generated without
reading historical preprocessed files. Their decompressed SHA256 hashes match
the six frozen analysis inputs, including the corrected no-23andMe Howard
file. The initial audit encountered the old Howard symlink in the original
project; that historical link is not the corrected comparison baseline.
The new workspace produces its own correct input and internal aliases.

6. **Verification scope.** Core LDSC was rerun from the rebuilt inputs:
six heritability rows and 15 genetic-correlation rows matched the published
exports. BH correction and key counts were independently checked from the
local, phenotype-extension, and annotation exports. The full LAVA, PLACO,
SuSiE/QTL, GenomicSEM-extension, and S-LDSC branches were not all re-executed.
See `reproducibility/VERIFICATION_STATUS.md` for the precise boundary.

Archived checksum manifests refer to their original snapshots. The current
public-file manifest and SHA256 list identify this revision; they do not
replace the historical manifests. This directory is prepared for publication,
not evidence that the target GitHub URL has already been made public.
