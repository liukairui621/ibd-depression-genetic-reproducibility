# Sources, software, and execution boundaries

## Input data

The repository does not redistribute source GWAS or LD panels. Download source
files under the providers' current terms, retain their identifiers and checksums,
and map their absolute locations in a copy of `assets.example.json`. The keys are
workspace destinations, not required names or locations on your own computer.

| Asset | Source and identity | Handling |
|---|---|---|
| Howard depression | Edinburgh DataShare https://doi.org/10.7488/ds/2458; public file excluding 23andMe | 170,756 cases + 329,443 controls, N=500,199; 8,483,301 source variants; broad depression meta-analysis, not strictly diagnosed MDD |
| de Lange IBD | GWAS Catalog GCST004131, PMID 28067908, harmonised summary statistics | Columns `hm_rsid`, `hm_effect_allele`, `hm_other_allele`, `hm_beta`, `p_value`; N=59,957 |
| FinnGen depression | Release 12, `F5_DEPRESSIO` | 59,333 cases, 434,831 controls |
| FinnGen IBD | Release 12, `IBD_STRICT` | 10,960 cases, 489,388 controls |
| FinnGen CD | Release 12; archived source alias `CD_STRICT2` | 2,489 cases, 497,622 controls; match the source release endpoint and checksum, not just the filename |
| FinnGen UC | Release 12; archived source alias `UC_STRICT2` | 7,220 cases, 492,160 controls; match the source release endpoint and checksum, not just the filename |
| UKB phenotype extension | Cai et al., https://doi.org/10.1038/s41588-020-0594-5 | Four Figshare file IDs, download URLs, MD5 values, and sample definitions in `manifests/extension_source_manifest.tsv`; per-variant NMISS is used |
| de Lange CD/UC extension | GCST004132 / GCST004133 | Harmonised files; URLs and total N are in the extension manifest |

FinnGen release access starts at https://www.finngen.fi/en/access_results.
Do not substitute a newer release or the larger Howard dataset containing
23andMe. The download scripts for the phenotype extension contain resumable
transfers and checksum checks. Core data are explicit user-supplied assets;
access conditions are not bypassed by this pipeline.

## Reference assets

Use the 1000 Genomes Phase 3 EUR resources corresponding to the LDSC reference
distribution, not a different ancestry or a GRCh38 replacement. Reference
resource descriptions are available at
https://github.com/bulik/ldsc/wiki and https://github.com/josefin-werme/LAVA.
The configured PLINK directory must contain `1000G.EUR.QC.1` through
`1000G.EUR.QC.22` (each with `.bed`, `.bim`, `.fam`) and the merged
`1000G.EUR.QC.ALL` files for fine-mapping. To create the latter, merge the same
chromosome files in numeric order with PLINK `--merge-list`, `--keep-allele-order`,
and `--make-bed`; preserve individual order and report any merge failure.
For chromosome 1 through 22, the staged drivers expect these prefixes:

- `LD_EUR/LDscore/LDscore.` for univariate/bivariate LDSC LD scores.
- `LD_EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.` for weights.
- `LD_EUR/1000G_Phase3_frq/1000G.EUR.QC.` for reference frequencies.
- `SLDSC/baselineLD.` for baseline-LD v2.2 LD scores and annotation metadata.

These prefixes are relative to `00_RawData/Reference/`. Standard provider
archives can use different directory/prefix names; create explicit aliases
to these names without changing file contents. Retain the accompanying
`.l2.M`, `.l2.M_5_50`, `.l2.ldscore.gz`, `.annot.gz`, and frequency files
required by the relevant resource. Do not mix reference versions.

`resources/LAVA_blocks_fixed.txt` is the frozen 2,495-block coordinate file used
by the analysis. It is provided as a small coordinate table; genotype data are
not included. Blocks follow the European LD partition described by Berisa and
Pickrell (2016), as used by LAVA. Its SHA256 is included in the revision manifest.

The legacy file `w_hm3_alleles.snplist` is NOT a HapMap3-only list. The archived
file has 9,997,231 reference variants. `rebuild_core_inputs.py` reconstructs it
from chromosome BIM columns 2, 5, and 6, with header `SNP A1 A2`, numeric
chromosome order, and one whitespace-delimited row per variant. LDSC's regression SNP
set is subsequently determined by its HapMap3 LD scores and weights. PLACO
uses the broader reference-matched set. Do not replace the legacy allele file
with the smaller public HapMap3 list when reproducing these results.

## Recorded environment

Use Linux with bash, curl, gzip, git, bgzip, tabix, PLINK, and Rscript on PATH.
Scripts are archived with their original absolute prefixes. The workspace
stager rewrites only paths, records every replacement, and never copies
precomputed statistical results into the analysis workspace.

| Component | Recorded version |
|---|---|
| LDSC | 1.0.1; commit `2fdeeb3b44379408794154993dbd6101b8946b7e` |
| LDSC Python | 2.7.18; numpy 1.16.5, pandas 0.24.2, scipy 1.2.1 |
| R | 4.3.3 |
| LAVA | 0.1.5 |
| GenomicSEM | 0.0.5; commit `0a63ac0ea01b61d28bd17e4a204e0fa561ce5040` |
| PLACO | v0.2.0; https://github.com/RayDebashree/PLACO; commit `3ba3cae1d323ad117fb4540e620bcefa79f70663` |
| PLINK | 1.9.0-b.7.7, 22 October 2024 |
| susieR / coloc | 0.14.2 / 5.2.3 |
| Molecular Python | 3.12.3; numpy 2.1.3, pandas 2.2.3, matplotlib 3.9.2, pysam 0.23.3, requests 2.32.5 |

The detailed R dependency versions are in the archived session-information
files. R packages include data.table, LAVA, GenomicSEM, susieR, coloc, ggplot2,
patchwork, and scales. `assets.example.json` points to existing software and an
R library; it is not an environment lockfile or a claim of a tested fresh
installation. Keep the Python 2 LDSC environment separate from Python 3.
Capture `sessionInfo()`, Python package versions, executable versions, and git
revisions before comparing results. Seeds are 42 for the core analysis and
20260817 for fine-mapping. S-LDSC uses 200 deterministic jackknife blocks.

## Run order and network-dependent branches

Follow the staged commands in the repository README. `rebuild_inputs` must
complete before `core`, `lava`, or any extension. `placo` must complete before
`finemap_qtl`; annotation uses the core munged inputs directly. The phenotype
extension downloads six additional source files. Fine-mapping retrieves
GENCODE/eQTL Catalogue annotations and queries Open Targets. Those live
resources can change independently of the code: a later query is not a
bit-identical recreation of the archived API responses. The public derived
QTL tables and query manifests preserve the evidence underlying this paper.

Historical scripts that only assembled manuscript reports are retained as
provenance, not required as statistical stages. `run_stage.py` records the
actual command, exit status, and log for each explicitly requested stage.
`VERIFICATION_STATUS.md` states which stages were actually re-executed for
this revision; syntax checks alone are not evidence of end-to-end execution.
