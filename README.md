# IBD-depression genetic reproducibility

Code, protocols, results and provenance for **Cross-cohort reproducibility and
molecular resolution of genetic overlap between inflammatory bowel disease and
depression**. Submission revision: 15 September 2026.

The reciprocal LDSC estimates were small and positive. No LAVA component met the
same-block, same-direction, dual-significance rule. Three non-MHC PLACO blocks
reproduced, but molecular resolution was uneven and no gene had support across
all four GWAS. Super-enhancer annotations carried positive covariance in both
pairs. Enrichment ratios are descriptive; departure from 1 was not tested.

## Source-data reproduction

Use Linux. The supported entry point stages archived drivers under an arbitrary,
empty work directory. All six core inputs are rebuilt from source GWAS; private
preprocessed files are not needed.

1. Obtain GWAS and reference resources using
   `reproducibility/SOURCE_AND_ENVIRONMENT_GUIDE.md` and `manifests/`.
2. Adapt `reproducibility/assets.example.json` to your downloaded files.
3. Stage a new workspace:

```bash
python3 reproducibility/prepare_workspace.py --work-root /work/ibd_reproduction --assets assets.json
```

Staging ends with `reproducibility/preflight_downstream.py`, which resolves all
stage entry points and interpreters, checks extension output directories, scans
for unresolved archived paths, and validates shell syntax without starting an
analysis. The resulting record is written to
`provenance/PORTABLE_PREFLIGHT.json`.

4. Execute stages in order:

```bash
PROJECT=/work/ibd_reproduction/IBD/20_Reproducibility_Ladder
python3 reproducibility/run_stage.py --project "$PROJECT" rebuild_inputs
python3 reproducibility/verify_rebuilt_core.py --project "$PROJECT"
python3 reproducibility/run_stage.py --project "$PROJECT" core
python3 reproducibility/run_stage.py --project "$PROJECT" lava
python3 reproducibility/run_stage.py --project "$PROJECT" phenotype_extension
python3 reproducibility/run_stage.py --project "$PROJECT" placo
python3 reproducibility/run_stage.py --project "$PROJECT" finemap_qtl
python3 reproducibility/run_stage.py --project "$PROJECT" annotation
```

The full analysis requires substantial time, RAM and disk space. Exact commands,
logs and return codes are retained. Staging records path substitutions and input
assets. Do not run the archived absolute-path drivers directly from the repository.

The historical `w_hm3_alleles.snplist` filename refers to the full European
reference allele list, regenerated from the BIM files. HapMap3 restriction is
imposed by the LDSC regression LD-score and weight files.

## Verify published tables and redraw figures

```bash
python3 reproducibility/verify_published_results.py
PLOS_FIGURE_OUT=/work/figure_check Rscript scripts/27_prepare_plos_figures.R
python scripts/28_refine_submission_figures.py
```

The R script redraws the two main statistical figures. The Python refinement
script then applies the submission layout for Fig 1, S1 Fig, and S2 Fig using
dedicated exterior legend regions and fixed panel-label anchors. It changes
layout only and does not recompute statistics. See
`reproducibility/VERIFICATION_STATUS.md` for the checks actually executed for
this revision and the remaining scope.

## Contents

- `results/derived/`: global, local and phenotype-sensitivity outputs.
- `extensions/`: PLACO, fine-mapping/QTL and annotation analyses.
- `scripts/`: archived statistical and figure code.
- `reproducibility/`: portable staging, source preprocessing and result checks.
- `manifests/`, `provenance/`, `reports/`: sources and correction records.

Historical protocols and numeric outputs remain unchanged. The legacy
`replicated_concentration` table label means significant covariance plus
enrichment point estimates greater than 1 in both pairs, not a formal enrichment
test. The fitted model used 84 binary annotations from the 97-column resource.

Howard data exclude 23andMe: 170,756 cases and 329,443 controls, N=500,199,
doi:10.7488/ds/2458. The sample-size correction history is retained. An absent
credible set is insufficient evidence of disease specificity.

## Data access

Source GWAS, LD panels, large fitted objects and provider-controlled QTL files
must be obtained from their providers and are not redistributed. Derived results
include positive and negative colocalisation tests. Live QTL queries may change;
dataset identifiers and query records identify the analyzed resources.

Target URL: https://github.com/liukairui621/ibd-depression-genetic-reproducibility

Code: MIT license. Original data retain their providers' terms.
Corresponding author: Youxing Huang (waiqike7@163.com).
