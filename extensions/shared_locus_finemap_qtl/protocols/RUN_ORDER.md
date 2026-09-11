# Reproducible run order

All commands are run on the analysis server. The project root is:

```bash
export ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_finemap_qtl
cd "$ROOT"
```

The frozen rules must be read before execution:

```bash
less provenance/PRESPECIFIED_RULES.md
```

## 1. Harmonise GWAS, construct LD and fine-map

```bash
bash scripts/00_run_gwas_stage.sh
```

This wrapper runs, in order:

```bash
python3 scripts/01_prepare_gwas_inputs.py --root "$ROOT"
bash scripts/02_build_ld.sh "$ROOT"
Rscript scripts/03_run_gwas_susie_coloc.R "$ROOT"
```

## 2. Map positional candidates

```bash
python3 scripts/04_map_candidate_genes.py --root "$ROOT"
```

## 3. Retrieve and prepare local QTL data

```bash
software/qtl_env/bin/python scripts/05_retrieve_prepare_qtl.py --root "$ROOT"
```

The retrieval stage uses the final harmonised GWAS variant universe and the frozen molecular-trait eligibility gate. CEDAR HTTP failures are written to the retrieval manifest and do not stop analyzable GTEx/BLUEPRINT datasets.

## 4. Run QTL-GWAS colocalisation

```bash
Rscript scripts/06_run_qtl_coloc.R "$ROOT"
```

## 5. Query pQTL evidence and check LD separation

```bash
software/qtl_env/bin/python scripts/07_query_opentargets_pqtl.py --root "$ROOT"
bash scripts/07b_check_mst1_pqtl_ld.sh "$ROOT"
```

The first script saves the complete Open Targets JSON responses. The `07b` script runs the selected MST1 lead-variant LD checks with PLINK 1.9 against the same 1000 Genomes EUR reference; its complete and selected outputs are retained in `results/opentargets/ld_checks/` and `logs/07b_MST1_pqtl_ld.log`.

## 6. Summarise, plot and validate

```bash
software/qtl_env/bin/python scripts/08_summarize_and_plot.py --root "$ROOT"
software/qtl_env/bin/python scripts/09_validate_results.py --root "$ROOT"
```

## 7. Finalise provenance

```bash
bash scripts/10_finalize_archive.sh
```

This records software versions, creates the complete file inventory and regenerates and verifies the SHA256 manifest. Re-running live QTL or Open Targets queries at a later date may retrieve a newer upstream snapshot; archived raw responses identify the snapshot used here.
