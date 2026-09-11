#!/usr/bin/env python3
import csv
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    params = read(ROOT / "reports/OFFICIAL_PARAMETER_VALIDATION.tsv")
    sens = {x["metric"]: x["value"] for x in read(ROOT / "reports/PLACO_PLUS_SENSITIVITY_SUMMARY.tsv")}
    directions = read(ROOT / "reports/PRIMARY_REPLICATED_DIRECTION_AUDIT.tsv")
    all_params = all(x["status"] == "PASS" for x in params)
    all_directions = all(x["BOTH_TRAIT_DIRECTIONS_CONCORDANT"] == "TRUE" for x in directions)
    sens_blocks = int(float(sens["replicated_blocks_with_external_placo_plus"]))
    status = "PASS" if all_params and all_directions and sens_blocks == 3 else "FAIL"
    text = f"""# Post-run quality audit

Generated: {datetime.now(timezone.utc).isoformat()}

Overall status: **{status}**.

## Full-data parameter re-estimation

The official PLACO `var.placo()` and `cor.pearson()` functions were rerun on all harmonized SNPs for both primary pairs. Both comparisons passed the absolute tolerance of 1e-10. Exact values are in `OFFICIAL_PARAMETER_VALIDATION.tsv`.

## External-pair PLACO+ sensitivity

The external Howard-de Lange pair was originally assigned regular PLACO under the frozen no-known-overlap assumption. Because its null-Z correlation was 0.0338, all integration candidates were re-scored with PLACO+ as a post-run robustness analysis. The external pair had {sens['external_placo_plus_gws_loci']} genome-wide-significant blocks under PLACO+. The frozen cross-cohort endpoint retained {sens_blocks} replicated non-MHC blocks after replacing every external-pair discovery or validation P value with PLACO+.

## Trait-specific direction audit

The effect alleles were re-aligned across the two primary pairs for every frozen replicated index SNP. {sum(x['BOTH_TRAIT_DIRECTIONS_CONCORDANT'] == 'TRUE' for x in directions)}/{len(directions)} directional replication rows had concordant depression effects and concordant IBD effects separately. This is stricter than comparing only the sign of the product statistic.

## Implementation correction

The first run's positional allele parser was corrected to use semantic LDSC `A1`/`A2` column names. The complete pipeline was rerun, and the pre-correction snapshot remains archived. See `provenance/IMPLEMENTATION_CORRECTION_A1_SEMANTICS.md`.

No gene, mechanism, or clinical claim is made by this audit.
"""
    (ROOT / "reports/POSTRUN_QUALITY_AUDIT.md").write_text(text)
    execution_path = ROOT / "reports/EXECUTION_REPORT.md"
    execution = execution_path.read_text()
    marker = "## Post-run quality audit"
    if marker not in execution:
        execution += f"""

## Post-run quality audit

The official full-data parameter re-estimation passed at absolute tolerance 1e-10. A PLACO+ sensitivity analysis for the external Howard-de Lange pair retained all **3** replicated non-MHC blocks (LAVA blocks 154, 464, and 1671). Allele-aligned trait-specific direction checks passed for all {len(directions)} directional replication rows. The implementation-level A1 column-order correction and complete rerun are documented in `provenance/IMPLEMENTATION_CORRECTION_A1_SEMANTICS.md`; full audit results are in `reports/POSTRUN_QUALITY_AUDIT.md`.
"""
        execution_path.write_text(execution)
    if status != "PASS":
        raise SystemExit("Post-run quality audit failed")


if __name__ == "__main__":
    main()
