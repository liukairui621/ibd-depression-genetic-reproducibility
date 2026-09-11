#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")

INPUT_META = [
    ("Howard_MDD_no23andMe", "/root/IBD/20_Reproducibility_Ladder/corrections/20260806_howard_no23andme/work/munged/MDD_Howard2019_no23andMe.sumstats.gz", "GRCh37", 170756, 329443, 500199),
    ("deLange2017_IBD", "/root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_deLange2017.sumstats.gz", "GRCh37", 25042, 34915, 59957),
    ("FinnGen_R12_DEP", "/root/IBD/20_Reproducibility_Ladder/data/munged_core/DEP_FinnGen_R12.sumstats.gz", "GRCh38_to_rsid_HapMap3", 59333, 434831, 494164),
    ("FinnGen_R12_IBD_STRICT", "/root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_FinnGen_R12.sumstats.gz", "GRCh38_to_rsid_HapMap3", 10960, 489388, 500348),
]


def rows(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_rows(path, fields, data):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def as_bool(value):
    return str(value).strip().upper() == "TRUE"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_numeric_rows(path):
    total = numeric = 0
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        z_idx = header.index("Z")
        for line in handle:
            total += 1
            fields = line.rstrip("\n\r").split("\t")
            try:
                float(fields[z_idx])
                numeric += 1
            except (ValueError, IndexError):
                pass
    return "|".join(header), total, numeric


def main():
    report_dir = ROOT / "reports"
    manifest_dir = ROOT / "manifests"
    report_dir.mkdir(exist_ok=True)
    manifest_dir.mkdir(exist_ok=True)

    scored = list(rows(ROOT / "results/loci/replication_tests_scored.tsv"))
    final = []
    for row in scored:
        n_loci = int(row["n_nonmhc_discovery_loci"])
        threshold = 0.025 / n_loci if n_loci else None
        p_val = float(row["VALIDATION_P_PLACO"]) if row["VALIDATION_P_PLACO"] else None
        p_dep = float(row["VALIDATION_P_DEP"]) if row["VALIDATION_P_DEP"] else None
        p_ibd = float(row["VALIDATION_P_IBD"]) if row["VALIDATION_P_IBD"] else None
        replicated = bool(
            row["MHC"] != "TRUE"
            and threshold is not None
            and as_bool(row["VALIDATION_FOUND"])
            and p_val is not None and p_val < threshold
            and p_dep is not None and p_dep < 0.05
            and p_ibd is not None and p_ibd < 0.05
            and row["DISCOVERY_SIGN"] == row["VALIDATION_SIGN"]
        )
        row["VALIDATION_BONF_THRESHOLD"] = "" if threshold is None else f"{threshold:.12g}"
        row["BOTH_VALIDATION_MARGINAL_P_LT_0.05"] = str(p_dep is not None and p_ibd is not None and p_dep < 0.05 and p_ibd < 0.05).upper()
        row["DIRECTION_CONCORDANT"] = str(row["DISCOVERY_SIGN"] == row["VALIDATION_SIGN"] and bool(row["VALIDATION_SIGN"])).upper()
        row["PRIMARY_REPLICATED"] = str(replicated).upper()
        final.append(row)

    final_fields = list(final[0]) if final else [
        "direction", "discovery_pair", "validation_pair", "n_nonmhc_discovery_loci", "LOC", "CHR",
        "BLOCK_START", "BLOCK_STOP", "MHC", "LEAD_SNP", "LEAD_BP", "DISCOVERY_P_PLACO",
        "DISCOVERY_Q_BH", "DISCOVERY_Z_PRODUCT", "DISCOVERY_SIGN", "VALIDATION_FOUND", "VALIDATION_EA",
        "VALIDATION_OA", "VALIDATION_Z_DEP", "VALIDATION_Z_IBD", "VALIDATION_P_DEP", "VALIDATION_P_IBD",
        "VALIDATION_Z_PRODUCT", "VALIDATION_SIGN", "VALIDATION_METHOD", "VALIDATION_VAR_Z_DEP",
        "VALIDATION_VAR_Z_IBD", "VALIDATION_COR_Z", "VALIDATION_P_PLACO", "VALIDATION_BONF_THRESHOLD",
        "BOTH_VALIDATION_MARGINAL_P_LT_0.05", "DIRECTION_CONCORDANT", "PRIMARY_REPLICATED"
    ]
    write_rows(ROOT / "results/loci/cross_cohort_replication.tsv", final_fields, final)

    replicated_blocks = sorted({int(x["LOC"]) for x in final if x["PRIMARY_REPLICATED"] == "TRUE"})
    decision = "GO" if replicated_blocks else "NO-GO"
    (report_dir / "FINAL_DECISION.txt").write_text(
        f"{decision}\nreplicated_non_MHC_blocks={len(replicated_blocks)}\n"
    )

    harmonization = list(rows(ROOT / "results/raw/pair_harmonization_stats.tsv"))
    placo_summary = list(rows(ROOT / "results/raw/placo_run_summary.tsv"))
    pair_summary = []
    for item in placo_summary:
        pair = item["pair_id"]
        locus_data = list(rows(ROOT / f"results/loci/{pair}.loci_p1e6.tsv"))
        pair_summary.append({
            **item,
            "loci_p_lt_1e6": sum(float(x["LEAD_P_PLACO"]) < 1e-6 for x in locus_data),
            "gws_loci": sum(x["GWS_LOCUS"] == "TRUE" for x in locus_data),
            "gws_nonmhc_loci": sum(x["GWS_LOCUS"] == "TRUE" and x["MHC"] != "TRUE" for x in locus_data),
            "gws_mhc_loci": sum(x["GWS_LOCUS"] == "TRUE" and x["MHC"] == "TRUE" for x in locus_data),
        })
    pair_fields = list(pair_summary[0])
    write_rows(report_dir / "PAIR_LEVEL_SUMMARY.tsv", pair_fields, pair_summary)

    convergence = list(rows(ROOT / "results/loci/exploratory_same_block_convergence.tsv"))
    same_sign_convergence = sum(as_bool(x["SAME_SIGN"]) for x in convergence)

    manifest = []
    for trait, path_text, build, cases, controls, total_n in INPUT_META:
        path = Path(path_text)
        schema, total_rows, numeric_rows = count_numeric_rows(path)
        manifest.append({
            "trait_id": trait, "input_path": str(path), "real_path": os.path.realpath(path),
            "sha256": sha256(path), "size_bytes": path.stat().st_size, "schema": schema,
            "rows_after_header": total_rows, "numeric_z_rows": numeric_rows, "build_or_mapping": build,
            "n_cases": cases, "n_controls": controls, "n_total": total_n,
        })
    write_rows(manifest_dir / "INPUT_MANIFEST.tsv", list(manifest[0]), manifest)

    metric_rows = [
        {"metric": "primary_replicated_non_MHC_blocks", "value": len(replicated_blocks)},
        {"metric": "primary_replication_decision", "value": decision},
        {"metric": "exploratory_same_block_p_lt_1e6", "value": len(convergence)},
        {"metric": "exploratory_same_block_same_sign", "value": same_sign_convergence},
    ]
    for item in pair_summary:
        metric_rows.extend([
            {"metric": f"{item['pair_id']}_gws_loci", "value": item["gws_loci"]},
            {"metric": f"{item['pair_id']}_loci_p_lt_1e6", "value": item["loci_p_lt_1e6"]},
            {"metric": f"{item['pair_id']}_minimum_p", "value": item["minimum_p"]},
        ])
    write_rows(report_dir / "RESULTS_SUMMARY.tsv", ["metric", "value"], metric_rows)

    harm_lines = []
    for row in harmonization:
        harm_lines.append(
            f"- `{row['pair_id']}`: {int(row['retained']):,} harmonized SNPs; "
            f"VarZ=({float(row['var_z_dep']):.4f}, {float(row['var_z_ibd']):.4f}); "
            f"null-Z correlation={float(row['cor_z']):.4f}; method={row['method']}."
        )
    result_lines = []
    for row in pair_summary:
        min_p = float(row["minimum_p"]) if row["minimum_p"] else float("nan")
        result_lines.append(
            f"- `{row['pair_id']}`: minimum P={min_p:.3e}; "
            f"{row['gws_loci']} genome-wide-significant LD block(s), including {row['gws_mhc_loci']} MHC block(s); "
            f"{row['loci_p_lt_1e6']} block(s) at P<1e-6."
        )

    if final:
        replication_lines = []
        for row in final:
            validation_p = row["VALIDATION_P_PLACO"] or "NA"
            replication_lines.append(
                f"- {row['direction']}, block {row['LOC']} ({row['LEAD_SNP']}): discovery P={float(row['DISCOVERY_P_PLACO']):.3e}, "
                f"validation P={validation_p}, threshold={row['VALIDATION_BONF_THRESHOLD'] or 'NA'}, "
                f"same direction={row['DIRECTION_CONCORDANT']}, replicated={row['PRIMARY_REPLICATED']}."
            )
    else:
        replication_lines = ["- Neither primary pair produced a genome-wide-significant index locus requiring formal validation."]

    report = f"""# Execution report: cross-cohort shared-locus analysis

Generated: {datetime.now(timezone.utc).isoformat()}

## Scope

This run tested variant-level depression-IBD pleiotropy with a cohort-separated primary comparison and two diagonal sensitivity pairs. It used the frozen protocol in `prespec/PRESPEC_SHARED_LOCUS_REPLICATION.md`. No manuscript, submitted figure, or submitted table was modified.

## Harmonization and model parameters

{chr(10).join(harm_lines)}

All four inputs retained the common European-reference-matched variant set after numeric-Z and allele checks. Exact paths, schemas, counts, and SHA256 values are in `manifests/INPUT_MANIFEST.tsv`.

## Pair-level results

{chr(10).join(result_lines)}

## Frozen primary replication endpoint

{chr(10).join(replication_lines)}

Primary replicated non-MHC blocks: **{len(replicated_blocks)}**.

Decision: **{decision}**.

## Secondary observations

The two primary pair analyses shared {len(convergence)} exploratory LD block(s) with lead PLACO P<1e-6 in both; {same_sign_convergence} had the same sign of Z_depression x Z_IBD. These are descriptive convergence signals only and do not meet the prespecified replication definition. Results at all primary lead SNPs across the four pairings are retained in `results/loci/four_pair_lead_snp_support.tsv`.

## Interpretation boundary

The endpoint concerns statistical sharing at a SNP/LD-block level. It does not establish a causal gene, molecular mechanism, direction of causation, or clinical utility. Pair-specific or P<1e-6 signals are not promoted to candidate genes. The FinnGen within-cohort analysis used PLACO+ to account for correlated Z scores; the external pair used original PLACO under the no-known-overlap assumption. Exact-SNP replication is deliberately stringent and may miss replication expressed through a different proxy SNP in the same LD block.

## Reproducibility

- Official PLACO source commit and checksum are retained under `software/PLACO`.
- Full harmonized pair files, integration candidates, exact PLACO outputs, locus tables, and validation tests remain on the server under this project.
- Software versions and R session information are under `provenance/`.
- A project-wide SHA256 manifest is generated after pipeline completion.
"""
    (report_dir / "EXECUTION_REPORT.md").write_text(report)


if __name__ == "__main__":
    main()
