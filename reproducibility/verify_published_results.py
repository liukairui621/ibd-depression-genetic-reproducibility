#!/usr/bin/env python3
"""Recompute BH correction and principal counts independently from exports."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    with (ROOT / relative).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh(rows, pkey="p"):
    order = sorted(range(len(rows)), key=lambda i: float(rows[i][pkey]))
    values = [1.0] * len(rows)
    previous = 1.0
    for rank in range(len(rows), 0, -1):
        index = order[rank - 1]
        previous = min(previous, float(rows[index][pkey]) * len(rows) / rank, 1.0)
        values[index] = previous
    return values


def main():
    checks = []
    families = [
        ("local", "results/derived/local/bivar_all_primary.tsv", "p", "q_global_primary", 110, 4),
        ("extension", "results/derived/global_extension/rg_ukb_definition_matrix.tsv", "p", "fdr_primary_extension", 24, 15),
        ("annotation_covariance", "extensions/annotation_stratified_gcov/results/PRIMARY_ANNOTATION_DECISION.tsv", "covariance_p_two_sided_discovery", "discovery_fdr_bh", 10, 2),
    ]
    for label, source, pkey, qkey, n, nsig in families:
        rows = read(source)
        q = bh(rows, pkey)
        delta = max(abs(float(row[qkey]) - computed) for row, computed in zip(rows, q))
        result = {"test": label, "n": len(rows), "significant": sum(x < .05 for x in q), "max_q_error": delta}
        assert len(rows) == n and result["significant"] == nsig and delta < 1e-10, result
        checks.append(result)
    extension = read(families[1][1])
    retained = [row for row in extension if row["gcov_intercept_flag"] == "0"]
    assert len(retained) == 18 and all(float(x["rg"]) > 0 for x in extension)
    assert sum(q < .05 for q in bh(retained)) == 9
    checks.append({"test": "unflagged_extension", "positive": len(retained), "significant": 9})
    same = read("results/derived/local/same_block_eligible_both.tsv")
    assert len(same) == 3
    reproduced = [x for x in same if float(x["pair1_q"]) < .05 and float(x["pair2_q"]) < .05 and float(x["pair1_rho"]) * float(x["pair2_rho"]) > 0]
    assert not reproduced
    annotations = read("extensions/annotation_stratified_gcov/results/discovery_annotation_covariance_all.tsv")
    assert len(annotations) == 84
    checks.append({"test": "estimability_and_model_size", "jointly_estimable_local_blocks": 3, "local_replications": 0, "fitted_annotations": 84})
    print(json.dumps({"status": "numerical_checks_passed", "scope": "published exports; no source GWAS rerun", "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
