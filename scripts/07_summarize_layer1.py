#!/usr/bin/env python3
import csv
import math
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
GLOBAL = ROOT / "results" / "global"


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rows = read_tsv(GLOBAL / "rg_summary.tsv")
lookup = {(r["trait1"], r["trait2"]): r for r in rows}

anchors = [
    ("MDD_Howard2019", "IBD_FinnGen_R12"),
    ("DEP_FinnGen_R12", "IBD_deLange2017"),
]

anchor_rows = [lookup[x] for x in anchors]
positive = all(float(r["rg"]) > 0 for r in anchor_rows)
strict = all(float(r["p"]) < 0.025 for r in anchor_rows)
nominal = all(float(r["p"]) < 0.05 for r in anchor_rows)

if positive and strict:
    status = "global_gate_passed_reciprocal_estimates_reported_separately"
elif positive and nominal:
    status = "positive_but_primary_familywise_gate_not_met"
elif positive:
    status = "inconclusive_or_power_limited"
else:
    status = "cohort_or_definition_dependent"

out = GLOBAL / "layer1_primary_decision.tsv"
with out.open("w", newline="") as handle:
    fields = ["decision", "trait1", "trait2", "rg", "se", "p", "gcov_intercept_z"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in anchor_rows:
        writer.writerow(
            {
                "decision": status,
                "trait1": r["trait1"],
                "trait2": r["trait2"],
                "rg": r["rg"],
                "se": r["se"],
                "p": r["p"],
                "gcov_intercept_z": r["gcov_intercept_z"],
            }
        )

print(status)
for row in anchor_rows:
    print(row)
