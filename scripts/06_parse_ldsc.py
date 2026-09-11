#!/usr/bin/env python3
import csv
import math
import re
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
GLOBAL = ROOT / "results" / "global"


def grab(pattern, text, flags=0):
    match = re.search(pattern, text, flags)
    return float(match.group(1)) if match else math.nan


def bh(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    adjusted = [math.nan] * len(values)
    running = 1.0
    n = len(values)
    for rank0 in reversed(range(n)):
        i = order[rank0]
        running = min(running, values[i] * n / (rank0 + 1))
        adjusted[i] = running
    return adjusted


h2_rows = []
for log in sorted((GLOBAL / "h2").glob("*.log")):
    text = log.read_text(errors="replace")
    h2 = grab(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text)
    se = grab(r"Total Observed scale h2:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text)
    h2_rows.append(
        {
            "trait": log.stem,
            "h2_obs": h2,
            "h2_se": se,
            "h2_z": h2 / se if se and not math.isnan(se) else math.nan,
            "intercept": grab(r"Intercept:\s*([-+0-9.eE]+)\s*\(", text),
            "intercept_se": grab(r"Intercept:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text),
            "ratio": grab(r"Ratio:\s*([-+0-9.eE]+)\s*\(", text),
            "mean_chisq": grab(r"Mean Chi\^2:\s*([-+0-9.eE]+)", text),
        }
    )

with (GLOBAL / "h2_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=h2_rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(h2_rows)

pair_meta = {}
with (ROOT / "provenance" / "layer1_pair_manifest.tsv").open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        pair_meta[(row["trait1"], row["trait2"])] = row["pair_class"]

rg_rows = []
for log in sorted((GLOBAL / "rg").glob("*.log")):
    text = log.read_text(errors="replace")
    parts = log.stem.split("__")
    if len(parts) != 3:
        continue
    cls, t1, t2 = parts
    rg = grab(r"Genetic Correlation:\s*([-+0-9.eE]+)\s*\(", text)
    se = grab(r"Genetic Correlation:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text)
    gcov = grab(r"Genetic Covariance.*?Intercept:\s*([-+0-9.eE]+)\s*\(", text, re.S)
    gcov_se = grab(
        r"Genetic Covariance.*?Intercept:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)",
        text,
        re.S,
    )
    rg_rows.append(
        {
            "pair_class": cls,
            "trait1": t1,
            "trait2": t2,
            "rg": rg,
            "se": se,
            "z": rg / se if se and not math.isnan(se) else math.nan,
            "p": grab(r"Genetic Correlation.*?\nP:\s*([-+0-9.eE]+)", text, re.S),
            "gcov_intercept": gcov,
            "gcov_intercept_se": gcov_se,
            "gcov_intercept_z": gcov / gcov_se
            if gcov_se and not math.isnan(gcov_se)
            else math.nan,
        }
    )

for cls in sorted(set(row["pair_class"] for row in rg_rows)):
    idx = [i for i, row in enumerate(rg_rows) if row["pair_class"] == cls]
    finite = [i for i in idx if not math.isnan(rg_rows[i]["p"])]
    q = bh([rg_rows[i]["p"] for i in finite])
    for i, value in zip(finite, q):
        rg_rows[i]["fdr_within_pair_class"] = value
for row in rg_rows:
    row.setdefault("fdr_within_pair_class", math.nan)

fields = list(rg_rows[0].keys())
with (GLOBAL / "rg_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rg_rows)

print("h2 rows", len(h2_rows))
print("rg rows", len(rg_rows))
