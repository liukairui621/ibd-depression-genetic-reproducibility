#!/usr/bin/env python3
import csv
import math
import re
from pathlib import Path


ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
GLOBAL = ROOT / "results" / "global_extension"
PAIR_FILE = ROOT / "provenance" / "extension_pair_manifest.tsv"


def grab(pattern, text, flags=0):
    match = re.search(pattern, text, flags)
    return float(match.group(1)) if match else math.nan


def grab_all(pattern, text, flags=0):
    return [float(value) for value in re.findall(pattern, text, flags)]


def final_snp_count(text):
    valid = re.findall(r"(\d+) SNPs with valid alleles", text)
    if valid:
        return int(valid[-1])
    merged = re.findall(
        r"After merging with (?:regression SNP LD|summary statistics), "
        r"(\d+) SNPs remain",
        text,
    )
    if merged:
        return int(merged[-1])
    raise RuntimeError("Could not parse the final LDSC SNP count")


def bh(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    adjusted = [math.nan] * n
    running = 1.0
    for reverse_rank in reversed(range(n)):
        index = order[reverse_rank]
        running = min(running, values[index] * n / (reverse_rank + 1))
        adjusted[index] = running
    return adjusted


def write_tsv(path, rows):
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


h2_rows = []
for log in sorted((GLOBAL / "h2").glob("*.log")):
    text = log.read_text(errors="replace")
    h2 = grab(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text)
    se = grab(
        r"Total Observed scale h2:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)",
        text,
    )
    h2_rows.append(
        {
            "trait": log.stem,
            "h2_obs": h2,
            "h2_se": se,
            "h2_z": h2 / se if math.isfinite(se) and se != 0 else math.nan,
            "intercept": grab(r"Intercept:\s*([-+0-9.eE]+)\s*\(", text),
            "intercept_se": grab(
                r"Intercept:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text
            ),
            "ratio": grab(r"Ratio:\s*([-+0-9.eE]+)\s*\(", text),
            "ratio_se": grab(
                r"Ratio:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text
            ),
            "mean_chisq": grab(r"Mean Chi\^2:\s*([-+0-9.eE]+)", text),
            "n_valid_snps": final_snp_count(text),
            "passes_h2_z4": int(
                math.isfinite(h2)
                and math.isfinite(se)
                and se > 0
                and h2 / se >= 4
            ),
        }
    )
write_tsv(GLOBAL / "h2_summary.tsv", h2_rows)

with PAIR_FILE.open(newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

rg_rows = []
for item in manifest:
    stem = f"{item['pair_class']}__{item['trait1']}__{item['trait2']}"
    log = GLOBAL / "rg" / f"{stem}.log"
    if not log.exists():
        raise FileNotFoundError(log)
    text = log.read_text(errors="replace")

    h2 = grab_all(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text)
    h2_se = grab_all(
        r"Total Observed scale h2:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)",
        text,
    )
    intercept = grab_all(r"^Intercept:\s*([-+0-9.eE]+)\s*\(", text, re.M)
    intercept_se = grab_all(
        r"^Intercept:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text, re.M
    )
    rg = grab(r"Genetic Correlation:\s*([-+0-9.eE]+)\s*\(", text)
    se = grab(
        r"Genetic Correlation:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)",
        text,
    )
    gcov = grab(r"Total Observed scale gencov:\s*([-+0-9.eE]+)\s*\(", text)
    gcov_se = grab(
        r"Total Observed scale gencov:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)",
        text,
    )
    row = dict(item)
    row.update(
        {
            "n_valid_snps": final_snp_count(text),
            "trait1_h2_pair": h2[0] if len(h2) > 0 else math.nan,
            "trait1_h2_se_pair": h2_se[0] if len(h2_se) > 0 else math.nan,
            "trait1_intercept_pair": (
                intercept[0] if len(intercept) > 0 else math.nan
            ),
            "trait1_intercept_se_pair": (
                intercept_se[0] if len(intercept_se) > 0 else math.nan
            ),
            "trait2_h2_pair": h2[1] if len(h2) > 1 else math.nan,
            "trait2_h2_se_pair": h2_se[1] if len(h2_se) > 1 else math.nan,
            "trait2_intercept_pair": (
                intercept[1] if len(intercept) > 1 else math.nan
            ),
            "trait2_intercept_se_pair": (
                intercept_se[1] if len(intercept_se) > 1 else math.nan
            ),
            "gcov": gcov,
            "gcov_se": gcov_se,
            "gcov_z": (
                gcov / gcov_se
                if math.isfinite(gcov_se) and gcov_se != 0
                else math.nan
            ),
            "gcov_intercept": (
                intercept[2] if len(intercept) > 2 else math.nan
            ),
            "gcov_intercept_se": (
                intercept_se[2] if len(intercept_se) > 2 else math.nan
            ),
            "gcov_intercept_z": (
                intercept[2] / intercept_se[2]
                if len(intercept) > 2
                and len(intercept_se) > 2
                and intercept_se[2] != 0
                else math.nan
            ),
            "rg": rg,
            "se": se,
            "z": (
                rg / se
                if math.isfinite(rg) and math.isfinite(se) and se != 0
                else math.nan
            ),
            "p": grab(
                r"Genetic Correlation.*?\nP:\s*([-+0-9.eE]+)",
                text,
                re.S,
            ),
        }
    )
    row["gcov_intercept_flag"] = int(
        math.isfinite(row["gcov_intercept_z"])
        and abs(row["gcov_intercept_z"]) >= 2
    )
    row["both_pair_h2_z4"] = int(
        len(h2) >= 2
        and len(h2_se) >= 2
        and h2_se[0] > 0
        and h2_se[1] > 0
        and h2[0] / h2_se[0] >= 4
        and h2[1] / h2_se[1] >= 4
    )
    rg_rows.append(row)

if len(h2_rows) != 12:
    raise RuntimeError(f"Expected 12 h2 rows, observed {len(h2_rows)}")
for row in h2_rows:
    for field in ("h2_obs", "h2_se", "intercept", "mean_chisq"):
        if not math.isfinite(row[field]):
            raise RuntimeError(
                f"Non-finite {field} in h2 log for {row['trait']}"
            )

if len(rg_rows) != 51:
    raise RuntimeError(f"Expected 51 rg rows, observed {len(rg_rows)}")
for row in rg_rows:
    for field in (
        "trait1_h2_pair",
        "trait2_h2_pair",
        "gcov",
        "gcov_intercept",
        "rg",
        "se",
        "p",
    ):
        if not math.isfinite(row[field]):
            raise RuntimeError(
                f"Non-finite {field} in rg log for "
                f"{row['trait1']} x {row['trait2']}"
            )

for row in rg_rows:
    row["fdr_all_cross_disease"] = math.nan
    row["fdr_primary_extension"] = math.nan
    row["fdr_within_pair_class"] = math.nan

cross_idx = [
    index
    for index, row in enumerate(rg_rows)
    if row["pair_class"] == "cross_disease" and math.isfinite(row["p"])
]
for index, q_value in zip(
    cross_idx, bh([rg_rows[index]["p"] for index in cross_idx])
):
    rg_rows[index]["fdr_all_cross_disease"] = q_value

primary_idx = [
    index
    for index, row in enumerate(rg_rows)
    if row["pair_class"] == "cross_disease"
    and row["primary_extension"] == "1"
    and math.isfinite(row["p"])
]
for index, q_value in zip(
    primary_idx, bh([rg_rows[index]["p"] for index in primary_idx])
):
    rg_rows[index]["fdr_primary_extension"] = q_value

classes = sorted({row["pair_class"] for row in rg_rows})
for pair_class in classes:
    indices = [
        index
        for index, row in enumerate(rg_rows)
        if row["pair_class"] == pair_class and math.isfinite(row["p"])
    ]
    for index, q_value in zip(
        indices, bh([rg_rows[index]["p"] for index in indices])
    ):
        rg_rows[index]["fdr_within_pair_class"] = q_value

write_tsv(GLOBAL / "rg_summary.tsv", rg_rows)
write_tsv(
    GLOBAL / "rg_cross_disease.tsv",
    [row for row in rg_rows if row["pair_class"] == "cross_disease"],
)
write_tsv(
    GLOBAL / "rg_ukb_definition_matrix.tsv",
    [
        row
        for row in rg_rows
        if row["pair_class"] == "cross_disease"
        and row["primary_extension"] == "1"
    ],
)

print(f"h2 rows: {len(h2_rows)}")
print(f"rg rows: {len(rg_rows)}")
print(f"cross-disease rows: {len(cross_idx)}")
print(f"primary extension rows: {len(primary_idx)}")
