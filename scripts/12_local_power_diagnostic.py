#!/usr/bin/env python3
"""Approximate normal-theory detectable-rho diagnostic for eligible LAVA tests."""

import csv
import math
from pathlib import Path
from statistics import NormalDist, median

BASE = Path(__file__).resolve().parent.parent
LOCAL = BASE / "results" / "derived" / "local"
OUT = BASE / "results" / "derived" / "local_power"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quantile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


rows = read_tsv(LOCAL / "bivar_all_primary.tsv")
if not rows:
    raise SystemExit("No combined eligible LAVA tests found")

n_tests = len(rows)
z_power_80 = NormalDist().inv_cdf(0.80)
z_nominal = NormalDist().inv_cdf(1 - 0.05 / 2)
family_alpha = 0.05 / n_tests
z_family = NormalDist().inv_cdf(1 - family_alpha / 2)

output_rows = []
for row in rows:
    lower = float(row["rho.lower"])
    upper = float(row["rho.upper"])
    se_from_ci = (upper - lower) / (2 * 1.96)
    output_rows.append(
        {
            "pair_id": row["pair_id"],
            "locus": row["locus"],
            "chr": row["chr"],
            "start": row["start"],
            "stop": row["stop"],
            "rho": row["rho"],
            "rho_lower": lower,
            "rho_upper": upper,
            "se_approx_from_95ci": se_from_ci,
            "mde_abs_rho_80pct_nominal_alpha_0_05": (z_nominal + z_power_80) * se_from_ci,
            "mde_abs_rho_80pct_bonferroni_family": (z_family + z_power_80) * se_from_ci,
            "bonferroni_family_alpha": family_alpha,
            "n_eligible_tests": n_tests,
        }
    )

fields = list(output_rows[0])
with (OUT / "approximate_detectable_rho.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(output_rows)

summary_rows = []
for pair_id in sorted({row["pair_id"] for row in output_rows}):
    values = [
        row["mde_abs_rho_80pct_bonferroni_family"]
        for row in output_rows
        if row["pair_id"] == pair_id
    ]
    summary_rows.append(
        {
            "pair_id": pair_id,
            "n_tests": len(values),
            "median_mde80_family": median(values),
            "q1_mde80_family": quantile(values, 0.25),
            "q3_mde80_family": quantile(values, 0.75),
            "min_mde80_family": min(values),
            "max_mde80_family": max(values),
            "n_mde80_gt_1": sum(value > 1 for value in values),
        }
    )

with (OUT / "approximate_detectable_rho_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_rows[0].keys(),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

report_lines = [
    "# Approximate local power diagnostic",
    "",
    f"Eligible bivariate tests in the jointly corrected family: {n_tests}.",
    f"Conservative Bonferroni alpha used for the familywise diagnostic: {family_alpha:.6g}.",
    "",
    "The standard error was approximated from each LAVA 95% confidence interval as (upper-lower)/(2x1.96). The minimum detectable absolute rho for 80% power was then calculated with a two-sided normal approximation. This is a descriptive power diagnostic, not a formal LAVA simulation and not the threshold used for the primary BH-FDR decision.",
    "",
]
for row in summary_rows:
    report_lines.append(
        f"- {row['pair_id']}: median familywise MDE80={row['median_mde80_family']:.3f} "
        f"(IQR {row['q1_mde80_family']:.3f}-{row['q3_mde80_family']:.3f}); "
        f"{row['n_mde80_gt_1']}/{row['n_tests']} block estimates had MDE80 > 1."
    )
(BASE / "reports" / "LOCAL_POWER_DIAGNOSTIC.md").write_text(
    "\n".join(report_lines) + "\n"
)
print("\n".join(report_lines))
