#!/usr/bin/env python3
"""Apply frozen cross-pair corrections and identify replicated local loci."""

import csv
import hashlib
import math
from pathlib import Path


ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
LOCAL = ROOT / "results" / "local"
OUT = LOCAL / "combined"
OUT.mkdir(parents=True, exist_ok=True)

PAIR_ORDER = ["HowardMDD__FinnGenIBD", "FinnGenDEP__deLangeIBD"]
PAIR_SOURCES = {
    "HowardMDD__FinnGenIBD": LOCAL / "HowardMDD__FinnGenIBD",
    "FinnGenDEP__deLangeIBD": LOCAL / "FinnGenDEP__deLangeIBD",
}


def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(p_values):
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [math.nan] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_from_end, (original_index, p_value) in enumerate(
        reversed(indexed), start=1
    ):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = running
    return adjusted


univ_rows = []
bivar_rows = []
for pair_id in PAIR_ORDER:
    source = PAIR_SOURCES[pair_id]
    if not (source / "DONE").exists():
        raise RuntimeError(f"Missing DONE marker for {pair_id}")
    univ_rows.extend(read_tsv(source / "univ_all.tsv"))
    bivar_rows.extend(read_tsv(source / "bivar_eligible.tsv"))

if bivar_rows:
    q_values = bh_adjust([float(row["p"]) for row in bivar_rows])
    for row, q_value in zip(bivar_rows, q_values):
        row["q_global_primary"] = f"{q_value:.12g}"
        row["passes_q05"] = str(q_value < 0.05).lower()

univ_fields = list(univ_rows[0]) if univ_rows else []
bivar_fields = list(dict.fromkeys(bivar_rows[0].keys())) if bivar_rows else []

write_tsv(OUT / "univ_all_primary.tsv", univ_rows, univ_fields)
write_tsv(OUT / "bivar_all_primary.tsv", bivar_rows, bivar_fields)

significant = [
    row for row in bivar_rows if float(row["q_global_primary"]) < 0.05
]
write_tsv(OUT / "bivar_q05_primary.tsv", significant, bivar_fields)

by_pair_locus = {
    (row["pair_id"], row["locus"]): row
    for row in bivar_rows
}
all_loci = sorted(
    {row["locus"] for row in bivar_rows},
    key=lambda value: int(float(value)),
)
replication_rows = []
for locus in all_loci:
    first = by_pair_locus.get((PAIR_ORDER[0], locus))
    second = by_pair_locus.get((PAIR_ORDER[1], locus))
    if first is None or second is None:
        status = "eligible_in_one_pair_only"
    else:
        first_sig = float(first["q_global_primary"]) < 0.05
        second_sig = float(second["q_global_primary"]) < 0.05
        same_direction = float(first["rho"]) * float(second["rho"]) > 0
        if first_sig and second_sig and same_direction:
            status = "replicated_same_direction"
        elif first_sig and second_sig:
            status = "significant_both_direction_discordant"
        else:
            status = "not_replicated_at_q05"

    template = first or second
    replication_rows.append(
        {
            "locus": locus,
            "chr": template["chr"],
            "start": template["start"],
            "stop": template["stop"],
            "pair1_rho": first["rho"] if first else "NA",
            "pair1_rho_lower": first["rho.lower"] if first else "NA",
            "pair1_rho_upper": first["rho.upper"] if first else "NA",
            "pair1_p": first["p"] if first else "NA",
            "pair1_q": first["q_global_primary"] if first else "NA",
            "pair2_rho": second["rho"] if second else "NA",
            "pair2_rho_lower": second["rho.lower"] if second else "NA",
            "pair2_rho_upper": second["rho.upper"] if second else "NA",
            "pair2_p": second["p"] if second else "NA",
            "pair2_q": second["q_global_primary"] if second else "NA",
            "replication_status": status,
        }
    )

replication_fields = [
    "locus", "chr", "start", "stop",
    "pair1_rho", "pair1_rho_lower", "pair1_rho_upper", "pair1_p", "pair1_q",
    "pair2_rho", "pair2_rho_lower", "pair2_rho_upper", "pair2_p", "pair2_q",
    "replication_status",
]
write_tsv(OUT / "local_reproducibility.tsv", replication_rows, replication_fields)

replicated = [
    row for row in replication_rows
    if row["replication_status"] == "replicated_same_direction"
]
write_tsv(OUT / "replicated_local_components.tsv", replicated, replication_fields)

eligible_in_both = [
    row for row in replication_rows
    if row["pair1_rho"] != "NA" and row["pair2_rho"] != "NA"
]
write_tsv(
    OUT / "same_block_eligible_both.tsv",
    eligible_in_both,
    replication_fields,
)

summary_rows = []
for pair_id in PAIR_ORDER:
    pair_univ = [row for row in univ_rows if row["pair_id"] == pair_id]
    pair_bivar = [row for row in bivar_rows if row["pair_id"] == pair_id]
    phenotypes = sorted({row["phen"] for row in pair_univ})
    counts = {
        pheno: sum(
            float(row["p"]) < float(row["univ_threshold"])
            for row in pair_univ if row["phen"] == pheno
        )
        for pheno in phenotypes
    }
    summary_rows.append(
        {
            "pair_id": pair_id,
            "univ_qualified_by_trait": ";".join(
                f"{pheno}:{counts[pheno]}" for pheno in phenotypes
            ),
            "n_bivar_eligible": len(pair_bivar),
            "n_bivar_q05": sum(
                float(row["q_global_primary"]) < 0.05 for row in pair_bivar
            ),
        }
    )

write_tsv(
    OUT / "local_pair_summary.tsv",
    summary_rows,
    ["pair_id", "univ_qualified_by_trait", "n_bivar_eligible", "n_bivar_q05"],
)

decision = (
    "layer2_pass_replicated_component"
    if replicated else
    "layer2_stop_no_replicated_component"
)
write_tsv(
    OUT / "layer2_decision.tsv",
    [{
        "decision": decision,
        "n_replicated_components": len(replicated),
        "n_same_block_eligible_both": len(eligible_in_both),
        "fine_mapping_permitted": str(bool(replicated)).lower(),
    }],
    [
        "decision", "n_replicated_components",
        "n_same_block_eligible_both", "fine_mapping_permitted",
    ],
)

with (OUT / "result_sha256.txt").open("w") as handle:
    for path in sorted(OUT.glob("*.tsv")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        handle.write(f"{digest}  {path}\n")

print(decision)
print(f"eligible_bivariate_tests={len(bivar_rows)}")
print(f"replicated_components={len(replicated)}")

old_summary = read_tsv(
    ROOT / "results" / "local" / "combined" / "local_pair_summary.tsv"
)
old_by_pair = {row["pair_id"]: row for row in old_summary}
new_by_pair = {row["pair_id"]: row for row in summary_rows}
old_same = read_tsv(
    ROOT / "results" / "local" / "combined" / "same_block_eligible_both.tsv"
)
old_replicated = read_tsv(
    ROOT / "results" / "local" / "combined" / "replicated_local_components.tsv"
)
report_lines = [
    "# LAVA pre/post Howard sample-size correction",
    "",
]
for pair_id in PAIR_ORDER:
    report_lines.extend(
        [
            f"## {pair_id}",
            "",
            f"- Eligible bivariate tests: {old_by_pair[pair_id]['n_bivar_eligible']} before; {new_by_pair[pair_id]['n_bivar_eligible']} after.",
            f"- Pooled-BH q < 0.05 tests: {old_by_pair[pair_id]['n_bivar_q05']} before; {new_by_pair[pair_id]['n_bivar_q05']} after.",
            "",
        ]
    )
report_lines.extend(
    [
        f"- Blocks estimable in both pairings: {len(old_same)} before; {len(eligible_in_both)} after.",
        f"- Replicated same-direction components: {len(old_replicated)} before; {len(replicated)} after.",
        f"- Corrected layer-2 decision: {decision}.",
        "- The FinnGen depression x de Lange IBD raw LAVA estimates were not rerun; their q values were recalculated jointly within the updated combined family.",
        "",
    ]
)
(CORR / "reports" / "LAVA_PREPOST_COMPARISON.md").write_text(
    "\n".join(report_lines)
)
