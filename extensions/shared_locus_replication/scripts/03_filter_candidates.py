#!/usr/bin/env python3
import csv
import gzip
import os
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")


def main():
    cutoffs = {}
    with (ROOT / "results/raw/placo_product_cutoffs.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cutoffs[row["pair_id"]] = (
                float(row["positive_product_cutoff"]),
                float(row["negative_product_abs_cutoff"]),
            )

    summary = []
    for pair_id, (pos_cut, neg_cut) in cutoffs.items():
        source = ROOT / "results/raw" / f"{pair_id}.harmonized.tsv.gz"
        target = ROOT / "results/raw" / f"{pair_id}.integration_candidates.tsv.gz"
        temp = target.with_suffix(target.suffix + ".tmp")
        total = selected = 0
        with gzip.open(source, "rt") as inp, gzip.open(temp, "wt", compresslevel=6) as out:
            header = inp.readline()
            columns = header.rstrip("\n\r").split("\t")
            product_idx = columns.index("Z_PRODUCT")
            out.write(header)
            for line in inp:
                total += 1
                fields = line.rstrip("\n\r").split("\t")
                product = float(fields[product_idx])
                if product >= pos_cut * (1 - 1e-12) or product <= -neg_cut * (1 - 1e-12):
                    out.write(line)
                    selected += 1
        os.replace(temp, target)
        summary.append(
            {
                "pair_id": pair_id,
                "harmonized_variants": total,
                "integration_candidates": selected,
                "positive_product_cutoff": pos_cut,
                "negative_product_abs_cutoff": neg_cut,
            }
        )
        print(f"{pair_id}: selected {selected:,} of {total:,}")

    with (ROOT / "results/raw/candidate_filter_stats.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
