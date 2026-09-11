#!/usr/bin/env python3
import gzip
from pathlib import Path

path = Path("/root/IBD/20_Reproducibility_Ladder/data/munged_core/MDD_Howard2019.sumstats.gz")
rows = 0
nonmissing = 0
n_values = set()
with gzip.open(path, "rt") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    n_idx = header.index("N")
    z_idx = header.index("Z")
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if not fields or len(fields) != len(header):
            continue
        rows += 1
        if fields[z_idx] not in {"", "NA", "nan", "NaN"}:
            nonmissing += 1
        if fields[n_idx] not in {"", "NA", "nan", "NaN"}:
            n_values.add(float(fields[n_idx]))

if n_values != {500199.0}:
    raise SystemExit(f"Unexpected nonmissing N values: {sorted(n_values)[:10]}")
if nonmissing != 6465253:
    raise SystemExit(f"Expected 6465253 nonmissing statistics, found {nonmissing}")
print(f"Template rows: {rows}")
print(f"Nonmissing statistics: {nonmissing}")
print("Verified constant nonmissing N: 500199")
