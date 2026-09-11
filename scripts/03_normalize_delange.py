#!/usr/bin/env python3
import csv
import gzip
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
RAW = ROOT / "data" / "raw"
PREP = ROOT / "data" / "prepared"
PREP.mkdir(parents=True, exist_ok=True)

FILES = {
    "IBD_deLange2017": RAW / "IBD_deLange2017.h.tsv.gz",
    "CD_deLange2017": RAW / "CD_deLange2017.h.tsv.gz",
    "UC_deLange2017": RAW / "UC_deLange2017.h.tsv.gz",
}


def choose(row, *names):
    for name in names:
        value = row.get(name)
        if value not in (None, "", "NA", "NaN", "nan"):
            return value
    return ""


for label, source in FILES.items():
    target = PREP / f"{label}.normalized.tsv.gz"
    kept = 0
    with gzip.open(source, "rt", errors="replace") as src, gzip.open(
        target, "wt", compresslevel=6
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(
            dst, fieldnames=["SNP", "A1", "A2", "BETA", "P"], delimiter="\t"
        )
        writer.writeheader()
        for row in reader:
            out = {
                "SNP": choose(row, "hm_rsid", "variant_id"),
                "A1": choose(row, "hm_effect_allele", "effect_allele"),
                "A2": choose(row, "hm_other_allele", "other_allele"),
                "BETA": choose(row, "hm_beta", "beta"),
                "P": choose(row, "p_value", "hm_p_value"),
            }
            if all(out.values()) and out["SNP"].startswith("rs"):
                writer.writerow(out)
                kept += 1
    print(label, kept, target)
