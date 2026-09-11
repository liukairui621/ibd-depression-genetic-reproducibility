#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder")

INPUTS = {
    "MDD_Howard2019": Path("/root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz"),
    "DEP_FinnGen_R12": Path("/root/IBD/19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz"),
    "IBD_deLange2017_raw": Path("/root/IBD/00_RawData/GWAS/IBD/deLange2017/ibd_build37_59957_20161107.txt.gz"),
    "IBD_deLange2017_harmonised": Path("/root/IBD/00_RawData/GWAS/IBD/deLange2017/deLange2017_harmonised.tsv.gz"),
    "IBD_FinnGen_R12": Path("/root/IBD/00_RawData/GWAS/IBD/IBD_STRICT_FinnGen_R12.gz"),
    "CD_FinnGen_R12": Path("/root/IBD/00_RawData/GWAS/IBD/CD_STRICT2_FinnGen_R12.gz"),
    "UC_FinnGen_R12": Path("/root/IBD/00_RawData/GWAS/IBD/UC_STRICT2_FinnGen_R12.gz"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_header(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        return handle.readline().rstrip("\n\r")


def main() -> None:
    rows = []
    for label, path in INPUTS.items():
        exists = path.exists()
        size = path.stat().st_size if exists else 0
        rows.append(
            {
                "trait_or_input": label,
                "path": str(path),
                "exists": exists,
                "bytes": size,
                "header": read_header(path) if exists and size else "",
                "sha256": sha256(path) if exists and size else "",
            }
        )

    out = ROOT / "provenance" / "existing_input_inventory.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
