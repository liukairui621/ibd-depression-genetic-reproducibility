#!/usr/bin/env python3
import csv
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sign(value):
    x = float(value)
    return 1 if x > 0 else -1 if x < 0 else 0


def main():
    replication = [x for x in read(ROOT / "results/loci/cross_cohort_replication.tsv") if x["PRIMARY_REPLICATED"] == "TRUE"]
    support = read(ROOT / "results/loci/four_pair_lead_snp_support.tsv")
    lookup = {(x["LEAD_SNP"], x["PAIR_ID"]): x for x in support if x["FOUND"] == "TRUE"}
    out = []
    for row in replication:
        discovery = lookup[(row["LEAD_SNP"], row["discovery_pair"])]
        validation = lookup[(row["LEAD_SNP"], row["validation_pair"])]
        val_dep = float(validation["Z_DEP"])
        val_ibd = float(validation["Z_IBD"])
        if discovery["EA"] == validation["EA"] and discovery["OA"] == validation["OA"]:
            alignment = "same"
        elif discovery["EA"] == validation["OA"] and discovery["OA"] == validation["EA"]:
            alignment = "swapped_validation_flipped"
            val_dep = -val_dep
            val_ibd = -val_ibd
        else:
            alignment = "mismatch"
        dep_ok = alignment != "mismatch" and sign(discovery["Z_DEP"]) == sign(val_dep)
        ibd_ok = alignment != "mismatch" and sign(discovery["Z_IBD"]) == sign(val_ibd)
        out.append({
            "LOC": row["LOC"], "LEAD_SNP": row["LEAD_SNP"], "direction": row["direction"],
            "discovery_pair": row["discovery_pair"], "validation_pair": row["validation_pair"],
            "DISCOVERY_EA": discovery["EA"], "DISCOVERY_OA": discovery["OA"],
            "VALIDATION_EA_ORIGINAL": validation["EA"], "VALIDATION_OA_ORIGINAL": validation["OA"],
            "ALLELE_ALIGNMENT": alignment, "DISCOVERY_Z_DEP": discovery["Z_DEP"],
            "VALIDATION_Z_DEP_ALIGNED": f"{val_dep:.8g}", "DEP_DIRECTION_CONCORDANT": str(dep_ok).upper(),
            "DISCOVERY_Z_IBD": discovery["Z_IBD"], "VALIDATION_Z_IBD_ALIGNED": f"{val_ibd:.8g}",
            "IBD_DIRECTION_CONCORDANT": str(ibd_ok).upper(),
            "BOTH_TRAIT_DIRECTIONS_CONCORDANT": str(dep_ok and ibd_ok).upper(),
        })
    fields = list(out[0]) if out else ["LOC", "LEAD_SNP", "direction", "BOTH_TRAIT_DIRECTIONS_CONCORDANT"]
    with (ROOT / "reports/PRIMARY_REPLICATED_DIRECTION_AUDIT.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main()
