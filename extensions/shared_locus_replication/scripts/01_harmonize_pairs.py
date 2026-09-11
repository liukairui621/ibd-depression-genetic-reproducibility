#!/usr/bin/env python3
import csv
import gzip
import math
import os
from itertools import zip_longest
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")
BASE = Path("/root/IBD/20_Reproducibility_Ladder")

INPUTS = {
    "Howard": BASE / "corrections/20260806_howard_no23andme/work/munged/MDD_Howard2019_no23andMe.sumstats.gz",
    "deLangeIBD": BASE / "data/munged_core/IBD_deLange2017.sumstats.gz",
    "FinnGenDEP": BASE / "data/munged_core/DEP_FinnGen_R12.sumstats.gz",
    "FinnGenIBD": BASE / "data/munged_core/IBD_FinnGen_R12.sumstats.gz",
}

PAIRS = [
    ("External_Howard_deLange", "Howard", "deLangeIBD", "placo"),
    ("FinnGen_DEP_IBD", "FinnGenDEP", "FinnGenIBD", "placo_plus"),
    ("Diagonal_Howard_FinnGenIBD", "Howard", "FinnGenIBD", "placo"),
    ("Diagonal_FinnGenDEP_deLange", "FinnGenDEP", "deLangeIBD", "placo"),
]


class OnlineBivariate:
    def __init__(self):
        self.n = 0
        self.mx = 0.0
        self.my = 0.0
        self.m2x = 0.0
        self.m2y = 0.0
        self.cxy = 0.0

    def update(self, x, y):
        self.n += 1
        dx = x - self.mx
        self.mx += dx / self.n
        dy = y - self.my
        self.my += dy / self.n
        self.m2x += dx * (x - self.mx)
        self.m2y += dy * (y - self.my)
        self.cxy += dx * (y - self.my)

    def var_x(self):
        return self.m2x / (self.n - 1)

    def var_y(self):
        return self.m2y / (self.n - 1)

    def cor(self):
        return self.cxy / math.sqrt(self.m2x * self.m2y)


def two_sided_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def parse_header(handle):
    header = handle.readline().rstrip("\n\r").split("\t")
    required = {"SNP", "A1", "A2", "Z"}
    if not required.issubset(header):
        raise ValueError(f"Unexpected LDSC schema: {header}")
    # LDSC Z is defined with respect to the semantic A1 column. Some FinnGen
    # outputs are ordered SNP,A2,A1,Z,N, so positional allele parsing is unsafe.
    return header, header.index("SNP"), header.index("A1"), header.index("A2"), header.index("Z")


def harmonize_pair(pair_id, dep_path, ibd_path, method):
    out_path = ROOT / "results/raw" / f"{pair_id}.harmonized.tsv.gz"
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    variance_null = OnlineBivariate()
    correlation_null = OnlineBivariate()
    counts = {
        "rows": 0,
        "numeric_both": 0,
        "allele_match": 0,
        "allele_flip": 0,
        "allele_mismatch": 0,
    }

    with gzip.open(dep_path, "rt") as dep, gzip.open(ibd_path, "rt") as ibd, gzip.open(
        tmp_path, "wt", compresslevel=6
    ) as out:
        _, ds, dea, doa, dz = parse_header(dep)
        _, is_, iea, ioa, iz = parse_header(ibd)
        out.write("SNP\tEA\tOA\tZ_DEP\tZ_IBD\tP_DEP\tP_IBD\tZ_PRODUCT\tPRODUCT_SIGN\n")

        for line_number, (dline, iline) in enumerate(zip_longest(dep, ibd), start=2):
            if dline is None or iline is None:
                raise RuntimeError(f"Input length mismatch at line {line_number} for {pair_id}")
            counts["rows"] += 1
            d = dline.rstrip("\n\r").split("\t")
            i = iline.rstrip("\n\r").split("\t")
            if d[ds] != i[is_]:
                raise RuntimeError(
                    f"SNP order mismatch for {pair_id} at line {line_number}: {d[ds]} vs {i[is_]}"
                )
            try:
                z_dep = float(d[dz])
                z_ibd = float(i[iz])
            except (ValueError, IndexError):
                continue
            if not (math.isfinite(z_dep) and math.isfinite(z_ibd)):
                continue
            counts["numeric_both"] += 1

            dep_ea, dep_oa = d[dea].upper(), d[doa].upper()
            ibd_ea, ibd_oa = i[iea].upper(), i[ioa].upper()
            if dep_ea == ibd_ea and dep_oa == ibd_oa:
                counts["allele_match"] += 1
            elif dep_ea == ibd_oa and dep_oa == ibd_ea:
                z_ibd = -z_ibd
                counts["allele_flip"] += 1
            else:
                counts["allele_mismatch"] += 1
                continue

            p_dep = two_sided_p(z_dep)
            p_ibd = two_sided_p(z_ibd)
            product = z_dep * z_ibd
            sign = "positive" if product > 0 else "negative" if product < 0 else "zero"

            # Official var.placo removes rows only when both traits pass the threshold.
            if not (p_dep < 1e-4 and p_ibd < 1e-4):
                variance_null.update(z_dep, z_ibd)
            # Official cor.pearson removes rows when either trait passes the threshold.
            if not (p_dep < 1e-4 or p_ibd < 1e-4):
                correlation_null.update(z_dep, z_ibd)

            out.write(
                f"{d[ds]}\t{dep_ea}\t{dep_oa}\t{z_dep:.8g}\t{z_ibd:.8g}\t"
                f"{p_dep:.12g}\t{p_ibd:.12g}\t{product:.12g}\t{sign}\n"
            )

    os.replace(tmp_path, out_path)
    return {
        "pair_id": pair_id,
        "method": method,
        "depression_input": str(dep_path),
        "ibd_input": str(ibd_path),
        "harmonized_output": str(out_path),
        **counts,
        "retained": counts["allele_match"] + counts["allele_flip"],
        "variance_n": variance_null.n,
        "var_z_dep": variance_null.var_x(),
        "var_z_ibd": variance_null.var_y(),
        "correlation_n": correlation_null.n,
        "cor_z": correlation_null.cor(),
    }


def main():
    (ROOT / "results/raw").mkdir(parents=True, exist_ok=True)
    rows = []
    for pair_id, dep_name, ibd_name, method in PAIRS:
        print(f"Harmonizing {pair_id}", flush=True)
        rows.append(harmonize_pair(pair_id, INPUTS[dep_name], INPUTS[ibd_name], method))

    out = ROOT / "results/raw/pair_harmonization_stats.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
