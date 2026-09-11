#!/usr/bin/env python3
import bisect
import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path("/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication")
LD_ROOT = Path("/root/IBD/00_RawData/Reference/LD_EUR/1000G_EUR_Phase3_plink")
BLOCK_FILE = Path("/root/IBD/00_RawData/Reference/LAVA/LAVA_blocks_fixed.txt")

PAIR_IDS = [
    "External_Howard_deLange",
    "FinnGen_DEP_IBD",
    "Diagonal_Howard_FinnGenIBD",
    "Diagonal_FinnGenDEP_deLange",
]
PRIMARY = ["External_Howard_deLange", "FinnGen_DEP_IBD"]


def read_tsv(path, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_blocks():
    by_chr = defaultdict(list)
    for row in read_tsv(BLOCK_FILE):
        block = {"LOC": int(row["LOC"]), "CHR": int(row["CHR"]), "START": int(row["START"]), "STOP": int(row["STOP"])}
        by_chr[block["CHR"]].append(block)
    starts = {}
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: x["START"])
        starts[chrom] = [x["START"] for x in by_chr[chrom]]
    return by_chr, starts


def map_positions(snps):
    positions = {}
    remaining = set(snps)
    for chrom in range(1, 23):
        if not remaining:
            break
        bim = LD_ROOT / f"1000G.EUR.QC.{chrom}.bim"
        with bim.open() as handle:
            for line in handle:
                fields = line.split()
                snp = fields[1]
                if snp in remaining:
                    positions[snp] = (int(fields[0]), int(fields[3]))
                    remaining.remove(snp)
    return positions, remaining


def assign_block(chrom, bp, blocks, starts):
    if chrom not in blocks:
        return None
    idx = bisect.bisect_right(starts[chrom], bp) - 1
    if idx < 0:
        return None
    block = blocks[chrom][idx]
    return block if bp <= block["STOP"] else None


def scan_harmonized(pair_id, targets):
    found = {}
    if not targets:
        return found
    path = ROOT / "results/raw" / f"{pair_id}.harmonized.tsv.gz"
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["SNP"] in targets:
                found[row["SNP"]] = row
                if len(found) == len(targets):
                    break
    return found


def main():
    stats = {row["pair_id"]: row for row in read_tsv(ROOT / "results/raw/pair_harmonization_stats.tsv")}
    selected_by_pair = {}
    needed_snps = set()
    for pair_id in PAIR_IDS:
        rows = []
        for row in read_tsv(ROOT / "results/raw" / f"{pair_id}.PLACO.tsv.gz", gz=True):
            if float(row["P_PLACO"]) < 1e-6:
                row["pair_id"] = pair_id
                rows.append(row)
                needed_snps.add(row["SNP"])
        selected_by_pair[pair_id] = rows

    positions, unmapped = map_positions(needed_snps)
    blocks, starts = load_blocks()
    loci_by_pair = {}
    all_variant_fields = [
        "pair_id", "SNP", "EA", "OA", "Z_DEP", "Z_IBD", "P_DEP", "P_IBD",
        "Z_PRODUCT", "PRODUCT_SIGN", "P_PLACO", "Q_BH_GENOMEWIDE", "METHOD",
        "CHR", "BP", "LOC", "BLOCK_START", "BLOCK_STOP", "MHC"
    ]
    locus_fields = [
        "pair_id", "LOC", "CHR", "BLOCK_START", "BLOCK_STOP", "MHC", "LEAD_SNP",
        "LEAD_BP", "LEAD_P_PLACO", "LEAD_Q_BH", "LEAD_Z_DEP", "LEAD_Z_IBD",
        "LEAD_Z_PRODUCT", "PRODUCT_SIGN", "N_VARIANTS_P_LT_1E6", "GWS_LOCUS"
    ]

    for pair_id, rows in selected_by_pair.items():
        annotated = []
        for row in rows:
            pos = positions.get(row["SNP"])
            if pos is None:
                row.update({"CHR": "", "BP": "", "LOC": "", "BLOCK_START": "", "BLOCK_STOP": "", "MHC": ""})
            else:
                chrom, bp = pos
                block = assign_block(chrom, bp, blocks, starts)
                row.update({
                    "CHR": chrom, "BP": bp,
                    "LOC": block["LOC"] if block else "",
                    "BLOCK_START": block["START"] if block else "",
                    "BLOCK_STOP": block["STOP"] if block else "",
                    "MHC": str(chrom == 6 and 25_000_000 <= bp <= 34_000_000).upper(),
                })
            annotated.append(row)
        write_tsv(ROOT / "results/loci" / f"{pair_id}.variants_p1e6.tsv", all_variant_fields, annotated)

        groups = defaultdict(list)
        for row in annotated:
            if row["LOC"] != "":
                groups[int(row["LOC"])].append(row)
        loci = []
        for loc, variants in groups.items():
            lead = min(variants, key=lambda x: float(x["P_PLACO"]))
            loci.append({
                "pair_id": pair_id, "LOC": loc, "CHR": lead["CHR"],
                "BLOCK_START": lead["BLOCK_START"], "BLOCK_STOP": lead["BLOCK_STOP"],
                "MHC": lead["MHC"], "LEAD_SNP": lead["SNP"], "LEAD_BP": lead["BP"],
                "LEAD_P_PLACO": lead["P_PLACO"], "LEAD_Q_BH": lead["Q_BH_GENOMEWIDE"],
                "LEAD_Z_DEP": lead["Z_DEP"], "LEAD_Z_IBD": lead["Z_IBD"],
                "LEAD_Z_PRODUCT": lead["Z_PRODUCT"], "PRODUCT_SIGN": lead["PRODUCT_SIGN"],
                "N_VARIANTS_P_LT_1E6": len(variants),
                "GWS_LOCUS": str(any(float(x["P_PLACO"]) < 5e-8 for x in variants)).upper(),
            })
        loci.sort(key=lambda x: float(x["LEAD_P_PLACO"]))
        loci_by_pair[pair_id] = loci
        write_tsv(ROOT / "results/loci" / f"{pair_id}.loci_p1e6.tsv", locus_fields, loci)

    discovery_loci = {pair: [x for x in loci_by_pair[pair] if x["GWS_LOCUS"] == "TRUE"] for pair in PRIMARY}
    combined_discovery = []
    for pair in PRIMARY:
        combined_discovery.extend(discovery_loci[pair])
    write_tsv(ROOT / "results/loci/primary_discovery_loci.tsv", locus_fields, combined_discovery)

    lead_targets = {x["LEAD_SNP"] for x in combined_discovery}
    harmonized_hits = {pair: scan_harmonized(pair, lead_targets) for pair in PAIR_IDS}

    replication_fields = [
        "direction", "discovery_pair", "validation_pair", "n_nonmhc_discovery_loci",
        "LOC", "CHR", "BLOCK_START", "BLOCK_STOP", "MHC", "LEAD_SNP", "LEAD_BP",
        "DISCOVERY_P_PLACO", "DISCOVERY_Q_BH", "DISCOVERY_Z_PRODUCT", "DISCOVERY_SIGN",
        "VALIDATION_FOUND", "VALIDATION_EA", "VALIDATION_OA", "VALIDATION_Z_DEP",
        "VALIDATION_Z_IBD", "VALIDATION_P_DEP", "VALIDATION_P_IBD", "VALIDATION_Z_PRODUCT",
        "VALIDATION_SIGN", "VALIDATION_METHOD", "VALIDATION_VAR_Z_DEP", "VALIDATION_VAR_Z_IBD",
        "VALIDATION_COR_Z"
    ]
    tests = []
    for discovery_pair, validation_pair, direction in [
        (PRIMARY[0], PRIMARY[1], "External_to_FinnGen"),
        (PRIMARY[1], PRIMARY[0], "FinnGen_to_External"),
    ]:
        n_nonmhc = sum(x["MHC"] != "TRUE" for x in discovery_loci[discovery_pair])
        for locus in discovery_loci[discovery_pair]:
            hit = harmonized_hits[validation_pair].get(locus["LEAD_SNP"])
            val_stats = stats[validation_pair]
            tests.append({
                "direction": direction, "discovery_pair": discovery_pair, "validation_pair": validation_pair,
                "n_nonmhc_discovery_loci": n_nonmhc, **{k: locus[k] for k in ["LOC", "CHR", "BLOCK_START", "BLOCK_STOP", "MHC", "LEAD_SNP", "LEAD_BP"]},
                "DISCOVERY_P_PLACO": locus["LEAD_P_PLACO"], "DISCOVERY_Q_BH": locus["LEAD_Q_BH"],
                "DISCOVERY_Z_PRODUCT": locus["LEAD_Z_PRODUCT"], "DISCOVERY_SIGN": locus["PRODUCT_SIGN"],
                "VALIDATION_FOUND": str(hit is not None).upper(),
                "VALIDATION_EA": hit["EA"] if hit else "", "VALIDATION_OA": hit["OA"] if hit else "",
                "VALIDATION_Z_DEP": hit["Z_DEP"] if hit else "", "VALIDATION_Z_IBD": hit["Z_IBD"] if hit else "",
                "VALIDATION_P_DEP": hit["P_DEP"] if hit else "", "VALIDATION_P_IBD": hit["P_IBD"] if hit else "",
                "VALIDATION_Z_PRODUCT": hit["Z_PRODUCT"] if hit else "", "VALIDATION_SIGN": hit["PRODUCT_SIGN"] if hit else "",
                "VALIDATION_METHOD": val_stats["method"], "VALIDATION_VAR_Z_DEP": val_stats["var_z_dep"],
                "VALIDATION_VAR_Z_IBD": val_stats["var_z_ibd"], "VALIDATION_COR_Z": val_stats["cor_z"],
            })
    write_tsv(ROOT / "results/loci/replication_tests_input.tsv", replication_fields, tests)

    support_fields = [
        "LEAD_SNP", "SOURCE_LOC", "SOURCE_PAIR", "PAIR_ID", "FOUND", "EA", "OA", "Z_DEP", "Z_IBD",
        "P_DEP", "P_IBD", "Z_PRODUCT", "PRODUCT_SIGN", "METHOD", "VAR_Z_DEP", "VAR_Z_IBD", "COR_Z"
    ]
    support = []
    source_by_snp = {x["LEAD_SNP"]: x for x in combined_discovery}
    for snp, source in source_by_snp.items():
        for pair in PAIR_IDS:
            hit = harmonized_hits[pair].get(snp)
            st = stats[pair]
            support.append({
                "LEAD_SNP": snp, "SOURCE_LOC": source["LOC"], "SOURCE_PAIR": source["pair_id"],
                "PAIR_ID": pair, "FOUND": str(hit is not None).upper(),
                "EA": hit["EA"] if hit else "", "OA": hit["OA"] if hit else "",
                "Z_DEP": hit["Z_DEP"] if hit else "", "Z_IBD": hit["Z_IBD"] if hit else "",
                "P_DEP": hit["P_DEP"] if hit else "", "P_IBD": hit["P_IBD"] if hit else "",
                "Z_PRODUCT": hit["Z_PRODUCT"] if hit else "", "PRODUCT_SIGN": hit["PRODUCT_SIGN"] if hit else "",
                "METHOD": st["method"], "VAR_Z_DEP": st["var_z_dep"], "VAR_Z_IBD": st["var_z_ibd"], "COR_Z": st["cor_z"],
            })
    write_tsv(ROOT / "results/loci/four_pair_support_input.tsv", support_fields, support)

    external_loci = {int(x["LOC"]): x for x in loci_by_pair[PRIMARY[0]]}
    finngen_loci = {int(x["LOC"]): x for x in loci_by_pair[PRIMARY[1]]}
    convergence = []
    for loc in sorted(set(external_loci) & set(finngen_loci)):
        a, b = external_loci[loc], finngen_loci[loc]
        convergence.append({
            "LOC": loc, "CHR": a["CHR"], "BLOCK_START": a["BLOCK_START"], "BLOCK_STOP": a["BLOCK_STOP"],
            "EXTERNAL_LEAD_SNP": a["LEAD_SNP"], "EXTERNAL_P": a["LEAD_P_PLACO"], "EXTERNAL_SIGN": a["PRODUCT_SIGN"],
            "FINNGEN_LEAD_SNP": b["LEAD_SNP"], "FINNGEN_P": b["LEAD_P_PLACO"], "FINNGEN_SIGN": b["PRODUCT_SIGN"],
            "SAME_SIGN": str(a["PRODUCT_SIGN"] == b["PRODUCT_SIGN"]).upper(),
        })
    convergence_fields = ["LOC", "CHR", "BLOCK_START", "BLOCK_STOP", "EXTERNAL_LEAD_SNP", "EXTERNAL_P", "EXTERNAL_SIGN", "FINNGEN_LEAD_SNP", "FINNGEN_P", "FINNGEN_SIGN", "SAME_SIGN"]
    write_tsv(ROOT / "results/loci/exploratory_same_block_convergence.tsv", convergence_fields, convergence)

    write_tsv(
        ROOT / "results/loci/position_mapping_summary.tsv",
        ["requested_snps", "mapped_snps", "unmapped_snps"],
        [{"requested_snps": len(needed_snps), "mapped_snps": len(positions), "unmapped_snps": len(unmapped)}],
    )
    if unmapped:
        (ROOT / "results/loci/unmapped_snps.txt").write_text("\n".join(sorted(unmapped)) + "\n")


if __name__ == "__main__":
    main()
