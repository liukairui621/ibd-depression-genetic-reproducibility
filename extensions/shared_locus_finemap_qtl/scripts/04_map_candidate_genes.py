#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import re
from pathlib import Path

import pandas as pd
import requests


GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
REGIONS = {
    "block154": (1, 200_134_006, 201_067_952, ["rs169850", "rs3861929"]),
    "block464": (3, 47_588_462, 50_387_742, ["rs9862080"]),
    "block1671": (11, 60_515_106, 61_717_117, ["rs174581"]),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    tmp = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with open(tmp, "wb") as out:
            for block in response.iter_content(1024 * 1024):
                if block:
                    out.write(block)
    tmp.replace(path)


def attr_value(attributes: str, key: str) -> str:
    match = re.search(rf'{key} "([^"]+)"', attributes)
    return match.group(1) if match else ""


def load_genes(gtf: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0].removeprefix("chr")
            if len(fields) != 9 or fields[2] != "gene" or chrom not in {"1", "3", "11"}:
                continue
            attrs = fields[8]
            rows.append({
                "chr": int(chrom), "gene_start": int(fields[3]), "gene_end": int(fields[4]),
                "strand": fields[6], "gene_id": attr_value(attrs, "gene_id").split(".")[0],
                "gene_symbol": attr_value(attrs, "gene_name"),
                "gene_type": attr_value(attrs, "gene_type"),
            })
    return pd.DataFrame(rows).drop_duplicates("gene_id")


def distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    return min(abs(pos - start), abs(pos - end))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    raw = args.root / "data" / "raw"
    out = args.root / "results" / "mapping"
    raw.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    gtf = raw / "gencode.v19.annotation.gtf.gz"
    download(GTF_URL, gtf)
    genes = load_genes(gtf)

    cs = pd.read_csv(args.root / "results" / "gwas" / "credible_set_members.tsv.gz", sep="\t")
    cs = cs[cs.ridge.eq(0)].copy()
    mapping_rows = []
    selected_variant_rows = []
    for locus, (chrom, start, stop, leads) in REGIONS.items():
        prepared = pd.read_csv(args.root / "data" / "prepared" / f"{locus}_four_gwas_common.tsv.gz", sep="\t")
        lead_table = prepared[prepared.rsid.isin(leads)][["rsid", "chr", "pos"]].copy()
        lead_table["variant_source"] = "PLACO_lead"
        cs_table = cs[cs.locus.eq(locus)][["rsid", "chr", "pos", "trait", "cs_id"]].drop_duplicates().copy()
        cs_table["variant_source"] = "SuSiE_95pct_CS"
        cs_table["variant_source"] += ":" + cs_table.trait.astype(str) + ":CS" + cs_table.cs_id.astype(str)
        cs_table = cs_table.drop(columns=["trait", "cs_id"])
        variants = pd.concat([lead_table, cs_table], ignore_index=True).drop_duplicates()
        selected_variant_rows.append(variants.assign(locus=locus))
        locus_genes = genes[genes.chr.eq(chrom)].copy()
        for variant in variants.itertuples(index=False):
            for gene in locus_genes.itertuples(index=False):
                distance = distance_to_interval(int(variant.pos), int(gene.gene_start), int(gene.gene_end))
                if distance <= 100_000:
                    mapping_rows.append({
                        "locus": locus, "variant": variant.rsid, "variant_pos37": variant.pos,
                        "variant_source": variant.variant_source, "gene_id": gene.gene_id,
                        "gene_symbol": gene.gene_symbol, "gene_type": gene.gene_type,
                        "gene_start37": gene.gene_start, "gene_end37": gene.gene_end,
                        "strand": gene.strand, "distance_bp": distance,
                        "variant_overlaps_gene": distance == 0,
                    })

    mapping = pd.DataFrame(mapping_rows).sort_values(["locus", "distance_bp", "gene_symbol", "variant"])
    mapping.to_csv(out / "candidate_gene_variant_map.tsv", sep="\t", index=False)
    selected = pd.concat(selected_variant_rows, ignore_index=True)
    selected.to_csv(out / "mapping_variant_manifest.tsv", sep="\t", index=False)
    candidates = (
        mapping.groupby(["locus", "gene_id", "gene_symbol", "gene_type", "gene_start37", "gene_end37"], as_index=False)
        .agg(min_distance_bp=("distance_bp", "min"), n_support_variants=("variant", "nunique"),
             overlaps_any_variant=("variant_overlaps_gene", "max"))
        .sort_values(["locus", "min_distance_bp", "gene_symbol"])
    )
    candidates.to_csv(out / "candidate_genes.tsv", sep="\t", index=False)
    with open(args.root / "provenance" / "gencode_v19_sha256.txt", "w") as handle:
        handle.write(f"{sha256(gtf)}  {gtf}\n{GTF_URL}\n")
    print(candidates.to_string(index=False))


if __name__ == "__main__":
    main()
