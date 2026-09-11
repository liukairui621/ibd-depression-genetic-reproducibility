#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


REGIONS = {
    "block154": (1, 200_134_006, 201_067_952),
    "block464": (3, 47_588_462, 50_387_742),
    "block1671": (11, 60_515_106, 61_717_117),
}
DATASETS = [
    ("GTEx_amygdala_eQTL", "QTS000015", "QTD000146", 129, "eQTL", "all"),
    ("GTEx_amygdala_sQTL", "QTS000015", "QTD000150", 129, "sQTL", "cc"),
    ("GTEx_colon_sigmoid_eQTL", "QTS000015", "QTD000226", 318, "eQTL", "all"),
    ("GTEx_colon_sigmoid_sQTL", "QTS000015", "QTD000230", 318, "sQTL", "cc"),
    ("GTEx_colon_transverse_eQTL", "QTS000015", "QTD000231", 368, "eQTL", "all"),
    ("GTEx_colon_transverse_sQTL", "QTS000015", "QTD000235", 368, "sQTL", "cc"),
    ("GTEx_small_intestine_eQTL", "QTS000015", "QTD000321", 174, "eQTL", "all"),
    ("GTEx_small_intestine_sQTL", "QTS000015", "QTD000325", 174, "sQTL", "cc"),
    ("GTEx_blood_eQTL", "QTS000015", "QTD000356", 670, "eQTL", "all"),
    ("GTEx_blood_sQTL", "QTS000015", "QTD000360", 670, "sQTL", "cc"),
    ("BLUEPRINT_monocyte_eQTL", "QTS000002", "QTD000021", 191, "eQTL", "all"),
    ("BLUEPRINT_monocyte_sQTL", "QTS000002", "QTD000025", 191, "sQTL", "cc"),
    ("BLUEPRINT_CD4_T_eQTL", "QTS000002", "QTD000031", 167, "eQTL", "all"),
    ("BLUEPRINT_CD4_T_sQTL", "QTS000002", "QTD000035", 167, "sQTL", "cc"),
    ("CEDAR_rectum_eQTL", "QTS000007", "QTD000072", 271, "eQTL", "all"),
    ("CEDAR_colon_transverse_eQTL", "QTS000007", "QTD000068", 286, "eQTL", "all"),
    ("CEDAR_ileum_eQTL", "QTS000007", "QTD000074", 180, "eQTL", "all"),
]
HEADER = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
    "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid",
]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TABIX_CACHE: dict[str, pysam.TabixFile] = {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def url(study: str, dataset: str, suffix: str) -> str:
    return f"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{study}/{dataset}/{dataset}.{suffix}.tsv.gz"


def retrieve_region(source: str, chrom: int, start: int, stop: int, gene_ids: set[str]) -> list[list[str]]:
    rows = []
    seen = set()
    if source not in TABIX_CACHE:
        TABIX_CACHE[source] = pysam.TabixFile(source)
    tabix = TABIX_CACHE[source]
    records = tabix.fetch(str(chrom), max(0, start - 500_001), stop + 500_000)
    for line in records:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != len(HEADER):
            continue
        gene = fields[16].split(".")[0]
        if gene not in gene_ids:
            continue
        key = tuple(fields)
        if key not in seen:
            seen.add(key)
            rows.append(fields)
    return rows


def retrieve_region_api(dataset: str, chrom: int, start: int, stop: int, gene_ids: set[str]) -> list[list[str]]:
    rows = []
    seen = set()
    region = f"{chrom}:{max(1, start - 500_000)}-{stop + 500_000}"
    for gene_id in sorted(gene_ids):
        offset = 0
        while True:
            query = urllib.parse.urlencode({
                "size": 1000, "start": offset, "pos": region, "gene_id": gene_id,
            })
            endpoint = f"https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset}/associations?{query}"
            try:
                with urllib.request.urlopen(endpoint, timeout=180) as response:
                    payload = json.load(response)
            except urllib.error.HTTPError as exc:
                if offset > 0 and exc.code in {400, 404}:
                    break
                raise
            if not payload:
                break
            for record in payload:
                fields = [
                    record.get("molecular_trait_id", ""), record.get("chromosome", ""),
                    record.get("position", ""), record.get("ref", ""), record.get("alt", ""),
                    record.get("variant", ""), "", record.get("maf", ""),
                    record.get("pvalue", ""), record.get("beta", ""), record.get("se", ""),
                    record.get("type", ""), record.get("ac", ""), record.get("an", ""),
                    record.get("r2", ""), "", record.get("gene_id", ""),
                    record.get("median_tpm", ""), record.get("rsid", ""),
                ]
                key = tuple(str(value) for value in fields)
                if key not in seen:
                    seen.add(key)
                    rows.append(list(key))
            if len(payload) < 1000:
                break
            offset += 1000
    return rows


def complement(value: str) -> str:
    return value.upper().translate(COMPLEMENT)


def align_qtl(data: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    x = data.copy()
    x["rsid"] = x.rsid.astype(str)
    x["ref"] = x.ref.astype(str).str.upper()
    x["alt"] = x.alt.astype(str).str.upper()
    for col in ["beta", "se", "maf", "pvalue"]:
        x[col] = pd.to_numeric(x[col], errors="coerce")
    x = x[np.isfinite(x.beta) & np.isfinite(x.se) & x.se.gt(0) & np.isfinite(x.pvalue)].copy()
    x = x.merge(ref[["rsid", "bim_a1", "bim_a2", "pos37"]], on="rsid", how="inner")
    direct_same = x.alt.eq(x.bim_a1) & x.ref.eq(x.bim_a2)
    direct_flip = x.alt.eq(x.bim_a2) & x.ref.eq(x.bim_a1)
    strand_same = x.alt.map(complement).eq(x.bim_a1) & x.ref.map(complement).eq(x.bim_a2)
    strand_flip = x.alt.map(complement).eq(x.bim_a2) & x.ref.map(complement).eq(x.bim_a1)
    compatible = direct_same | direct_flip | strand_same | strand_flip
    x = x[compatible].copy()
    flip = (direct_flip | strand_flip)[compatible]
    x["aligned_beta"] = np.where(flip, -x.beta, x.beta)
    x["aligned_se"] = x.se
    x["aligned_effect_af"] = np.where(flip, 1 - x.maf, x.maf)
    # The same variant is expected to recur across genes or splice junctions.
    # Exclude only duplicate records within a single molecular trait.
    counts = x.groupby(["molecular_trait_id", "rsid"])["rsid"].transform("size")
    x = x[counts.eq(1)].copy()
    return x


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    raw_dir = args.root / "data" / "qtl_raw"
    prepared_dir = args.root / "data" / "qtl_prepared"
    out_dir = args.root / "results" / "qtl"
    raw_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.root / "results" / "mapping" / "candidate_genes.tsv", sep="\t")
    retrieval_manifest = []
    all_prepared = []
    trait_manifest = []

    for label, study, dataset, sample_size, qtl_type, suffix in DATASETS:
        source = url(study, dataset, suffix)
        for locus, (chrom, start, stop) in REGIONS.items():
            gene_ids = set(candidates.loc[candidates.locus.eq(locus), "gene_id"].astype(str))
            output = raw_dir / f"{locus}__{label}.tsv.gz"
            status = "ok"
            message = ""
            try:
                if not label.startswith("CEDAR") and output.exists():
                    with gzip.open(output, "rt") as handle:
                        reader = csv.reader(handle, delimiter="\t")
                        existing_header = next(reader)
                        if existing_header != HEADER:
                            raise RuntimeError("existing slice has an unexpected header")
                        rows = list(reader)
                    status = "reused_complete_slice"
                elif label.startswith("CEDAR"):
                    rows = retrieve_region_api(dataset, chrom, start, stop, gene_ids)
                    status = "ok_api"
                else:
                    rows = retrieve_region(source, chrom, start, stop, gene_ids)
            except Exception as exc:
                rows = []
                status = "retrieval_error"
                message = str(exc)
            with gzip.open(output, "wt", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(HEADER)
                writer.writerows(rows)
            retrieval_manifest.append({
                "locus": locus, "dataset_label": label, "study_id": study,
                "dataset_id": dataset, "sample_size": sample_size, "qtl_type": qtl_type,
                "source_url": source, "n_rows": len(rows), "status": status,
                "message": message, "output": str(output), "sha256": sha256(output),
            })
            time.sleep(1.5)

            if not rows:
                continue
            data = pd.DataFrame(rows, columns=HEADER)
            data["gene_id_base"] = data.gene_id.astype(str).str.split(".").str[0]
            data = data[data.gene_id_base.isin(gene_ids)].copy()
            bim = pd.read_csv(
                args.root / "data" / "ld" / f"{locus}_1000G_EUR.bim",
                sep=r"\s+", header=None,
                names=["chr", "rsid", "cm", "pos37", "bim_a1", "bim_a2"],
            )
            bim["bim_a1"] = bim.bim_a1.str.upper()
            bim["bim_a2"] = bim.bim_a2.str.upper()
            aligned = align_qtl(data, bim)
            if aligned.empty:
                continue

            stats = (
                aligned.groupby(["gene_id_base", "molecular_trait_id"], as_index=False)
                .agg(n_aligned=("rsid", "nunique"), min_p=("pvalue", "min"))
            )
            stats["eligible"] = stats.n_aligned.ge(200) & stats.min_p.lt(1e-5)
            for gene_id, group in stats.groupby("gene_id_base"):
                eligible = group[group.eligible].sort_values(["min_p", "n_aligned"], ascending=[True, False])
                if eligible.empty:
                    continue
                if qtl_type == "sQTL":
                    selected_traits = [eligible.iloc[0].molecular_trait_id]
                else:
                    selected_traits = list(eligible.molecular_trait_id)
                for molecular_trait in selected_traits:
                    subset = aligned[
                        aligned.gene_id_base.eq(gene_id)
                        & aligned.molecular_trait_id.eq(molecular_trait)
                    ].copy()
                    qtl_key = f"{locus}__{dataset}__{gene_id}__{molecular_trait}"
                    subset["qtl_key"] = qtl_key
                    subset["locus"] = locus
                    subset["dataset_label"] = label
                    subset["dataset_id"] = dataset
                    subset["qtl_type"] = qtl_type
                    subset["sample_size"] = sample_size
                    all_prepared.append(subset[[
                        "qtl_key", "locus", "dataset_label", "dataset_id", "qtl_type",
                        "sample_size", "gene_id_base", "molecular_trait_id", "rsid", "pos37",
                        "aligned_beta", "aligned_se", "aligned_effect_af", "pvalue",
                    ]])
                    stat = stats[
                        stats.gene_id_base.eq(gene_id)
                        & stats.molecular_trait_id.eq(molecular_trait)
                    ].iloc[0]
                    symbol = candidates.loc[
                        candidates.locus.eq(locus) & candidates.gene_id.eq(gene_id), "gene_symbol"
                    ].iloc[0]
                    trait_manifest.append({
                        "qtl_key": qtl_key, "locus": locus, "dataset_label": label,
                        "dataset_id": dataset, "qtl_type": qtl_type, "sample_size": sample_size,
                        "gene_id": gene_id, "gene_symbol": symbol,
                        "molecular_trait_id": molecular_trait,
                        "n_aligned": int(stat.n_aligned), "min_qtl_p": float(stat.min_p),
                        "n_molecular_traits_screened_for_gene": len(group),
                    })

    pd.DataFrame(retrieval_manifest).to_csv(out_dir / "qtl_retrieval_manifest.tsv", sep="\t", index=False)
    trait_df = pd.DataFrame(trait_manifest)
    trait_df.to_csv(out_dir / "eligible_qtl_traits.tsv", sep="\t", index=False)
    if all_prepared:
        pd.concat(all_prepared, ignore_index=True).to_csv(
            prepared_dir / "qtl_prepared_long.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    else:
        raise RuntimeError("No QTL molecular traits passed the frozen eligibility rules")
    with open(args.root / "provenance" / "qtl_dataset_manifest.json", "w") as handle:
        json.dump(retrieval_manifest, handle, indent=2)
    print(trait_df.to_string(index=False))


if __name__ == "__main__":
    os.environ.setdefault("SSL_CERT_FILE", "/etc/ssl/certs/ca-certificates.crt")
    os.environ.setdefault("CURL_CA_BUNDLE", "/etc/ssl/certs/ca-certificates.crt")
    main()
