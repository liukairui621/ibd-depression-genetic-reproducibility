#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


REGIONS = {
    "block154": (1, 200_134_006, 201_067_952, ["rs169850", "rs3861929"]),
    "block464": (3, 47_588_462, 50_387_742, ["rs9862080"]),
    "block1671": (11, 60_515_106, 61_717_117, ["rs174581"]),
}
REF_DIR = Path("/root/IBD/00_RawData/Reference/LD_EUR/1000G_EUR_Phase3_plink")
HOWARD = Path("/root/IBD/00_RawData/GWAS/Psychiatric/MDD_Howard2019_new.txt.gz")
DELANGE = Path("/root/IBD/00_RawData/GWAS/IBD/deLange2017/deLange2017_harmonised.tsv.gz")
FINNGEN_DEP = Path("/root/IBD/19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz")
FINNGEN_IBD = Path("/root/IBD/00_RawData/GWAS/IBD/IBD_STRICT_FinnGen_R12.gz")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}


def complement(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.translate(COMPLEMENT)


def load_reference(chrom: int, start: int, stop: int) -> pd.DataFrame:
    path = REF_DIR / f"1000G.EUR.QC.{chrom}.bim"
    ref = pd.read_csv(
        path, sep=r"\s+", header=None,
        names=["chr", "rsid", "cm", "pos", "bim_a1", "bim_a2"],
        dtype={"rsid": str, "bim_a1": str, "bim_a2": str},
    )
    ref = ref[
        ref["pos"].between(start, stop)
        & ref["rsid"].str.startswith("rs", na=False)
        & ref["bim_a1"].str.fullmatch("[ACGTacgt]", na=False)
        & ref["bim_a2"].str.fullmatch("[ACGTacgt]", na=False)
    ].copy()
    ref["bim_a1"] = ref["bim_a1"].str.upper()
    ref["bim_a2"] = ref["bim_a2"].str.upper()
    ref = ref[
        [frozenset((a1, a2)) not in PALINDROMIC for a1, a2 in zip(ref.bim_a1, ref.bim_a2)]
    ]
    counts = ref["rsid"].value_counts()
    ref = ref[ref["rsid"].map(counts).eq(1)].copy()
    ref["reference_order"] = np.arange(len(ref))
    return ref


def align(
    data: pd.DataFrame,
    ref: pd.DataFrame,
    rsid_col: str,
    effect_col: str,
    other_col: str,
    beta_col: str,
    se_col: str,
    af_col: str | None,
) -> tuple[pd.DataFrame, dict]:
    cols = [rsid_col, effect_col, other_col, beta_col, se_col]
    if af_col:
        cols.append(af_col)
    x = data[cols].copy()
    x.columns = ["rsid", "effect", "other", "beta", "se"] + (["af"] if af_col else [])
    x["rsid"] = x["rsid"].astype(str)
    x["effect"] = x["effect"].astype(str).str.upper()
    x["other"] = x["other"].astype(str).str.upper()
    for col in ["beta", "se"] + (["af"] if af_col else []):
        x[col] = pd.to_numeric(x[col], errors="coerce")
    x = x[np.isfinite(x.beta) & np.isfinite(x.se) & x.se.gt(0)].copy()
    n_input = len(x)
    x = x.merge(ref[["rsid", "bim_a1", "bim_a2"]], on="rsid", how="inner")
    n_joined = len(x)
    direct_same = x.effect.eq(x.bim_a1) & x.other.eq(x.bim_a2)
    direct_flip = x.effect.eq(x.bim_a2) & x.other.eq(x.bim_a1)
    strand_same = complement(x.effect).eq(x.bim_a1) & complement(x.other).eq(x.bim_a2)
    strand_flip = complement(x.effect).eq(x.bim_a2) & complement(x.other).eq(x.bim_a1)
    compatible = direct_same | direct_flip | strand_same | strand_flip
    x = x[compatible].copy()
    flip = (direct_flip | strand_flip)[compatible]
    x["aligned_beta"] = np.where(flip, -x.beta, x.beta)
    x["aligned_se"] = x.se
    if "af" in x:
        x["aligned_effect_af"] = np.where(flip, 1 - x.af, x.af)
    keep = ["rsid", "aligned_beta", "aligned_se"] + (["aligned_effect_af"] if "af" in x else [])
    x = x[keep].drop_duplicates()
    duplicate_ids = set(x.loc[x.rsid.duplicated(keep=False), "rsid"])
    x = x[~x.rsid.isin(duplicate_ids)].copy()
    return x, {
        "n_numeric_input": n_input,
        "n_reference_joined": n_joined,
        "n_allele_compatible_unique": len(x),
        "n_duplicate_rsids_removed": len(duplicate_ids),
    }


def read_chunks(path: Path, sep, usecols, predicate) -> pd.DataFrame:
    frames = []
    for chunk in pd.read_csv(path, sep=sep, usecols=usecols, chunksize=500_000, low_memory=False):
        selected = predicate(chunk)
        if not selected.empty:
            frames.append(selected)
    if not frames:
        return pd.DataFrame(columns=usecols)
    return pd.concat(frames, ignore_index=True)


def load_howard(rsids: set[str]) -> pd.DataFrame:
    cols = ["MarkerName", "A1", "A2", "Freq", "LogOR", "StdErrLogOR"]
    return read_chunks(HOWARD, r"\s+", cols, lambda x: x[x.MarkerName.astype(str).isin(rsids)])


def load_delange(rsids: set[str]) -> pd.DataFrame:
    cols = ["hm_rsid", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "hm_effect_allele_frequency"]
    return read_chunks(DELANGE, "\t", cols, lambda x: x[x.hm_rsid.astype(str).isin(rsids)])


def load_finngen(path: Path, rsids: set[str]) -> pd.DataFrame:
    cols = ["#chrom", "pos", "ref", "alt", "rsids", "beta", "sebeta", "af_alt"]
    frames = []
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=500_000, low_memory=False):
        expanded = chunk.copy()
        expanded["rsid"] = expanded["rsids"].fillna("").astype(str).str.split(r"[,;]")
        expanded = expanded.explode("rsid")
        expanded["rsid"] = expanded["rsid"].str.strip()
        selected = expanded[expanded.rsid.isin(rsids)]
        if not selected.empty:
            frames.append(selected)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    prepared = args.root / "data" / "prepared"
    results = args.root / "results" / "gwas"
    provenance = args.root / "provenance"
    prepared.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    provenance.mkdir(parents=True, exist_ok=True)

    sample_rows = [
        ["Howard_DEP", "cc", 500_199, 170_756, 329_443],
        ["FinnGen_DEP", "cc", 494_164, 59_333, 434_831],
        ["deLange_IBD", "cc", 59_957, 25_042, 34_915],
        ["FinnGen_IBD", "cc", 500_348, 10_960, 489_388],
    ]
    pd.DataFrame(sample_rows, columns=["trait", "type", "N", "cases", "controls"]).to_csv(
        prepared / "trait_sample_sizes.tsv", sep="\t", index=False
    )

    qc_rows = []
    manifest = []
    for locus, (chrom, start, stop, leads) in REGIONS.items():
        ref = load_reference(chrom, start, stop)
        rsids = set(ref.rsid)
        raw = {
            "Howard_DEP": load_howard(rsids),
            "FinnGen_DEP": load_finngen(FINNGEN_DEP, rsids),
            "deLange_IBD": load_delange(rsids),
            "FinnGen_IBD": load_finngen(FINNGEN_IBD, rsids),
        }
        specs = {
            "Howard_DEP": ("MarkerName", "A1", "A2", "LogOR", "StdErrLogOR", "Freq"),
            "FinnGen_DEP": ("rsid", "alt", "ref", "beta", "sebeta", "af_alt"),
            "deLange_IBD": ("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "hm_effect_allele_frequency"),
            "FinnGen_IBD": ("rsid", "alt", "ref", "beta", "sebeta", "af_alt"),
        }
        aligned = {}
        for trait, frame in raw.items():
            aligned[trait], qc = align(frame, ref, *specs[trait])
            qc_rows.append({"locus": locus, "trait": trait, "n_reference_region": len(ref), **qc})

        common = set(ref.rsid)
        for frame in aligned.values():
            common &= set(frame.rsid)
        wide = ref[ref.rsid.isin(common)].sort_values("reference_order").copy()
        for trait, frame in aligned.items():
            rename = {
                "aligned_beta": f"{trait}_beta",
                "aligned_se": f"{trait}_se",
                "aligned_effect_af": f"{trait}_effect_af",
            }
            wide = wide.merge(frame.rename(columns=rename), on="rsid", how="inner", validate="one_to_one")
            wide[f"{trait}_z"] = wide[f"{trait}_beta"] / wide[f"{trait}_se"]
        wide = wide.sort_values("reference_order")
        wide.drop(columns=["reference_order"], inplace=True)
        wide.to_csv(prepared / f"{locus}_four_gwas_common.tsv.gz", sep="\t", index=False, compression="gzip")
        wide[["rsid"]].to_csv(prepared / f"{locus}_common_snps.txt", index=False, header=False)

        lead_status = {lead: bool((wide.rsid == lead).any()) for lead in leads}
        manifest.append({
            "locus": locus, "chrom": chrom, "start37": start, "stop37": stop,
            "n_reference_nonpal_unique": len(ref), "n_four_gwas_common": len(wide),
            "leads": leads, "lead_present": lead_status,
        })
        if len(wide) < 200:
            raise RuntimeError(f"{locus}: only {len(wide)} four-GWAS common variants")

    pd.DataFrame(qc_rows).to_csv(results / "variant_alignment_qc.tsv", sep="\t", index=False)
    with open(provenance / "gwas_input_manifest.json", "w") as handle:
        json.dump({"regions": manifest, "inputs": [str(HOWARD), str(FINNGEN_DEP), str(DELANGE), str(FINNGEN_IBD)]}, handle, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
