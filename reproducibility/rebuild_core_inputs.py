#!/usr/bin/env python3
"""Rebuild the six core LDSC inputs from source GWAS, without legacy caches."""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import subprocess
from pathlib import Path


def content_hash(path):
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_delange(source, target):
    # These named columns reproduce the archived five-column extraction.
    fields = ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "p_value"]
    with gzip.open(source, "rt", newline="") as src, gzip.open(target, "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        missing = set(fields) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"de Lange schema missing {sorted(missing)}")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "BETA", "P"])
        for row in reader:
            writer.writerow([row[x] for x in fields])


def build_alleles(ref, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as dst:
        dst.write("SNP A1 A2\n")
        for chrom in range(1, 23):
            with (ref / f"1000G.EUR.QC.{chrom}.bim").open() as src:
                for line in src:
                    fields = line.split()
                    dst.write(f"{fields[1]} {fields[4]} {fields[5]}\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ibd-root", required=True, type=Path)
    p.add_argument("--python2", required=True)
    p.add_argument("--ldsc-dir", required=True, type=Path)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--traits", nargs="*")
    args = p.parse_args()
    root = args.ibd_root.resolve()
    project = root / "20_Reproducibility_Ladder"
    out = project / "data/munged_core"
    out.mkdir(parents=True, exist_ok=True)
    logs = project / "logs/rebuilt_inputs"
    logs.mkdir(parents=True, exist_ok=True)
    reference = root / "00_RawData/Reference/LD_EUR"
    alleles = reference / "w_hm3_alleles.snplist"
    build_alleles(reference / "1000G_EUR_Phase3_plink", alleles)
    raw = root / "00_RawData/GWAS"
    specs = {
        "MDD_Howard2019": (raw / "Psychiatric/MDD_Howard2019_new.txt.gz", ["--snp", "MarkerName", "--a1", "A1", "--a2", "A2", "--p", "P", "--frq", "Freq", "--signed-sumstats", "LogOR,0", "--N-cas", "170756", "--N-con", "329443"]),
        "DEP_FinnGen_R12": (root / "19_DENND1B_IndependentValidation/data/raw/finngen_R12_F5_DEPRESSIO.gz", None),
        "IBD_FinnGen_R12": (raw / "IBD/IBD_STRICT_FinnGen_R12.gz", None),
        "CD_FinnGen_R12": (raw / "IBD/CD_STRICT2_FinnGen_R12.gz", None),
        "UC_FinnGen_R12": (raw / "IBD/UC_STRICT2_FinnGen_R12.gz", None),
        "IBD_deLange2017": (raw / "IBD/deLange2017/deLange2017_harmonised.tsv.gz", ["--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P", "--signed-sumstats", "BETA,0", "--N", "59957"]),
    }
    counts = {"DEP_FinnGen_R12": (59333, 434831), "IBD_FinnGen_R12": (10960, 489388), "CD_FinnGen_R12": (2489, 497622), "UC_FinnGen_R12": (7220, 492160)}
    selected = args.traits or list(specs)
    if set(selected) - set(specs):
        raise ValueError("Unknown trait requested")

    def run(trait):
        source, flags = specs[trait]
        if not source.is_file():
            raise FileNotFoundError(source)
        target = out / f"{trait}.sumstats.gz"
        if target.exists() or target.is_symlink():
            raise FileExistsError(f"Refusing to replace existing analysis input: {target}")
        if trait == "IBD_deLange2017":
            prepared = out / "IBD_deLange2017.normalized.tsv.gz"
            normalize_delange(source, prepared)
            source = prepared
        if flags is None:
            cases, controls = counts[trait]
            flags = ["--snp", "rsids", "--a1", "alt", "--a2", "ref", "--p", "pval", "--frq", "af_alt", "--signed-sumstats", "beta,0", "--N-cas", str(cases), "--N-con", str(controls)]
        cmd = [args.python2, str(args.ldsc_dir / "munge_sumstats.py"), "--sumstats", str(source), *flags, "--merge-alleles", str(alleles), "--chunksize", "500000", "--out", str(out / trait)]
        with (logs / f"{trait}.stdout.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        record = {"trait": trait, "source": str(source), "command": cmd, "output": str(target), "uncompressed_sha256": content_hash(target)}
        (logs / f"{trait}.json").write_text(json.dumps(record, indent=2) + "\n")
        print(json.dumps(record), flush=True)
        return record

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        records = list(pool.map(run, selected))
    aliases = {
        "IBD_deLange2017": "16_QualityUpgrade/01_LDSC_deLange_MDD/munged/deLange2017_IBD.sumstats.gz",
        "IBD_FinnGen_R12": "03_LDSC/munged/IBD_STRICT_FinnGen.sumstats.gz",
        "CD_FinnGen_R12": "03_LDSC/munged/CD_STRICT2_FinnGen.sumstats.gz",
        "UC_FinnGen_R12": "03_LDSC/munged/UC_STRICT2_FinnGen.sumstats.gz",
        "MDD_Howard2019": "20_Reproducibility_Ladder/corrections/20260806_howard_no23andme/work/munged/MDD_Howard2019_no23andMe.sumstats.gz",
    }
    # Aliases satisfy archived drivers; every target was rebuilt above.
    for trait, relative in aliases.items():
        if trait not in selected:
            continue
        alias = root / relative
        alias.parent.mkdir(parents=True, exist_ok=True)
        if alias.exists() or alias.is_symlink():
            raise FileExistsError(alias)
        alias.symlink_to(out / f"{trait}.sumstats.gz")
    (logs / "REBUILD_MANIFEST.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
