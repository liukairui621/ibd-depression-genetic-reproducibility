#!/usr/bin/env python3
import csv
import gzip
import hashlib
import subprocess
from pathlib import Path


ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
REQUIRED = [
    ROOT / "provenance/extension_prespecified_protocol.md",
    ROOT / "provenance/extension_source_manifest.tsv",
    ROOT / "reports/MANUSCRIPT_GLOBAL_ROBUSTNESS_INTEGRATION.md",
    ROOT / "provenance/extension_download_checksums.tsv",
    ROOT / "provenance/extension_munged_manifest.tsv",
    ROOT / "provenance/extension_pair_manifest.tsv",
    ROOT / "results/global_extension/h2_summary.tsv",
    ROOT / "results/global_extension/rg_cross_disease.tsv",
    ROOT / "results/global_extension/rg_ukb_definition_matrix.tsv",
    ROOT / "results/global_extension/genomicsem/genomicsem_ldsc.rds",
    ROOT / "results/global_extension/genomicsem/parameter_index.tsv",
    ROOT / "results/global_extension/attribution/dispersion_ratio.tsv",
    ROOT / "results/global_extension/attribution/leave_one_subtype.tsv",
    ROOT
    / "results/global_extension/attribution/cd_uc_direction_by_cohort.tsv",
    ROOT
    / "results/global_extension/attribution/instability_attribution_decision.tsv",
    ROOT
    / "results/global_extension/attribution/sampling_covariance_aware_contrasts.tsv",
    ROOT / "figures/extension/FigExt1A_UKB_definition_by_IBD_subtype_rg_heatmap.pdf",
    ROOT / "figures/extension/FigExt1B_adjusted_marginal_rg_forest.pdf",
    ROOT / "figures/extension/FigExt1C_instability_dispersion_comparison.pdf",
    ROOT / "reports/EXTENSION_UKB_DEPRESSION_CD_UC_REPORT.md",
]

for path in REQUIRED:
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required output: {path}")

with (ROOT / "results/global_extension/rg_cross_disease.tsv").open(
    newline=""
) as handle:
    cross_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(cross_rows) != 36:
    raise RuntimeError(f"Expected 36 cross-disease rows, observed {len(cross_rows)}")

with (ROOT / "results/global_extension/rg_ukb_definition_matrix.tsv").open(
    newline=""
) as handle:
    primary_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(primary_rows) != 24:
    raise RuntimeError(f"Expected 24 UKB extension rows, observed {len(primary_rows)}")

raw_gzip = [
    ROOT / "data/raw/UKB_LifetimeMDD.covararraypcs.assoc.logistic.gz",
    ROOT / "data/raw/UKB_MDDRecur.covararraypcs.assoc.logistic.gz",
    ROOT / "data/raw/UKB_GPpsy.covararraypcs.assoc.logistic.gz",
    ROOT / "data/raw/UKB_ICD10Dep.covararraypcs.assoc.logistic.gz",
    ROOT / "data/raw/CD_deLange2017.h.tsv.gz",
    ROOT / "data/raw/UC_deLange2017.h.tsv.gz",
]
for path in raw_gzip:
    with gzip.open(path, "rb") as handle:
        while handle.read(1024 * 1024):
            pass

pdfinfo = []
for path in REQUIRED:
    if path.suffix == ".pdf":
        result = subprocess.run(
            ["pdfinfo", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        page_line = next(
            line for line in result.stdout.splitlines() if line.startswith("Pages:")
        )
        pdfinfo.append((path.name, path.stat().st_size, page_line.split(":")[1].strip()))

manifest = ROOT / "provenance/extension_final_sha256.tsv"
versions = ROOT / "provenance/extension_software_versions.txt"
commands = [
    ["date", "--iso-8601=seconds"],
    ["uname", "-a"],
    ["R", "--version"],
    ["/root/anaconda3/envs/ldsc/bin/python", "--version"],
]
lines = []
for command in commands:
    completed = subprocess.run(
        command, check=True, capture_output=True, text=True
    )
    lines.append((completed.stdout or completed.stderr).strip())
lines.extend(
    [
        "LDSC_script_md5="
        + hashlib.md5(Path("/root/ldsc/ldsc.py").read_bytes()).hexdigest(),
        "munge_sumstats_script_md5="
        + hashlib.md5(
            Path("/root/ldsc/munge_sumstats.py").read_bytes()
        ).hexdigest(),
        "GenomicSEM_version=0.0.5",
        "GenomicSEM_RemoteSha=0a63ac0ea01b61d28bd17e4a204e0fa561ce5040",
        "random_seed=42",
        "GenomicSEM_n_blocks=200",
    ]
)
versions.write_text("\n".join(lines) + "\n")

audit = ROOT / "reports/EXTENSION_FINAL_AUDIT.md"
audit_lines = [
    "# Extension final audit",
    "",
    f"- Required outputs checked: {len(REQUIRED)}",
    "- Cross-disease LDSC rows: 36",
    "- Primary UKB-definition x IBD-phenotype rows: 24",
    "- Raw gzip streams: all passed full decompression",
    f"- Checksum manifest: `{manifest}`",
    "",
    "## PDF validation",
    "",
    "| File | Bytes | Pages |",
    "|---|---:|---:|",
]
for name, size, pages in pdfinfo:
    audit_lines.append(f"| {name} | {size} | {pages} |")
audit.write_text("\n".join(audit_lines) + "\n")

(ROOT / "provenance/extension_completed_at.txt").write_text(
    subprocess.run(
        ["date", "--iso-8601=seconds"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
)

manifest_paths = [ROOT / "README.md"]
for directory in [
    ROOT / "scripts",
    ROOT / "provenance",
    ROOT / "results/global_extension",
    ROOT / "figures/extension",
    ROOT / "reports",
    ROOT / "logs",
]:
    for path in directory.rglob("*"):
        if (
            path.is_file()
            and "partial_downloads" not in str(path)
            and path != manifest
        ):
            manifest_paths.append(path)

with manifest.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["relative_path", "bytes", "sha256"])
    for path in sorted(set(manifest_paths)):
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        writer.writerow(
            [str(path.relative_to(ROOT)), path.stat().st_size, digest.hexdigest()]
        )

print(audit)
print(manifest)
