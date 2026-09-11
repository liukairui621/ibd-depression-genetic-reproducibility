#!/usr/bin/env python3
"""Stage archived drivers in an isolated Linux workspace with explicit assets."""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

PROJECT = "20_Reproducibility_Ladder"
EXTENSIONS = {
    "shared_locus_replication": "20260817_shared_locus_replication",
    "shared_locus_finemap_qtl": "20260817_shared_locus_finemap_qtl",
    "annotation_stratified_gcov": "20260823_annotation_stratified_gcov",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--assets", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.assets.read_text())
    repo = Path(__file__).resolve().parents[1]
    work = args.work_root.resolve()
    if work.exists() and any(work.iterdir()):
        raise SystemExit("Use a new empty work directory; archived results will not be overwritten")
    ibd = work / "IBD"
    project = ibd / PROJECT
    required_tools = ["python2", "python3", "qtl_python", "ldsc_dir", "plink", "placo_dir", "r_library"]
    for key in required_tools:
        if key not in config or not Path(config[key]).exists():
            raise SystemExit(f"Missing configured software path: {key}")
    for relative, source in config["assets"].items():
        if Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise ValueError(f"Asset destination must be relative: {relative}")
        if not Path(source).exists():
            raise FileNotFoundError(f"Source asset {relative}: {source}")
    project.mkdir(parents=True)
    replacements = {
        "/root/anaconda3/envs/ldsc/bin/python": config["python2"],
        "/root/anaconda3/bin/python": config["python3"],
        "/root/anaconda3/bin/plink": config["plink"],
        "/root/ldsc": config["ldsc_dir"],
        "software/qtl_env/bin/python": config["qtl_python"],
        "/root/IBD": str(ibd),
    }
    # Protect configured values from later substitutions. In particular, the
    # configured QTL interpreter may itself live below /root/IBD.
    tokens = {original: f"__PORTABLE_PATH_{index}__" for index, original in enumerate(replacements)}

    def copy_tree(source, destination):
        for path in sorted(source.rglob("*")):
            if not path.is_file() or ".git" in path.parts or "__pycache__" in path.parts:
                continue
            target = destination / path.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            if path.suffix in {".py", ".R", ".sh", ".md", ".json", ".tsv", ".txt"}:
                text = path.read_text().replace("\r\n", "\n")
                for original, token in tokens.items():
                    text = text.replace(original, token)
                for original, token in tokens.items():
                    text = text.replace(token, replacements[original])
                target.write_text(text)
            else:
                shutil.copy2(path, target)

    for subdir in ["scripts", "protocols", "manifests", "reproducibility"]:
        copy_tree(repo / subdir, project / subdir)
    for name, dated in EXTENSIONS.items():
        dest = project / "extensions" / dated
        for subdir in ["scripts", "protocols", "config"]:
            source = repo / "extensions" / name / subdir
            if source.is_dir():
                copy_tree(source, dest / subdir)
        for subdir in ["logs", "reports", "results", "provenance", "data", "figures", "software"]:
            (dest / subdir).mkdir(parents=True, exist_ok=True)
        # The analysis drivers use these two protocol locations.
        for name2 in ["PRESPECIFIED_RULES.md", "RUN_ORDER.md"]:
            source = dest / "protocols" / name2
            if source.exists():
                shutil.copy2(source, dest / "provenance" / name2)
    for subdir in ["logs", "reports", "results", "provenance", "data/raw", "data/prepared", "figures", "software"]:
        (project / subdir).mkdir(parents=True, exist_ok=True)
    for path in (repo / "manifests").glob("*.tsv"):
        shutil.copy2(path, project / "provenance" / path.name)
    for relative, source in config["assets"].items():
        target = ibd / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(Path(source).resolve(), target_is_directory=Path(source).is_dir())
    ref = ibd / "00_RawData/Reference"
    (ref / "LAVA").mkdir(parents=True, exist_ok=True)
    blocks = repo / "reproducibility/resources/LAVA_blocks_fixed.txt"
    if blocks.exists():
        shutil.copy2(blocks, ref / "LAVA/LAVA_blocks_fixed.txt")
    (project / "software/R_lib").symlink_to(Path(config["r_library"]).resolve())
    placo = project / "extensions" / EXTENSIONS["shared_locus_replication"] / "software/PLACO"
    placo.parent.mkdir(parents=True, exist_ok=True)
    placo.symlink_to(Path(config["placo_dir"]).resolve())
    # Keep rebuilt files rather than replacing them with legacy aliases.
    munging = project / "scripts/04_munge_core.sh"
    lines = munging.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("ln -sfn "):
            target = line.rsplit(" ", 1)[-1]
            lines[i] = f"[[ -s {target} ]] || {line}"
    munging.write_text("\n".join(lines) + "\n")
    # Record the transformation; archived scripts remain unchanged in the repository.
    record = {"workspace": str(work), "project": str(project), "path_replacements": replacements, "assets": config["assets"], "source_repository": str(repo)}
    (project / "provenance/PORTABLE_WORKSPACE.json").write_text(json.dumps(record, indent=2) + "\n")
    rebuild = [config["python3"], str(project / "reproducibility/rebuild_core_inputs.py"), "--ibd-root", str(ibd), "--python2", config["python2"], "--ldsc-dir", config["ldsc_dir"]]
    command_file = project / "reproducibility/run_order.json"
    commands = {
        "rebuild_inputs": rebuild,
        "core": ["bash", str(project / "scripts/00_run_phase1_core.sh")],
        "lava": ["bash", str(project / "scripts/12_run_layer2_primary.sh")],
        "phenotype_extension": ["bash", str(project / "scripts/00_run_extension.sh")],
        "placo": ["bash", str(project / "extensions" / EXTENSIONS["shared_locus_replication"] / "scripts/00_run_pipeline.sh")],
        "finemap_qtl": ["bash", str(project / "extensions" / EXTENSIONS["shared_locus_finemap_qtl"] / "scripts/run_all.sh")],
        "annotation": ["bash", str(project / "extensions" / EXTENSIONS["annotation_stratified_gcov"] / "scripts/run_pipeline.sh")],
    }
    command_file.write_text(json.dumps(commands, indent=2) + "\n")
    print(json.dumps({"project": str(project), "run_order": str(command_file), "status": "staged_no_analysis_run"}, indent=2))
    for path in project.rglob("*.sh"):
        subprocess.run(["bash", "-n", str(path)], check=True)
    subprocess.run(
        [config["python3"], str(project / "reproducibility/preflight_downstream.py"), "--project", str(project)],
        check=True,
    )


if __name__ == "__main__":
    main()
