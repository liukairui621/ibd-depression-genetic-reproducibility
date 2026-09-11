#!/usr/bin/env python3
"""Validate staged downstream entry points without running analyses."""
import argparse
import json
import shutil
import subprocess
from pathlib import Path


def resolve_executable(value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        if not path.exists():
            raise FileNotFoundError(f"Configured executable does not exist: {path}")
        return str(path)
    resolved = shutil.which(value)
    if not resolved:
        raise FileNotFoundError(f"Executable is not on PATH: {value}")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    commands = json.loads((project / "reproducibility/run_order.json").read_text())
    checks = []

    for stage, command in commands.items():
        if not command:
            raise ValueError(f"Empty command for stage {stage}")
        executable = resolve_executable(command[0])
        for argument in command[1:]:
            if argument.startswith("-"):
                continue
            candidate = Path(argument)
            if candidate.is_absolute() and candidate.suffix in {".py", ".sh", ".R"} and not candidate.exists():
                raise FileNotFoundError(f"Missing stage script for {stage}: {candidate}")
        checks.append({"stage": stage, "executable": executable, "status": "resolved"})

    for shell_script in project.rglob("*.sh"):
        subprocess.run(["bash", "-n", str(shell_script)], check=True)

    extensions = project / "extensions"
    for extension in extensions.iterdir():
        if not extension.is_dir():
            continue
        for directory in ("logs", "reports", "results", "provenance", "data", "figures", "software"):
            path = extension / directory
            if not path.is_dir():
                raise FileNotFoundError(f"Missing staged extension directory: {path}")

    provenance = json.loads((project / "provenance/PORTABLE_WORKSPACE.json").read_text())
    qtl_key = "/".join(("software", "qtl_env", "bin", "python"))
    qtl_python = provenance["path_replacements"][qtl_key]
    finemap_driver = project / "extensions/20260817_shared_locus_finemap_qtl/scripts/run_all.sh"
    qtl_invocations = [
        line.strip().split()[0]
        for line in finemap_driver.read_text().splitlines()
        if "qtl_env/bin/python" in line
    ]
    if not qtl_invocations or any(invocation != qtl_python for invocation in qtl_invocations):
        raise RuntimeError(
            f"QTL interpreter substitution failed: expected {qtl_python}, observed {qtl_invocations}"
        )

    report = {
        "project": str(project),
        "status": "PASS",
        "scope": "entry-point, interpreter, directory, path-substitution, and shell-syntax checks only",
        "stages": checks,
    }
    output = project / "provenance/PORTABLE_PREFLIGHT.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
