#!/usr/bin/env python3
"""Create an LF-delimited checksum manifest for the manuscript audit package."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path


ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
OUTPUT = ROOT / "provenance" / "MANUSCRIPT_AUDIT_PACKAGE_SHA256_20260728.tsv"
TOP_LEVEL_FILES = [
    ROOT / "README.md",
]
DIRECTORIES = [
    ROOT / "manuscript",
    ROOT / "provenance",
    ROOT / "reports",
    ROOT / "results",
    ROOT / "scripts",
    ROOT / "figures",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_entries() -> list[Path]:
    entries = list(TOP_LEVEL_FILES)
    for directory in DIRECTORIES:
        for dirpath, dirnames, filenames in os.walk(directory, followlinks=False):
            dirnames.sort()
            for filename in sorted(filenames):
                entries.append(Path(dirpath) / filename)
    return sorted(set(entries), key=lambda path: str(path.relative_to(ROOT)))


def main() -> None:
    rows = ["relative_path\tkind\tbytes\tsha256\tlink_target"]
    for path in iter_entries():
        if path == OUTPUT:
            continue
        relative = path.relative_to(ROOT).as_posix()
        if path.is_symlink():
            target = os.readlink(path)
            if path.exists() and path.is_file():
                rows.append(
                    f"{relative}\tsymlink\t{path.stat().st_size}\t"
                    f"{sha256_file(path)}\t{target}"
                )
            else:
                rows.append(f"{relative}\tbroken_symlink\tNA\tNA\t{target}")
        elif path.is_file():
            rows.append(
                f"{relative}\tfile\t{path.stat().st_size}\t"
                f"{sha256_file(path)}\t"
            )
    OUTPUT.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(rows) - 1} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
