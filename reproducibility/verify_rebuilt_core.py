#!/usr/bin/env python3
"""Compare a reconstructed workspace with frozen uncompressed input hashes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--project", type=Path, required=True)
args = p.parse_args()
expected = json.loads((Path(__file__).parent / "FROZEN_CORE_INPUT_HASHES.json").read_text())
checks = []
for entry in expected:
    path = args.project / "data/munged_core" / (entry["trait"] + ".sumstats.gz")
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    actual = digest.hexdigest()
    checks.append({"trait": entry["trait"], "sha256": actual, "matches": actual == entry["uncompressed_sha256"]})
print(json.dumps({"scope": "Decompressed input-content comparison", "checks": checks}, indent=2))
if not all(row["matches"] for row in checks):
    raise SystemExit(1)
