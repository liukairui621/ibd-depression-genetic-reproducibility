#!/usr/bin/env python3
"""Execute one explicitly selected stage from an isolated workspace."""
import argparse
import datetime
import json
import os
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--project", required=True, type=Path)
p.add_argument("stage")
args = p.parse_args()
project = args.project.resolve()
commands = json.loads((project / "reproducibility/run_order.json").read_text())
if args.stage not in commands:
    raise SystemExit(f"Choose a stage from {list(commands)}")
logdir = project / "logs/portable_execution"
logdir.mkdir(parents=True, exist_ok=True)
stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
record = {"stage": args.stage, "command": commands[args.stage], "started_utc": stamp}
env = os.environ.copy()
env["R_LIBS_USER"] = str(project / "software/R_lib")
provenance = json.loads((project / "provenance/PORTABLE_WORKSPACE.json").read_text())
python3 = provenance["path_replacements"]["/root/anaconda3/bin/python"]
env["PATH"] = str(Path(python3).parent) + os.pathsep + env.get("PATH", "")
with (logdir / f"{args.stage}_{stamp}.log").open("w") as log:
    result = subprocess.run(commands[args.stage], cwd=project, env=env, stdout=log, stderr=subprocess.STDOUT)
record["exit_code"] = result.returncode
record["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
(logdir / f"{args.stage}_{stamp}.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
raise SystemExit(result.returncode)
