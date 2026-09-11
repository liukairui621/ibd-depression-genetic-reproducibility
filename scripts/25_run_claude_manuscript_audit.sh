#!/usr/bin/env bash
set -euo pipefail

ROOT="/root/IBD/20_Reproducibility_Ladder"
REPORT="$ROOT/reports/CLAUDE_CODE_MANUSCRIPT_AUDIT_REPORT.md"
LOG="$ROOT/logs/claude_manuscript_audit_$(date +%Y%m%d_%H%M%S).log"

cd "$ROOT"

logged_in="$(
  claude auth status 2>/dev/null |
    python3 -c 'import json,sys; print(str(json.load(sys.stdin).get("loggedIn", False)).lower())'
)"
if [[ "$logged_in" != "true" ]]; then
  echo "Claude Code is not authenticated. Run 'claude /login' interactively, then rerun this script." >&2
  exit 2
fi

rm -f "$REPORT"

set +e
claude -p \
  "Read CLAUDE_CODE_MANUSCRIPT_REVIEW_HANDOFF.md and follow it exactly. Perform the complete read-only manuscript-level audit. You may create or overwrite only reports/CLAUDE_CODE_MANUSCRIPT_AUDIT_REPORT.md. Do not modify any other file. After writing the report, print a concise completion summary." \
  --model opus \
  --effort high \
  --permission-mode acceptEdits \
  --tools=Read,Bash,Write,Glob,Grep \
  </dev/null 2>&1 | tee "$LOG"
status="${PIPESTATUS[0]}"
set -e

if [[ "$status" -ne 0 ]]; then
  echo "Claude Code exited with status $status. See $LOG." >&2
  exit "$status"
fi
if [[ ! -s "$REPORT" ]]; then
  echo "Claude Code returned successfully but did not create $REPORT." >&2
  exit 3
fi

sha256sum "$REPORT"
echo "Audit report: $REPORT"
echo "Execution log: $LOG"
