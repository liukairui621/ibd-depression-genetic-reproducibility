#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov
R_SCRIPT="$ROOT/scripts/01_run_sldsc_pair.R"

Rscript "$R_SCRIPT" \
  discovery \
  /root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_deLange2017.sumstats.gz \
  /root/IBD/20_Reproducibility_Ladder/data/munged_core/DEP_FinnGen_R12.sumstats.gz \
  IBD_deLange2017 \
  DEP_FinnGen_R12 \
  > "$ROOT/logs/discovery.stdout.log" 2> "$ROOT/logs/discovery.stderr.log"

Rscript "$R_SCRIPT" \
  replication \
  /root/IBD/20_Reproducibility_Ladder/data/munged_core/IBD_FinnGen_R12.sumstats.gz \
  /root/IBD/20_Reproducibility_Ladder/corrections/20260806_howard_no23andme/work/munged/MDD_Howard2019_no23andMe.sumstats.gz \
  IBD_FinnGen_R12 \
  DEP_Howard2019_no23andMe \
  > "$ROOT/logs/replication.stdout.log" 2> "$ROOT/logs/replication.stderr.log"

Rscript "$ROOT/scripts/02_summarize_results.R" \
  > "$ROOT/logs/summarize.stdout.log" 2> "$ROOT/logs/summarize.stderr.log"

