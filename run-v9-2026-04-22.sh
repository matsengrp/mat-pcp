#!/usr/bin/env bash
# v9 production run: the gisaid-and-public_2026-04-22 tree only (EPIC #131 step 11).
#
# Deliberately NOT `set -e` -- the snakemake exit code is captured and reported
# rather than silently ending the script, so a detached tmux pane still shows
# what happened after the fact.
set -uo pipefail

cd /fh/fast/matsen_e/shared/viral-dasm/mat-pcp || exit 1

CONFIG=config/config-2026-04-22-gisaid-and-public-six-genes.yaml
STAMP=$(date +%Y%m%d-%H%M%S)
LOG="run-v9-2026-04-22-${STAMP}.log"

{
  echo "=== START      $(date -Is)"
  echo "=== host       $(hostname)"
  echo "=== branch     $(git rev-parse --abbrev-ref HEAD)"
  echo "=== commit     $(git rev-parse --short HEAD)"
  echo "=== VERSION    $(cat VERSION)"
  echo "=== config     ${CONFIG}"
  echo "=== cores      16"
  echo "=== expect     results/gisaid-and-public_2026-04-22_v9/"
} | tee -a "$LOG"

pixi run snakemake -j 16 --use-conda --configfile "$CONFIG" 2>&1 | tee -a "$LOG"
EXIT=${PIPESTATUS[0]}

echo "=== SNAKEMAKE EXIT: ${EXIT}   $(date -Is) ===" | tee -a "$LOG"
echo "${EXIT}" > "run-v9-2026-04-22-${STAMP}.exit"
echo "=== log: ${LOG}" | tee -a "$LOG"
