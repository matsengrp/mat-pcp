#!/usr/bin/env bash
# v10 production run: the gisaid-and-public_2026-04-22 tree only (EPIC #138).
#
# Deliberately NOT `set -e` -- the snakemake exit code is captured and reported
# rather than silently ending the script, so a detached tmux pane still shows
# what happened after the fact.
set -uo pipefail

cd /fh/fast/matsen_e/shared/viral-dasm/mat-pcp || exit 1

CONFIG=config/config-2026-04-22-gisaid-and-public-six-genes.yaml
# Derived from VERSION rather than written out, so the "expect" line below cannot
# disagree with the directory Snakemake actually writes. The v9 copy of this
# script hardcoded it, which is exactly the drift #88 is about.
VERSION=$(cat VERSION)
STAMP=$(date +%Y%m%d-%H%M%S)
LOG="run-${VERSION}-2026-04-22-${STAMP}.log"

{
  echo "=== START      $(date -Is)"
  echo "=== host       $(hostname)"
  echo "=== branch     $(git rev-parse --abbrev-ref HEAD)"
  echo "=== commit     $(git rev-parse --short HEAD)"
  echo "=== VERSION    ${VERSION}"
  echo "=== config     ${CONFIG}"
  echo "=== cores      16"
  echo "=== expect     results/gisaid-and-public_2026-04-22_${VERSION}/"
  echo "=== new in v10 glm/master_table_{pre_omicron,omicron}.csv,"
  echo "===            rates_table_{pre_omicron,omicron}.csv, glm/era_fit_summary.csv"
} | tee -a "$LOG"

pixi run snakemake -j 16 --use-conda --configfile "$CONFIG" 2>&1 | tee -a "$LOG"
EXIT=${PIPESTATUS[0]}

echo "=== SNAKEMAKE EXIT: ${EXIT}   $(date -Is) ===" | tee -a "$LOG"
echo "${EXIT}" > "run-${VERSION}-2026-04-22-${STAMP}.exit"
echo "=== log: ${LOG}" | tee -a "$LOG"

# Post-run sanity check: mat-pcp's per-era rates against upstream's published
# per-era master tables (#137 success criterion 4). Era against matching era --
# the pooled-vs-both comparison in analysis/137-era-offsets/upstream/README.md
# was the pre-Task-2 stand-in for this.
if [ "${EXIT}" -eq 0 ]; then
  RUN="results/gisaid-and-public_2026-04-22_${VERSION}"
  echo "=== per-era rate sanity check vs SARS2-mut-fitness-v2 ===" | tee -a "$LOG"
  pixi run python workflow/scripts/compare_rates_to_upstream_eras.py \
    --pre-omicron "${RUN}/glm/master_table_pre_omicron.csv" \
    --omicron     "${RUN}/glm/master_table_omicron.csv" \
    --outdir      analysis/137-era-offsets/upstream 2>&1 | tee -a "$LOG"
fi
