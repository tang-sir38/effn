#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_all.sh -- rerun every stage that produces a number in the paper.
#
# Prerequisites
#   * effn is installed          R CMD INSTALL ../effn   (or effn_1.1.0.tar.gz)
#   * the GEO files are in place $EFFN_ROOT/GSE96583 (see data/README.md)
#   * EFFN_ROOT points at your working directory (default /root/effn)
#
# Usage
#   bash run_all.sh            # full pipeline
#   bash run_all.sh fast       # skip the calibration arm and the Hierarchicell run
# ---------------------------------------------------------------------------
set -u

MODE="${1:-full}"
cd "$(dirname "$0")" || exit 1

: "${EFFN_ROOT:=/root/effn}"
: "${EFFN_OUT:=$EFFN_ROOT/out}"
export EFFN_ROOT EFFN_OUT
mkdir -p "$EFFN_OUT" "$EFFN_ROOT/logs"

RSCRIPT="${RSCRIPT:-Rscript}"
if ! command -v "$RSCRIPT" >/dev/null 2>&1; then
  echo "Rscript not found. Set RSCRIPT=/full/path/to/Rscript" >&2
  exit 1
fi

if ! "$RSCRIPT" -e 'quit(status = !requireNamespace("effn", quietly = TRUE))'; then
  echo "package 'effn' is not installed. Run: R CMD INSTALL ../effn" >&2
  exit 1
fi

run () {                      # run <label> <script> [args...]
  local label="$1"; shift
  echo "=== $label ==="
  "$RSCRIPT" "$@" > "$EFFN_ROOT/logs/${label}.log" 2>&1
  local rc=$?
  echo "    exit=$rc  (log: $EFFN_ROOT/logs/${label}.log)"
  [ "$rc" -ne 0 ] && echo "    !! stage failed, stopping" && exit "$rc"
  return 0
}

run 01_build_gse96583      01_build_gse96583.R
run 02_effn_main           02_effn_main.R
run 03_depth_and_export    03_depth_and_export.R
run 07_unpaired_check      07_unpaired_check.R
run 10_small_dataset       10_small_dataset.R
run 13_quantile_gap        13_quantile_gap.R "$EFFN_OUT/tables"

if [ "$MODE" != "fast" ]; then
  run 04_sim_calibration   04_sim_calibration.R
  run 05_subsample_recheck 05_subsample_recheck.R
  run 12_small_ci          12_small_ci.R
  # The next two need Hierarchicell 1.0.0 and are skipped if it is absent.
  if "$RSCRIPT" -e 'quit(status = !requireNamespace("Hierarchicell", quietly = TRUE))'; then
    run 08_hierarchicell   08_hierarchicell_compare.R
    run 09_hc_power_curve  09_hc_power_curve.R
  else
    echo "=== 08/09 skipped: Hierarchicell not installed ==="
  fi
fi

# Figure 1 needs the stage 02 and 03 outputs.
if [ -f "$EFFN_OUT/gse96583_summary.json" ]; then
  run 06_figure1           06_figure1.R
fi

echo "ALL STAGES DONE -- outputs in $EFFN_OUT"
