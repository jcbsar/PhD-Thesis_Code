
#!/usr/bin/env bash
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00

set -euo pipefail

# ============================================================
# Chapter 3 — Multivariate GWAS
# 03. Submit chromosome-result combination
#
# Usage:
#   bash 03_submit_combine.sh ROOT ANALYSIS --dry-run
#   bash 03_submit_combine.sh ROOT ANALYSIS
#
# ANALYSIS: GGE | FE | GGE_no_CAE
#
# Run after all chromosome-level GWAS jobs finish.
# Uses 03_combine_results.R to validate and combine outputs.
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMBINER="$SCRIPT_DIR/03_combine_results.R"

# ------------------------------------------------------------
# 1. Execution inside the SLURM job
# ------------------------------------------------------------

if [[ -n "${SLURM_JOB_ID:-}" ]]; then

  : "${ROOT:?Missing ROOT environment variable}"
  : "${ANALYSIS:?Missing ANALYSIS environment variable}"

  export R_LIBS_USER="$HOME/R/4.6"

  echo "========================================"
  echo "Combining Genomic SEM results"
  echo "========================================"
  echo "Analysis: $ANALYSIS"
  echo "Root: $ROOT"
  echo "Started: $(date)"

  Rscript "$COMBINER" "$ROOT" "$ANALYSIS"

  echo "Completed: $(date)"
  exit 0
fi

# ------------------------------------------------------------
# 2. Validate submission arguments
# ------------------------------------------------------------

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash $0 ROOT ANALYSIS [--dry-run]"
  exit 1
fi

if [[ ! -d "$1" ]]; then
  echo "ERROR: Root directory does not exist: $1"
  exit 1
fi

ROOT="$(realpath "$1")"
ANALYSIS="$2"
MODE="${3:-submit}"

case "$ANALYSIS" in
  GGE|FE|GGE_no_CAE)
    ;;
  *)
    echo "ERROR: Invalid analysis: $ANALYSIS"
    exit 1
    ;;
esac

if [[ ! -f "$COMBINER" ]]; then
  echo "ERROR: Combination script not found: $COMBINER"
  exit 1
fi

if [[ "$MODE" != "submit" && "$MODE" != "--dry-run" ]]; then
  echo "ERROR: Invalid option: $MODE"
  exit 1
fi

RESULT_DIR="$ROOT/results/consolidated/$ANALYSIS"
FINAL_DIR="$ROOT/final_consolidated/$ANALYSIS"

echo "========================================"
echo "Genomic SEM result combination"
echo "========================================"
echo "Analysis: $ANALYSIS"
echo "Root: $ROOT"
echo "Input: $RESULT_DIR"
echo "Output: $FINAL_DIR"
echo "Memory: 24G"
echo "Time limit: 8 hours"
echo "Script: $COMBINER"

# ------------------------------------------------------------
# 3. Dry run or submit
# ------------------------------------------------------------

if [[ "$MODE" == "--dry-run" ]]; then
  echo
  echo "DRY RUN: No job submitted."
  exit 0
fi

if [[ ! -d "$RESULT_DIR" ]]; then
  echo "ERROR: Chromosome-result directory does not exist."
  exit 1
fi

if [[ -e "$FINAL_DIR" ]]; then
  echo "ERROR: Final output already exists."
  echo "Refusing to overwrite: $FINAL_DIR"
  exit 1
fi

mkdir -p "$ROOT/logs"

sbatch \
  --job-name="combine_${ANALYSIS}" \
  --output="$ROOT/logs/combine_${ANALYSIS}_%j.out" \
  --error="$ROOT/logs/combine_${ANALYSIS}_%j.err" \
  --export="ALL,ROOT=$ROOT,ANALYSIS=$ANALYSIS" \
  "$0"
