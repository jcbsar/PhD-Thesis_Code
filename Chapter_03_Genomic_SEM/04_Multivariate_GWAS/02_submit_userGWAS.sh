
#!/usr/bin/env bash
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --time=24:00:00

set -euo pipefail

# ============================================================
# Chapter 3 — Multivariate GWAS
# Submit chromosome-level Genomic SEM analyses
#
# Usage:
#   bash 02_submit_userGWAS.sh ROOT ANALYSIS --dry-run
#   bash 02_submit_userGWAS.sh ROOT ANALYSIS
#
# ANALYSIS: GGE | FE | GGE_no_CAE
#
# GGE / FE: 22 chromosome jobs, 4 CPUs each
# GGE_no_CAE: 23 jobs, 1 CPU each
#   Chromosome 2 is processed in two chunks.
#
# The launcher calls 02_run_userGWAS.R.
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="$SCRIPT_DIR/02_run_userGWAS.R"

# ------------------------------------------------------------
# 1. Execution inside a SLURM array task
# ------------------------------------------------------------

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

  : "${ROOT:?Missing ROOT environment variable}"
  : "${ANALYSIS:?Missing ANALYSIS environment variable}"

  task="$SLURM_ARRAY_TASK_ID"

  export R_LIBS_USER="$HOME/R/4.6"

  if [[ "$ANALYSIS" == "GGE_no_CAE" ]]; then

    cores=1

    if (( task == 2 || task == 3 )); then

      # Chromosome 2 is divided into two equal chunks.
      chr=2
      part=$((task - 1))
      parts=2

      exec Rscript "$RUNNER" \
        "$ROOT" "$ANALYSIS" "$chr" "$cores" \
        "$part" "$parts"

    elif (( task >= 4 )); then
      chr=$((task - 1))
    else
      chr=1
    fi

  else
    cores=4
    chr="$task"
  fi

  exec Rscript "$RUNNER" \
    "$ROOT" "$ANALYSIS" "$chr" "$cores"
fi

# ------------------------------------------------------------
# 2. Submission from the HPC login node
# ------------------------------------------------------------

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash $0 ROOT ANALYSIS [--dry-run]"
  exit 1
fi

ROOT="$(realpath "$1")"
ANALYSIS="$2"
MODE="${3:-submit}"

if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: Root directory does not exist."
  exit 1
fi

if [[ ! -f "$RUNNER" ]]; then
  echo "ERROR: Missing runner: $RUNNER"
  exit 1
fi

case "$ANALYSIS" in
  GGE|FE)
    array="1-22%4"
    cpus=4
    memory="48G"
    ;;
  GGE_no_CAE)
    array="1-23%4"
    cpus=1
    memory="48G"
    ;;
  *)
    echo "ERROR: Invalid analysis: $ANALYSIS"
    exit 1
    ;;
esac

echo "========================================"
echo "Genomic SEM multivariate GWAS"
echo "========================================"
echo "Analysis: $ANALYSIS"
echo "Root: $ROOT"
echo "Array: $array"
echo "CPUs per task: $cpus"
echo "Memory per task: $memory"
echo "Time limit: 24 hours"
echo "Runner: $RUNNER"

# ------------------------------------------------------------
# 3. Safe preview or actual submission
# ------------------------------------------------------------

if [[ "$MODE" == "--dry-run" ]]; then
  echo
  echo "DRY RUN: No jobs submitted."
  exit 0
elif [[ "$MODE" != "submit" ]]; then
  echo "ERROR: Unknown option: $MODE"
  exit 1
fi

mkdir -p "$ROOT/logs"

sbatch \
  --job-name="userGWAS_${ANALYSIS}" \
  --array="$array" \
  --cpus-per-task="$cpus" \
  --mem="$memory" \
  --output="$ROOT/logs/${ANALYSIS}_%A_%a.out" \
  --error="$ROOT/logs/${ANALYSIS}_%A_%a.err" \
  --export="ALL,ROOT=$ROOT,ANALYSIS=$ANALYSIS" \
  "$0"
