#!/usr/bin/env bash
# Reconstructed from archived August 2026 standalone LDSC v1.0.1 execution logs.
# Performs 13 ordered multi-trait --rg commands, comprising 91 unique pairs.
# This is NOT a verbatim recovered original submission script.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MANIFEST="$ROOT/config/trait_prevalence_specification.tsv"
: "${LDSC_PY:?Set LDSC_PY to the full path of ldsc.py}"
: "${MUNGED_DIR:?Set MUNGED_DIR to the directory containing .hm3.sumstats.gz files}"
: "${EUR_LD_DIR:?Set EUR_LD_DIR to the European eur_w_ld_chr directory}"
: "${ANALYSIS_OUT_DIR:?Set ANALYSIS_OUT_DIR to a private working output directory}"
PYTHON_BIN="${PYTHON_BIN:-python}"

[[ -f "$LDSC_PY" && -d "$EUR_LD_DIR" ]] || { echo 'Missing LDSC script or reference directory' >&2; exit 1; }
mkdir -p "$ANALYSIS_OUT_DIR/rg"
traits=(); files=(); sample=(); population=()
while IFS=$'\t' read -r trait filename sample_prev population_prev; do
  [[ "$trait" == "Trait" || -z "$trait" ]] && continue
  file="$MUNGED_DIR/$filename"
  [[ -f "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
  traits+=("$trait"); files+=("$file")
  sample+=("$sample_prev"); population+=("$population_prev")
done < "$MANIFEST"
[[ "${#traits[@]}" -eq 14 ]] || { echo 'Expected exactly 14 phenotypes' >&2; exit 1; }

# For each phenotype, LDSC estimates correlations with subsequent phenotypes.
# The final (14th) phenotype requires no separate log because it has no later pair.
for ((i=0; i<${#traits[@]}-1; i++)); do
  input_files=$(IFS=,; echo "${files[*]:i}")
  input_sample=$(IFS=,; echo "${sample[*]:i}")
  input_pop=$(IFS=,; echo "${population[*]:i}")
  echo "Running rg from ${traits[i]} (${#traits[@]} total phenotypes)"
  "$PYTHON_BIN" "$LDSC_PY" \
    --rg "$input_files" \
    --samp-prev "$input_sample" --pop-prev "$input_pop" \
    --ref-ld-chr "$EUR_LD_DIR/" --w-ld-chr "$EUR_LD_DIR/" \
    --out "$ANALYSIS_OUT_DIR/rg/${traits[i]}_rg"
done
echo "Completed. See $ANALYSIS_OUT_DIR/rg"
