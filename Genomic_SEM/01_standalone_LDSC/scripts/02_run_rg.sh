
#!/usr/bin/env bash
# Standalone LDSC v1.0.1: pairwise genetic correlation analyses.
# Command parameters and phenotype ordering follow the recorded
# August 2026 LDSC analyses.
#
# Includes seven epilepsy subtypes and seven brain-volume phenotypes.
# Broad GGE and FE are excluded from the correlation analyses.
# Performs 13 ordered --rg commands covering 91 unique pairs.

set -euo pipefail

# Locate the configuration file.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MANIFEST="$ROOT/config/trait_prevalence_specification.tsv"

# Define software, input and output locations.
: "${LDSC_PY:?Set LDSC_PY to the full path of ldsc.py}"
: "${MUNGED_DIR:?Set MUNGED_DIR to the directory containing .hm3.sumstats.gz files}"
: "${EUR_LD_DIR:?Set EUR_LD_DIR to the European LD-score directory}"
: "${ANALYSIS_OUT_DIR:?Set ANALYSIS_OUT_DIR to the output directory}"

PYTHON_BIN="${PYTHON_BIN:-python}"

# Check required inputs.
[[ -f "$LDSC_PY" ]] || {
  echo "Missing LDSC script: $LDSC_PY" >&2
  exit 1
}

[[ -d "$EUR_LD_DIR" ]] || {
  echo "Missing LD-score directory: $EUR_LD_DIR" >&2
  exit 1
}

[[ -f "$MANIFEST" ]] || {
  echo "Missing phenotype manifest: $MANIFEST" >&2
  exit 1
}

mkdir -p "$ANALYSIS_OUT_DIR/rg"

# Read the fourteen phenotypes included in genetic correlations.
traits=()
files=()
sample=()
population=()

while IFS=$'\t' read -r trait filename sample_prev population_prev; do

  [[ "$trait" == "Trait" || -z "$trait" ]] && continue

  # Broad GGE and FE are included in the heritability analyses only.
  if [[ "$trait" == "GGE" || "$trait" == "FE" ]]; then
    continue
  fi

  file="$MUNGED_DIR/$filename"

  [[ -f "$file" ]] || {
    echo "Missing GWAS input: $file" >&2
    exit 1
  }

  traits+=("$trait")
  files+=("$file")
  sample+=("$sample_prev")
  population+=("$population_prev")

done < "$MANIFEST"

# Confirm the expected number of phenotypes.
[[ "${#traits[@]}" -eq 14 ]] || {
  echo "Expected exactly 14 phenotypes for genetic correlations" >&2
  exit 1
}

# Estimate pairwise genetic correlations.
# Each command estimates correlations between the current
# phenotype and the subsequent phenotypes.

for ((i=0; i<${#traits[@]}-1; i++)); do

  input_files=$(IFS=,; echo "${files[*]:i}")
  input_sample=$(IFS=,; echo "${sample[*]:i}")
  input_pop=$(IFS=,; echo "${population[*]:i}")

  echo "Running LDSC genetic correlations: ${traits[i]}"

  "$PYTHON_BIN" "$LDSC_PY" \
    --rg "$input_files" \
    --samp-prev "$input_sample" \
    --pop-prev "$input_pop" \
    --ref-ld-chr "$EUR_LD_DIR/" \
    --w-ld-chr "$EUR_LD_DIR/" \
    --out "$ANALYSIS_OUT_DIR/rg/${traits[i]}_rg"

done

echo "Genetic correlation analyses completed."
echo "Results directory: $ANALYSIS_OUT_DIR/rg"
