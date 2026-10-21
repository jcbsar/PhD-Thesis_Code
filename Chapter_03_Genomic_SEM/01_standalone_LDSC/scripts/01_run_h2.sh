
#!/usr/bin/env bash
# Standalone LDSC v1.0.1: SNP-based heritability analyses.
# Command parameters follow the recorded August 2026 LDSC analyses.
# Requires munged, HapMap3-aligned GWAS summary statistics and
# European-ancestry reference LD scores.

set -euo pipefail

# Locate the phenotype configuration file.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MANIFEST="$ROOT/config/trait_prevalence_specification.tsv"

# Define software, input and output locations.
: "${LDSC_PY:?Set LDSC_PY to the full path of ldsc.py}"
: "${MUNGED_DIR:?Set MUNGED_DIR to the directory containing .hm3.sumstats.gz files}"
: "${EUR_LD_DIR:?Set EUR_LD_DIR to the European eur_w_ld_chr directory}"
: "${ANALYSIS_OUT_DIR:?Set ANALYSIS_OUT_DIR to the output directory}"

PYTHON_BIN="${PYTHON_BIN:-python}"

# Check that required files and directories exist.
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

mkdir -p "$ANALYSIS_OUT_DIR/h2"

# Estimate heritability for seven epilepsy subtypes and seven
# subcortical brain-volume phenotypes.
#
# Broad GGE and FE are included in Figure 3.1A but are not
# represented in this execution manifest.

while IFS=$'\t' read -r trait filename sample_prev population_prev; do

  [[ "$trait" == "Trait" || -z "$trait" ]] && continue

  file="$MUNGED_DIR/$filename"

  [[ -f "$file" ]] || {
    echo "Missing GWAS input: $file" >&2
    exit 1
  }

  # Set LDSC heritability arguments.
  args=(
    "$PYTHON_BIN" "$LDSC_PY"
    --h2 "$file"
    --ref-ld-chr "$EUR_LD_DIR/"
    --w-ld-chr "$EUR_LD_DIR/"
    --out "$ANALYSIS_OUT_DIR/h2/$trait"
  )

  # Apply liability-scale transformation for epilepsy phenotypes.
  # Sample prevalence is 0.50 when using effective sample sizes.
  # Quantitative brain-volume traits do not require this correction.
  if [[ "$sample_prev" != "nan" && "$sample_prev" != "NA" ]]; then
    args+=(
      --samp-prev "$sample_prev"
      --pop-prev "$population_prev"
    )
  fi

  echo "Running LDSC heritability analysis: $trait"
  "${args[@]}"

done < "$MANIFEST"

echo "Heritability analyses completed."
echo "Results directory: $ANALYSIS_OUT_DIR/h2"

echo "Completed. See $ANALYSIS_OUT_DIR/h2"
