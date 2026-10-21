#!/usr/bin/env bash
# Reconstructed from archived August 2026 standalone LDSC v1.0.1 execution logs.
# This is NOT a verbatim recovered original submission script.
# Requires an LDSC Python environment, European reference LD scores and the
# previously munged, HapMap3-aligned GWAS summary statistics.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MANIFEST="$ROOT/config/trait_prevalence_specification.tsv"

: "${LDSC_PY:?Set LDSC_PY to the full path of ldsc.py}"
: "${MUNGED_DIR:?Set MUNGED_DIR to the directory containing .hm3.sumstats.gz files}"
: "${EUR_LD_DIR:?Set EUR_LD_DIR to the European eur_w_ld_chr directory}"
: "${ANALYSIS_OUT_DIR:?Set ANALYSIS_OUT_DIR to a private working output directory}"
PYTHON_BIN="${PYTHON_BIN:-python}"

[[ -f "$LDSC_PY" ]] || { echo "Missing LDSC script: $LDSC_PY" >&2; exit 1; }
[[ -d "$EUR_LD_DIR" ]] || { echo "Missing European LD-score directory: $EUR_LD_DIR" >&2; exit 1; }
mkdir -p "$ANALYSIS_OUT_DIR/h2"

# The manifest covers the seven epilepsy subtypes and seven brain volumes.
# Broad GGE and FE are not included because their original --h2 logs and the
# corresponding input-file specifications have not been recovered.
while IFS=$'\t' read -r trait filename sample_prev population_prev; do
  [[ "$trait" == "Trait" || -z "$trait" ]] && continue
  file="$MUNGED_DIR/$filename"
  [[ -f "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }

  args=("$PYTHON_BIN" "$LDSC_PY" --h2 "$file"
        --ref-ld-chr "$EUR_LD_DIR/" --w-ld-chr "$EUR_LD_DIR/"
        --out "$ANALYSIS_OUT_DIR/h2/$trait")
  if [[ "$sample_prev" != "nan" && "$sample_prev" != "NA" ]]; then
    args+=(--samp-prev "$sample_prev" --pop-prev "$population_prev")
  fi

  printf 'Running h2 for %s\n' "$trait"
  "${args[@]}"
done < "$MANIFEST"
echo "Completed. See $ANALYSIS_OUT_DIR/h2"
