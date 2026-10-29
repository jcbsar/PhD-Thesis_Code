
#!/bin/bash
#SBATCH --job-name=GGE_combine_corrected
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --array=1-2%1
#SBATCH --output=/scratch/jcbsar002/gsem/logs/GGE_combine_corrected_%A_%a.out
#SBATCH --error=/scratch/jcbsar002/gsem/logs/GGE_combine_corrected_%A_%a.err

set -euo pipefail

export R_LIBS_USER="$HOME/R/4.6"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

cd /scratch/jcbsar002/gsem/scripts

# Only combine the corrected GGE model outputs
case "${SLURM_ARRAY_TASK_ID}" in
  1)
    ANALYSIS="GGE"
    PATH_INDEX="1"
    PATH_LABEL="GGE"
    ;;
  2)
    ANALYSIS="GGE"
    PATH_INDEX="2"
    PATH_LABEL="Brain"
    ;;
  *)
    echo "Invalid array task ID."
    exit 1
    ;;
esac

OUTPUT_DIR="/scratch/jcbsar002/gsem/final_streamed_corrected"
OUTPUT_FILE="${OUTPUT_DIR}/${ANALYSIS}_brain_${PATH_LABEL}_userGWAS_full_streamed.csv"

# Protect existing corrected results from accidental overwrite
if [[ -e "${OUTPUT_FILE}" || -e "${OUTPUT_FILE}.gz" ]]; then
  echo "Output already exists. Refusing to overwrite."
  exit 1
fi

echo "========================================"
echo "Corrected GGE GWAS combination"
echo "Factor: ${PATH_LABEL}"
echo "========================================"

Rscript 06_stream_combine_one_path_corrected.R \
  "${ANALYSIS}" \
  "${PATH_INDEX}" \
  "${PATH_LABEL}"

# Verify output before compression
if [[ ! -s "${OUTPUT_FILE}" ]]; then
  echo "Expected combined CSV is missing or empty."
  exit 1
fi

gzip -n "${OUTPUT_FILE}"

echo "Combination completed successfully:"
echo "${OUTPUT_FILE}.gz"
