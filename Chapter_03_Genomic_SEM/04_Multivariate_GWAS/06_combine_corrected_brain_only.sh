
#!/bin/bash
#SBATCH --job-name=Brain_combine_corrected
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/jcbsar002/gsem/logs/Brain_combine_corrected_%j.out
#SBATCH --error=/scratch/jcbsar002/gsem/logs/Brain_combine_corrected_%j.err

set -euo pipefail

export R_LIBS_USER="$HOME/R/4.6"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /scratch/jcbsar002/gsem/scripts

OUT="/scratch/jcbsar002/gsem/final_streamed_corrected/GGE_brain_Brain_userGWAS_full_streamed.csv"

if [[ -e "$OUT" || -e "${OUT}.gz" ]]; then
    echo "Brain output already exists. Refusing to overwrite."
    exit 1
fi

Rscript 06_stream_combine_one_path_corrected.R GGE 2 Brain

echo "Brain-factor combination completed."
echo "Compression intentionally skipped."
