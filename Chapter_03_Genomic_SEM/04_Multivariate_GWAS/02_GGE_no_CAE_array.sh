
#!/bin/bash
#SBATCH --job-name=GGE_no_CAE
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-22%4
#SBATCH --output=/scratch/jcbsar002/gsem/logs/GGE_no_CAE_chr%a_%A.out
#SBATCH --error=/scratch/jcbsar002/gsem/logs/GGE_no_CAE_chr%a_%A.err

set -euo pipefail

export R_LIBS_USER="$HOME/R/4.6"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export TMPDIR="/scratch/jcbsar002/gsem/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT

cd /scratch/jcbsar002/gsem/scripts

Rscript 02_run_userGWAS_chr_no_CAE.R \
  "${SLURM_ARRAY_TASK_ID}" \
  "${SLURM_CPUS_PER_TASK}"
