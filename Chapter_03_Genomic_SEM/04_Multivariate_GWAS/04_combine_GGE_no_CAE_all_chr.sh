
#!/bin/bash
#SBATCH --job-name=combine_GGE_no_CAE
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/jcbsar002/gsem/logs/combine_GGE_no_CAE_all_%j.out
#SBATCH --error=/scratch/jcbsar002/gsem/logs/combine_GGE_no_CAE_all_%j.err

set -euo pipefail

export R_LIBS_USER="$HOME/R/4.6"

cd /scratch/jcbsar002/gsem/scripts

Rscript 04_combine_GGE_no_CAE_all_chr.R
