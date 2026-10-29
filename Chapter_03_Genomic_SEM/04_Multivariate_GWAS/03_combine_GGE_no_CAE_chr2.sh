
#!/bin/bash
#SBATCH --job-name=combine_chr2
#SBATCH --account=psychiatry
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/jcbsar002/gsem/logs/combine_GGE_no_CAE_chr2_%j.out
#SBATCH --error=/scratch/jcbsar002/gsem/logs/combine_GGE_no_CAE_chr2_%j.err

set -euo pipefail

export R_LIBS_USER="$HOME/R/4.6"

cd /scratch/jcbsar002/gsem/scripts

Rscript 03_combine_GGE_no_CAE_chr2.R
