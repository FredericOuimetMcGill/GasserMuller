#!/bin/sh

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=63
#SBATCH --mem=240G
#SBATCH -o log/%x_%j.out
#SBATCH --mail-user=frederic.ouimet.23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=message/%x_%j-result.txt
#SBATCH --error=message/%x_%j-error.txt
#SBATCH --exclusive

module load StdEnv/2023 r/4.3.1
Rscript GM_vs_LL_vs_NW.R
