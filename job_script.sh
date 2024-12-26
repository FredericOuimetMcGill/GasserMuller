#!/bin/sh

#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=39
#SBATCH --mem=180G
#SBATCH -o log/%x_%j.out
#SBATCH --mail-user=frederic.ouimet.23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=message/%x_%j-result.txt
#SBATCH --error=message/%x_%j-error.txt
#SBATCH --exclusive
#SBATCH --account=def-cgenest

module load StdEnv/2023 r/4.3.1
Rscript LL_vs_NW_random_design.R
