## Run the program in a HPC cluster Managed by Slurm
SLURMconfig = """#!/bin/bash

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 200G
#SBATCH -t 5-20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ferran.pegenaute@upf.edu

module purge
module load modulepath/noarch
module spider BLAST+/2.4.0
module load AlphaFold/2.1.0-Miniconda3-4.7.10
source activate alphafold-2.1.0
module load CUDA/11.1.0
"""