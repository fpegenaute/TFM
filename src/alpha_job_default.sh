#!/bin/bash

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 200G
#SBATCH -t 5-20:00:00
#SBATCH --output=reports/alpha_ptm_monomer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ferran.pegenaute@upf.edu

module purge
module load modulepath/noarch
module load AlphaFold/2.1.0-Miniconda3-4.7.10
source activate alphafold-2.1.0
module load CUDA/11.1.0

# bash $alphafold_path/run_alphafold.sh -d $ALPHAFOLD_DATA_DIR -o $outdir -m monomer_ptm -f $filename -t 2021-11-23 --is_prokaryote_list=false -g true
