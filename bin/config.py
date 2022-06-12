# This file contains the location of the database for BLAST to use, and the
# text to be placed on the batch scripts for SLURM when running AlphaFold or 
# RoseTTaFold

blastconfig = {
    "blastdb" : "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa",
    "eval_cutoff" : 0.000005, 
    "best_hit_score_edge" : 0.1, 
    "best_hit_overhang" : 0.25,
    "min_length_match" : 5
}

AF2config = {
    "AF2command" : "bash $alphafold_path/run_alphafold.sh", 
    "AF2datadir" : "$ALPHAFOLD_DATA_DIR", 
    "AF2preset_monomer" : "monomer",
    "AF2preset_monomer_ptm" : "monomer_ptm", 
    "AF2preset_multimer" : "multimer", 
    "AF2_useGPU" : "true", 
    "AF2_prokaryote" : "false", 
    }

SLURMconfig_AF = """#!/bin/bash

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
module load AlphaFold/2.1.0-Miniconda3-4.7.10
source activate alphafold.F-2.1.0
module load CUDA/11.1.0

module load CUDA/11.1.0

alphafold_path="/homes/users/fpegenaute/opt/alphafold.F-2.1.0/alphafold/"
"""

SLURMconfig_RF = """#!/bin/bash

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
module load RoseTTAFold/v1.1.0-Miniconda3-4.7.10
"""

PACKMANconfig ={
    "alpha" : 4.5,
    "alpha_range" : None
}

CCTBXconfig = {
    "AF2_maximum_rmsd" : 1.5,
    "AF2_domain_size" : 15,
    "RF_maximum_rmsd" : 2, 
    "RF_domain_size" : 15
}
if __name__ == "__main__":
    print("config.py file not executable!")