
blastconfig = {
    "blastdb" : "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"

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
SLURMconfig = """
#!/bin/bash

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
source activate alphafold-2.1.0
module load CUDA/11.1.0
"""
SSHconfig = {
    "HPCCluster" : "marvin",
    "Public_key": "", 
}










if __name__ == "__main__":
    print("config.py file not executable!")