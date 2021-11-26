import subprocess

multimer_preset_str= """#!/bin/bash\n
    #SBATCH -N 1
    #SBATCH -n 2
    #SBATCH -p normal
    #SBATCH --gres=gpu:1
    #SBATCH --gres-flags=enforce-binding
    #SBATCH --mem 200G
    #SBATCH -t 5-20:00:00
    #SBATCH --output=reports/alpha_SC1_%j.out
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=ferran.pegenaute@upf.edu


    filename="SC1_exocyst.fa"
    outdir="./output_multimer/"

    module purge
    module load modulepath/noarch
    module load AlphaFold/2.1.0-Miniconda3-4.7.10
    source activate alphafold-2.1.0
    module load CUDA/11.1.0

    bash $alphafold_path/run_alphafold.sh -d $ALPHAFOLD_DATA_DIR -o $outdir -m multimer -f $filename -t $

    conda deactivate

"""

monomer_preset_str= """#!/bin/bash\n
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 200G
#SBATCH -t 5-20:00:00
#SBATCH --output=reports/alpha_monomer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ferran.pegenaute@upf.edu


filename="SC1_exocyst.fa"
outdir="./output_multimer/"

module purge
module load modulepath/noarch
module load AlphaFold/2.1.0-Miniconda3-4.7.10
source activate alphafold-2.1.0
module load CUDA/11.1.0

bash $alphafold_path/run_alphafold.sh -d $ALPHAFOLD_DATA_DIR -o $outdir -m multimer -f $filename -t $

conda deactivate
"""

test_slurm = """#!/bin/bash\n
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p normal
#SBATCH --mem 4G
#SBATCH -t 5:00:00
#SBATCH --output=reports/test_slurm_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ferran.pegenaute@upf.edu


filename="SC1_exocyst.fa"
outdir="./output_multimer/"

module purge
module load modulepath/noarch
module load AlphaFold/2.1.0-Miniconda3-4.7.10
source activate alphafold-2.1.0
module load CUDA/11.1.0 

bash $alphafold_path/run_alphafold.sh -d $ALPHAFOLD_DATA_DIR -o $outdir -m multimer -f $filename -t $

conda deactivate
"""

af_presets = {
    "monomer":monomer_preset_str,
    "multimer": multimer_preset_str,
}

def write_batch_script(slurm_dir, model_preset_str):
    """
    Given a preset name for AlphaFold2, generate a batch script

    -usage write_batch_script(slurm_dir, af_presets["preset"])
    """

    with open(f"{slurm_dir}/testslurm.txt", "w") as slurmfile:
        slurmfile.write(model_preset_str)
    
    return



if __name__ == "__main__":
    write_batch_script(test_slurm)