import os
import subprocess
import time
from pathlib import Path, PurePosixPath
import packman
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

## CONFIG ##
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

alphafold_path="/homes/users/fpegenaute/opt/alphafold.F-2.1.0/alphafold"
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


## FUNCTIONS ##
def get_filename_ext(filepath):
    """
    Given a file path, split it by dots and get the extension as the last element
    Also, take the basename in the path and split it by the point, get the first
    element as the name of the file
    """

    extension = filepath.split(".")[-1]
    filename = PurePosixPath(filepath).stem
    
    return filename, extension

def choose_parser(pdbfile):
    """
    Check the extension of a file and return the according Bio.PDB structure 
    parser, and the file name (assumed identifier)
    """
    # get ID and extension
    identifier, extension= get_filename_ext(pdbfile)
    # identifier = identifier.upper()

    if extension == "pdb" or extension == "ent" :
        parser = PDBParser(QUIET=True)
    elif extension == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise NameError("""Your file must have \"pdb\", \"ent\" or \"cif\" as 
            an extension""")
    return parser, identifier

def number_of_jobs_in_queue(squeue="squeue"):
    """
    This functions returns the number of jobs in queue for a given
    user.
    """

    # Initialize #
    user_name = "$USER"

    process = subprocess.check_output([squeue, "-u", user_name])

    return process


def submit_command_to_queue(command, queue="normal", max_jobs_in_queue=None, queue_file=None, dummy_dir=".", submit="sbatch", squeue="squeue"):
    """
    This function submits any {command} to a cluster {queue}.
    @input:
    command {string}
    queue {string} by default it submits to any queue (partition)
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue
    """
    import hashlib

    if max_jobs_in_queue is not None:
        while number_of_jobs_in_queue(squeue) >= max_jobs_in_queue: time.sleep(5)

    cwd = os.path.join(dummy_dir)
    if not os.path.exists(cwd): 
        os.makedirs(cwd)

    script= os.path.join(cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")

    # Use the config for SLURM given in a file and add your command
    if queue_file is not None:
      fd=open(script,"w")

      with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
      fd.close()
      queue_standard.close()
      
      # choose which workload manager to use
      if queue is not None:
       if  submit=="qsub":
          os.system("%s -q %s %s" % (submit, queue,script))
       elif submit=="sbatch":
          os.system("%s -p %s %s" % (submit, queue,script))
       else:
          os.system("%s %s"% (submit,script))
      else:
        os.system("%s %s"% (submit,script))
    else:
      if queue is not None:
        os.system("echo \"%s\" | %s -q %s" % (command, submit, queue))
      else:
        os.system("echo \"%s\" | %s" % (submit,command))

def submit_AF_to_SLURM(query_fasta, outdir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None ):
    """
    This function submits AF2 command to a cluster via {workload_manager}.
    
    @input:
    query_fasta: fasta file to be analyzed
    outdir: output for *AlphaFold*, not SLURM
    workload_manager: SLURM by default
    dummy_dir=default directory for SLURM to look at, defaul: the current
    max_jobs_in_queue {int} limits the number of jobs in queue

    """
    import hashlib
    from datetime import date

    today = str(date.today())
    

    process = number_of_jobs_in_queue("squeue")
    print(process)
    
    cwd = os.path.join(dummy_dir)
    if not os.path.exists(cwd): 
        os.makedirs(cwd)

    script = os.path.join(cwd,"alpha_test.sh")

    command = "${alphafold_path}run_alphafold.sh -d"+AF2config["AF2datadir"]+" -o "+outdir+" -m "+AF2config["AF2preset_monomer"]+" -f "+query_fasta+" -t "+today+" --is_prokaryote_list="+AF2config["AF2_prokaryote"]+" -g "+AF2config["AF2_useGPU"]

    # Use the config for SLURM given in a file and add your command
    with open(script,"w") as batch_script:
        batch_script.write(SLURMconfig_AF)
        batch_script.write("%s\n\n"%(command))
        batch_script.write("conda deactivate\n")

    # Connect with the HPC Cluster:
    # command = "ssh -T marvin"
    # result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
    # Better like this: https://kb.iu.edu/d/aews, https://stackoverflow.com/questions/8382847/how-to-ssh-connect-through-python-paramiko-with-ppk-public-key
      
    os.system("%s %s" % (workload_manager, script))


def submit_RF_to_SLURM(query_fasta, outdir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None ):
    """
    This function submits AF2 command to a cluster via {workload_manager}.
    
    @input:
    query_fasta: fasta file to be analyzed
    outdir: output for *AlphaFold*, not SLURM
    workload_manager: SLURM by default
    dummy_dir=default directory for SLURM to look at, defaul: the current
    max_jobs_in_queue {int} limits the number of jobs in queue

    """
    from datetime import date

    today = str(date.today())
    

    process = number_of_jobs_in_queue("squeue")
    print(process)
    
    cwd = os.path.join(dummy_dir)
    if not os.path.exists(cwd): 
        os.makedirs(cwd)

    script = os.path.join(cwd,"alpha_test.sh")

    command = f"""run_pyrosetta_ver.sh {query_fasta} {outdir}"""

    # Use the config for SLURM given in a file and add your command
    with open(script,"w") as batch_script:
        batch_script.write(SLURMconfig_RF)
        batch_script.write("%s\n\n"%(command))
        batch_script.write("conda deactivate\n")

    os.system("%s %s" % (workload_manager, script))


def write_hng_file(pdbfile, hinges, outfile):
    """
    Generate a hinge file
    """
    filename = Path(pdbfile).stem
    ALL_RESIDUES = {}
    Protein = packman.molecule.load_structure(pdbfile)
    for i in Protein[0].get_chains():
        try:
            ALL_RESIDUES[i.get_id()] = sorted([i.get_id() for i in \
                            Protein[0][i.get_id()].get_residues() if i!=None])
        except:
            None

    select_count = 0
    last_hinge_end = 0
    fh = open(outfile, 'w')
    for numi, i in enumerate(hinges):
        current_hinge = hinges[numi]
        ChainOfHinge = current_hinge.get_elements()[0].get_parent().get_id()

        if select_count==0:
            hinge_res_ids = sorted([j.get_id() for j in \
                                                current_hinge.get_elements()])
            
            select_count += 1
            if(ALL_RESIDUES[ChainOfHinge][0]!=hinge_res_ids[0]):
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) \
                    +'\t'+ str(ALL_RESIDUES[ChainOfHinge][0])+':'\
                        +str(hinge_res_ids[0]-1)+'\n' )
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) \
                    +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            else:
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) \
                    +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            last_hinge_end = hinge_res_ids[-1]
        else:
            hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
            select_count += 1
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) \
                +'\t'+ str(last_hinge_end+1)+':'+str(hinge_res_ids[0]-1)+'\n' )
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) \
                +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            last_hinge_end = hinge_res_ids[-1]
        try:
            if(ChainOfHinge != hinges[numi+1].get_elements()[0].get_parent().get_id()):
                select_count += 1
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) \
                    +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
                last_hinge_end = 0
        except:
            None

        if(last_hinge_end != ALL_RESIDUES[ChainOfHinge][-1]):
            select_count += 1
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) \
                +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
    fh.flush()
    fh.close()

def pdb_to_fasta(pdbfile, outdir):
    """
    Given a PDB/mmCif file, make a fasta file in the outdir
    """
    # get ID and extension
    parser, identifier = choose_parser(pdbfile)
    outfile = os.path.join(outdir, f"{identifier}_covered.fasta")
    structure = parser.get_structure(identifier, pdbfile) 
    model = structure[0]
    with open(outfile, "w") as out_fasta:
        for chain in model.get_chains():
            out_fasta.write('>' + identifier + "\n")
            for res in chain.get_residues(): 
                if res.id[0] == " ": 
                    # seq1 changes 3 letter symbls to 1
                    residue = seq1(str(res.get_resname()))              
                    out_fasta.write(residue)
    
    return outfile
    
def get_chain_names(structure_file):
    """
    Given a PDB/mmCif file, return a list with the names 
    of the chains
    """

    # get ID and extension
    parser, identifier = choose_parser(structure_file)

    # Get the structure
    structure = parser.get_structure(identifier, structure_file) 
    model = structure[0]
    chains = []
    
    for chain in model.get_chains():
        chains.append(chain.get_id())
   
    return chains

def get_residue_range(structure_file, chain=None):
    """
    Given a PDB/mmCif file, extract the first and last aa positions.
    Return them as a tuple

    Optionally, specify a chain for doing this operation
    """

    parser, identifier = choose_parser(structure_file)

    # Get the structure
    structure = parser.get_structure(identifier, structure_file) 
    model = structure[0]
    chains = []
    
    if chain is not None:
        chain = model[str(chain)]
        i = 0
        for residue in chain.get_residues():
            if residue.get_full_id()[3][0] == " ": # Exclude hetatm and h20
                if i == 0:
                    first = residue.get_full_id()[3][1]
                    i += 1
                else:
                    last = residue.get_full_id()[3][1]
        
    else:
        for chain in model.get_chains():
            i = 0
            for residue in chain.get_residues():
                if residue.get_full_id()[3][0] == " ": # Exclude hetatm and h20
                    if i == 0:
                        first = residue.get_full_id()[3][1]
                        i += 1
                    else:
                        last = residue.get_full_id()[3][1]
   
    return (first, last)




## CLASSES
    
    





if __name__ == "__main__":
    pass
