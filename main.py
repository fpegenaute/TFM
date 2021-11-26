###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

import Bio.SeqIO as IO
import argparse
from bin.blast import *
from bin.PDB_retriever import *
import  bin.config as cfg
import logging as l
from pathlib import Path
import os
import shutil
from bin.extract_flexible_residues import extract_residue_list
from bin.slurm_utilities import write_batch_script


parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="main.py input_fasta outdir AF_preset [options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")

parser.add_argument("outdir", 
                    help="Output directory to store the retrieved PDBs", 
                    default=".")
parser.add_argument("model_preset", 
                    help="model preset for AlphaFold2", 
                    default="monomer")
parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")
parser.add_argument("-c", "--custom", 
                    help="Use custom slurm batch script for AlphaFold2")
parser.add_argument("-n", "--noslurm",
                    help="""run AlphaFold2 locally, without using SLURM.
                    Only recommended if this python script is submitted on a 
                    batch script on its own""", 
                    default=False, 
                    action="store_true")




args = parser.parse_args()

## Info about the query
fasta = args.FASTA
query_name = Path(fasta).stem

record_dict = IO.to_dict(IO.parse(fasta, "fasta"))
if len(record_dict.keys()) > 1:
    raise NotImplemented("This program only works with one sequence at a time")
if len(record_dict.keys()) == 1:
    if list(record_dict.keys())[0] != query_name:
        raise NameError(f"""Please, make sure your filename and fasta identifier 
        coincide. filename: {query_name} / ID name: {record_dict.keys()}""")
    
    query_length = len(record_dict[query_name].seq)
    print(f"{query_length}")


# For now, it only forks with one squence at a time
# for key in record_dict.items():
#     print(key[0],"\n ",len(key[1].seq))






### Initializing the LOG system ###


l.basicConfig(format = "%(levelname)s:%(message)s", 
                        filename = f"{query_name}.log", level = l.DEBUG)
l.debug("...STARTING...\n")		

# If verbose is set, the LOG file is also printed in STDOUT
if args.verbose:		
	l.getLogger().addHandler(l.StreamHandler())		


## 1. Check if the input sequence is already in the PDB  

# Locate the Database
blastdb = cfg.blastconfig["blastdb"]
l.info(f"BLAST database is located at: {blastdb}")

# Run BLAST
l.info(f"Starting BLAST. Query: {fasta}, Database: {Path(blastdb).stem}")
outblast = run_blast_local(fasta, blastdb)
l.info(f"BLAST results stored in : {outblast}")

# Catch exact matches
exact_matches = exact_match_retriever(outblast)
l.info(f" The target sequence has close homologs in the PDB with code/s: {exact_matches.keys()}")


# Retrieve exact matches from the PDB

# Create folders
pdb_dir = f"{args.outdir}/PDB/{query_name}"
fasta_dir = f"{args.outdir}/FASTA/{query_name}"

Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)

# Retrieve
if exact_matches:
    retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
    # Check lengths of the actual PDBs and store them accordingly
    for file in os.listdir(pdb_dir):
        l.info(f"File: {file}")
        current = os.path.join(pdb_dir, file)
        if os.path.isfile(current):
            identifier = file.split(".")[0].upper()
            pdb_len = check_PDB_len(current, exact_matches[identifier])
            l.info(f"Length of the template {identifier}: {pdb_len}")
            # Store partial matches (<90% of the query length)
            if pdb_len < 10 and pdb_len < (0.9*query_length):
                print(f"{identifier} has length {pdb_len}, it will be stored as a partial match")
                try:
                    shutil.move(file, f"./partial/{file}")
                except Exception:
                    l.info(f"./partial/{file} does not exist, it will be created")
                    os.mkdir("./partial/")
                    shutil.move(file, f"./partial/{file}")

### Submit a sob in Slurm with the AlphaFold run

# Make folder for the AF2 output
af_dir = f"{args.outdir}/ALPHAFOLD/{query_name}"
Path(af_dir).mkdir(parents=True, exist_ok=True)

if args.custom:
    slurm_script = args.custom
else:
    slurm_dir = f"{af_dir}"
    Path(slurm_dir).mkdir(parents=True, exist_ok=True)
    write_batch_script(slurm_dir, model_preset)





# Extract confident regions

domains_dir = f"{af_dir}/DOMAINS/"
Path(domains_dir).mkdir(parents=True, exist_ok=True)



