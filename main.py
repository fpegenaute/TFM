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
parser.add_argument("-m", "--model_preset", 
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
fasta = args.FASTA
query_name = Path(fasta).stem.split('.')[0]

### Initializing the LOG system ###

logdir = os.path.join(args.outdir, query_name, "LOG", "")
Path(logdir).mkdir(parents=True, exist_ok=True)

l.basicConfig(format = "%(levelname)s:%(message)s", 
                        filename = os.path.join(logdir, f"{query_name}.log"), 
                        level = l.DEBUG)

l.debug("...STARTING...\n")		

# If verbose is set, the LOG file is also printed in STDOUT
if args.verbose:		
	l.getLogger().addHandler(l.StreamHandler())		


## Info about the query

record_dict = IO.to_dict(IO.parse(fasta, "fasta"))
if len(record_dict.keys()) > 1:
    raise NotImplemented("This program only works with one sequence at a time")
if len(record_dict.keys()) == 1:
    if list(record_dict.keys())[0] != query_name:
        raise NameError(f"""Please, make sure your filename and fasta identifier 
        coincide. filename: {query_name} / ID name: {record_dict.keys()}""")
    
    query_length = len(record_dict[query_name].seq)
    l.info(f"Query length: {query_length}")


# For now, it only forks with one squence at a time
# for key in record_dict.items():
#     print(key[0],"\n ",len(key[1].seq))

# Create folders. last path element is empty to add a slash
pdb_dir = os.path.join(args.outdir, query_name, "PDB", "" )
fasta_dir = os.path.join(args.outdir, query_name, "FASTA", "" )

Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)

## 1. Check if the input sequence is already in the PDB  

# Locate the Database
blastdb = cfg.blastconfig["blastdb"]
l.info(f"BLAST database is located at: {blastdb}")

# Create folders. last path element is empty to add a slash
blast_dir = os.path.join(args.outdir, query_name,  "BLAST", "" )
Path(blast_dir).mkdir(parents=True, exist_ok=True)
l.info(f"The BLAST output will be stored in:{blast_dir}")

# Run BLAST
l.info(f"Starting BLAST. Query: {fasta}, Database: {Path(blastdb).stem}")
outblast = run_blast_local(fasta, blastdb, blast_dir)
l.info(f"BLAST results stored in : {outblast}")

# Catch exact matches
exact_matches = exact_match_retriever(outblast)
l.info(f""" The target sequence has close homologs in the PDB with 
    code/s: {exact_matches.keys()}""")


# Retrieve exact matches from the PDB

# Retrieve
if exact_matches:
    retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
    # Check lengths of the actual PDBs and store them accordingly
    for file in os.listdir(pdb_dir):
        l.info(f"File: {file}")
        current = os.path.join(pdb_dir, file)
        if os.path.isfile(current):
            identifier = file.split(".")[0].upper()
            # Check which chain of the hit should it choose
            
            # Make the directory for the chains
            chain_dir = os.path.join(pdb_dir, "CHAINS", "" )
            Path(chain_dir).mkdir(parents=True, exist_ok=True)
            # Extract the desired chain
            splitter = ChainSplitter(mmcif=True, out_dir=os.path.join(pdb_dir, "CHAINS"))
            splitter.make_pdb(os.path.join(pdb_dir, file), exact_matches[identifier[:4]] )


            pdb_len = check_PDB_len(current, exact_matches[identifier[:4]])
            l.info(f"Length of the template {identifier[:4]}: {pdb_len}")
            # Store partial matches (<95% of the query length)
            l.info(f"PDB_LEN: {pdb_len} . QUERY_LEN: {0.95*query_length}")
            if pdb_len > 10 and pdb_len < (query_length):
                l.info(f"""{identifier[:4]} has length {pdb_len}, it will be stored 
                    as a partial match""")
                try:
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"partial", file))
                except Exception:
                    directory = os.path.join(pdb_dir,"partial",file)
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    os.mkdir(os.path.join(pdb_dir,"partial"))
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"partial", file))
            if pdb_len < 10 and pdb_len > (0.95*query_length):
                l.info(f"""{identifier[:4]} has length {pdb_len}, it will be stored 
                    as a full-length match""")
                try:
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"total", file))
                except Exception:
                    directory = os.path.join(pdb_dir,"total",file)
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    os.mkdir(os.path.join(pdb_dir,"partial"))
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"total", file))




### Submit a sob in Slurm with the AlphaFold run

print("Alphafold will not run, this is a test")


af_dir = os.path.join(args.outdir,  query_name, "ALPHAFOLD", "" )
Path(af_dir).mkdir(parents=True, exist_ok=True)
# Make folder for the AF2 output
l.info(f"Creating folder for AF2 output in:{af_dir}")

if args.custom:
    slurm_script = args.custom
else:
    slurm_dir = f"{af_dir}"
    Path(slurm_dir).mkdir(parents=True, exist_ok=True)





# Extract confident regions

l.info("Extracting high confidence domains")
domains_dir = os.path.join(af_dir, "DOMAINS", "")
Path(domains_dir).mkdir(parents=True, exist_ok=True)
l.info(f"Domains will be stored in:{domains_dir}")
for filename in os.listdir(af_dir):
    if os.path.isfile(join(af_dir, filename)):
        l.info(f"Processing file: {filename}")
        conf_domains = extract_residue_list(os.path.join(af_dir, filename), 
        domains_dir)
        l.info(f"Residue list of confident domains: {conf_domains}")



