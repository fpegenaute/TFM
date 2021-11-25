###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

import argparse
from bin.blast import *
from bin.PDB_retriever import *
import  bin.config as cfg
import logging as l
from pathlib import Path


parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="[options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")

parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")
parser.add_argument("-p", "--pdbdir", 
                    help="Directory to store the retrieved PDBs", 
                    default=".")
parser.add_argument("-f", "--fastadir", 
                    help="Directory to store the retrieved FASTA files", 
                    default=".")

args = parser.parse_args()




### Initializing the LOG system ###

fasta = args.FASTA
query_name = Path(fasta).stem
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
pdb_dir = f"{args.pdbdir}/{query_name}"
fasta_dir = f"{args.fastadir}/{query_name}"

Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)

# Retrieve
if exact_matches:
    req = retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
    # Check lengths of the actual PDBs






if not exact_matches:
    l.info("No templates were found in the PDB")


