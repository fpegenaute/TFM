###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

import argparse
from bin.blast import *
from bin.PDB_retriever import *
import logging as l


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
blastdb = "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"
l.info(f"BLAST database is located at: {blastdb}")

# Run BLAST
l.info(f"Starting BLAST. Query: {fasta}, Database: {Path(blastdb).stem}")
outblast = run_blast_local(fasta, blastdb)
l.info(f"BLAST results stored in : {outblast}")

# Catch exact matches
exact_matches = exact_match_retriever(outblast)
l.info(f" The target sequence is already in the PDB with code/s: {exact_matches.keys()}")


# Retrieve exact matches from the PDB
pdb_dir = args.pdbdir
fasta_dir = args.fastadir

if exact_matches:
    retrieve_pdb_info(exact_matches.keys(), pdb_dir, fasta_dir)
if not exact_matches:
    l.info("Your full sequence is not on the PDB")


