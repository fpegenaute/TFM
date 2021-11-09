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

args = parser.parse_args()




### Initializing the LOG system ###

fasta= args.FASTA
l.basicConfig(format = "%(levelname)s:%(message)s", 
                        filename = f"{Path(fasta).stem}.log", level = l.DEBUG)
l.debug("...STARTING...\n")		

# If verbose is set, the LOG file is also printed in STDOUT
if args.verbose:		
	l.getLogger().addHandler(l.StreamHandler())		


## 1. Check if the input sequence is already in the PDB  

# Locate the Database
blastdb = "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"
l.info(f"BLAST database is located at: {blastdb}")



# Run BLAST
outblast = run_blast_local(fasta, blastdb)
print(outblast)
print(fasta)