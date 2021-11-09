###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

import argparse
from bin.blast import *
from bin.PDB_retriever import *

parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="[options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")

parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")

args = parser.parse_args()


if args.verbose:
    print("verbosity turned on")

## Check if the input sequence is already in the PDB ##

# Set vars for BLAST
blastdb = "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"

fasta= "test.fa"

# Run BLAST
outblast =run_blast_local(fasta, blastdb)
print(outblast)