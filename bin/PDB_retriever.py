from Bio.PDB import *
from src.dicts import aa_symbol_to_letter

# Retrieve the templates used by alphafold

def retrieve_templates(pdblist):
    pdbl = PDBList()
    for template in pdblist:
        pdbl.retrieve_pdb_file(template)
    



# Just an example input pdb
record = '1pa2.pdb'

# run parser
parser = PDBParser(QUIET=True)
structure = parser.get_structure('struct', record)    

# iterate each model, chain, and residue
# printing out the sequence for each chain

for model in structure:
    for chain in model:
        seq = []
        for residue in chain:
            seq.append(aa_symbol_to_letter[residue.resname])
        print('>some_header\n',''.join(seq))