from bin.process_predicted_model import *

from Bio import *
from Bio.PDB import MMCIFParser

structure = MMCIFParser().get_structure('1a7f', '1a7f.cif')    


model = structure[0]
chain = model['A']

for i in chain.get_residues():
    print(f"{i.get_full_id()[2]}\t{i.get_full_id()[3][1]}" )