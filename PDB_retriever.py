from posixpath import join
from Bio import PDB
from Bio.PDB import *
import requests

  
template_dir = "./templates/PDB/"
fasta_dir="./templates/FASTA/"
  
# Example
pdblist = ["1pa2","5tnt"]

def retrieve_pdb_info(pdblist):
    """
    Given a list of PDB codes, it downloads the PDB/MMcif files, and gets the data
    from the pdb entry, from where it extracts the sequence and writes it into a fasta file

    Inputs:

    - pdb_list: List of PDB codes
    - pdb_dir: Directory to store the PDB files
    - fasta_dir: Directory to store the Fasta files
    """
    # Initialize object PDBList
    pdbl = PDBList() 

    for template in pdblist:
        # Dowload the PDB file from the web
        pdbl.retrieve_pdb_file(template, pdir=template_dir)
        req = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{template}').json()[template.lower()]
        with open(fasta_dir+template+".fa", "w")as f:
            # Print chain A
            header = [template,req[0]["in_chains"][0]]
            headerstr = "_".join(header)
            sequence = req[0]['sequence']
            f.write(headerstr+"\n")
            f.write(sequence)
            
    
        
retrieve_pdb_info(pdblist)




