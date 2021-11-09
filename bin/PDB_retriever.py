from posixpath import join
from Bio import PDB
from Bio.PDB import *
import requests

  
def retrieve_pdb_info(pdb_list, pdb_dir, fasta_dir):
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

    for template in pdb_list:
        # Dowload the PDB file from the web
        pdbl.retrieve_pdb_file(template, pdir=pdb_dir)
        req = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{template}').json()[template.lower()]
        with open(fasta_dir+template+".fa", "w") as f:
            # Print chain (usually A)
            header = [template,req[0]["in_chains"][0]]
            headerstr = "_".join(header)
            sequence = req[0]['sequence']
            f.write(">"+headerstr+"\n")
            f.write(sequence)
    return

if __name__ == "__main__":
    
    # Example
    pdb_list = ["1pa2","5tnt"]
    pdb_dir = "./templates/PDB/"
    fasta_dir="./templates/FASTA/"
    retrieve_pdb_info(pdb_list, pdb_dir, fasta_dir)


