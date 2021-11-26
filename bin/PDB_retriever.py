from posixpath import join
from Bio import PDB
from Bio.PDB import *
from iotbx import pdb
import requests
import os
import logging as l
from utilities import get_filename_ext
  
def retrieve_pdb_info(hit_dict, pdb_dir, fasta_dir):
    """
    Given a list of PDB codes, it downloads the PDB/MMcif files, and gets the data
    from the pdb entry, from where it extracts the sequence and writes it into a fasta file

    Inputs:

    - hit_dict: List of PDB codes
    - pdb_dir: Directory to store the PDB files
    - fasta_dir: Directory to store the Fasta files

    It returns the request
    """

    # Initialize object PDBList
    pdbl = PDBList() 
    pdb_list = hit_dict.keys()

    input(f"PDB List = {len(pdb_list)} Press Enter to continue...")
    for template in pdb_list:
        # Dowload the PDB file from the web
        pdbl.retrieve_pdb_file(template, pdir=pdb_dir)     
        req = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{template}').json()[template.lower()]
        fasta_template_name = f"{fasta_dir}/{template}.fa"
        with open(fasta_template_name, "w") as f:
            # Print chain (atm A)
            header = [template,req[0]["in_chains"][0]]
            headerstr = "_".join(header)
            sequence = req[0]['sequence']
            f.write(">"+headerstr+"\n")
            f.write(sequence)

    return




def check_PDB_len(pdbfile, chain):
    """
    Given a PDB/mmCIF file, check the length of the desired chain
    """

    # get ID and extension

    identifier, extension= get_filename_ext(pdbfile)
    identifier = identifier.upper()
    l.info(f"ID: {identifier} . Extension: {extension}, ")

    if extension == "pdb":
        parser = PDB.PDBParser(QUIET=True)
    elif extension == "cif":
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        raise NameError("Your file has to have \"pdb\" or \"cif\" as an extension")

    # Get the structure
    structure = parser.get_structure(identifier, pdbfile) 
    model = structure[0]
    
    res_number = 0
    non_resi = 0

    for r in model[chain].get_residues():
        if r.id[0] == ' ':
            res_number +=1
        else:
            non_resi +=1

    return res_number

if __name__ == "__main__":
    
    # Example
    pdb_list = ["1pa2","5tnt"]
    pdb_dir = "./templates/PDB/"
    fasta_dir="./templates/FASTA/"
    req = retrieve_pdb_info(pdb_list, pdb_dir, fasta_dir)


