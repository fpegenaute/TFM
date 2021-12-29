
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
import logging as l
from Bio import SeqIO
import os
from bin.utilities import get_filename_ext
import mplcursors


# def PDB_get_resid_set(structure_file):
#     """
#     Given a PDB or MMCif File, return a set of [Res_ID] ONE CHAIN
#     """
#     identifier, extension = get_filename_ext(structure_file)

#     if extension == "pdb":
#         parser = PDBParser(QUIET=True)
#     elif extension == "cif":
#         parser = MMCIFParser(QUIET=True)
#     else:
#         raise NameError("Your file has to have \"pdb\" or \"cif\" as an extension")

#     structure = parser.get_structure(identifier, structure_file)    

#     model = structure[0]
#     resid_set = set()
#     chain_num = 0
#     for chain in model:
#         if chain_num == 0: 
#             for i in chain.get_residues():               
#                 resid_set.add(i.get_full_id()[3][1])
#             chain_num += 1
#         else:
#             l.info(f"Only chain {chain} in {structure_file} selected.")
#     return resid_set

def PDB_get_resid_set(structure_file):
    """
    Given a PDB or MMCif File, return a set of [Res_ID] ALL CHAINS
    """
    identifier, extension = get_filename_ext(structure_file)

    if extension == "pdb":
        parser = PDBParser(QUIET=True)
    elif extension == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise NameError("Your file has to have \"pdb\" or \"cif\" as an extension")

    structure = parser.get_structure(identifier, structure_file)    

    model = structure[0]
    resid_set = set()
    chain_num = 0
    for chain in model:
        for i in chain.get_residues():               
            resid_set.add(i.get_full_id()[3][1])
    return resid_set


def string_to_dict(string):
    """
    Given a str, it retunrs a dict with {position : character}

    starts at position 0 
    """
    string_dict = dict()
    position = 0
    for letter in string:
        string_dict.update({position : str(letter)})
        position += 1
    return string_dict



def FASTA_get_resid_dict(fastafile):
    """
    Given a fasta file, returns  a dict with {position : character}
    starts at 0
    """
    records = list(SeqIO.parse(fastafile, "fasta"))
    seq = records[0].seq
    seq_dict = string_to_dict(seq)
    return seq_dict

def compare_dict_set(dict, set):
    """
    Given a dict and a set, it checks if the keys are inthe set. If they are, 
    it updates the value of the given key with a 1, if not, with a 0. It
    returns the updated dict {key : 0/1}
    """
    for key in dict.keys():
        if key in set:
            dict.update({key : 1})
        if key not in set:
            dict.update({key : 0})

    return dict



def extract_coincident_positions(reference_fasta, pdbfile):
    """
    Plot the coincident resid positions of a PDB wrt a FASTA file.
    """
    # Get the dictionary of {position: aa} from the reference fasta file
    fasta_dict = FASTA_get_resid_dict(reference_fasta)
    
    # get the residue positions of the solved residues
    pdb_set = PDB_get_resid_set(pdbfile)

    # extract the solved positions of the fasta file represented in the PDB
    solved_dict = compare_dict_set(fasta_dict, pdb_set)

    reference_array = np.array(list(solved_dict.keys()))
    covered_array = np.array(list(solved_dict.values()))

    return reference_array, covered_array


def generate_plots(fasta_reference, pdb_chains, ):
    
    for chain in pdb_chains:
        x, y = extract_coincident_positions(fasta_reference, chain) 

        

        pdbs = pdb_chains
        rows = ['Template {}'.format(template) for template in pdbs]

        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(12, 8))



        for ax, row in zip(axes, rows):
            ax.set_ylabel(row, rotation=0, size='large', labelpad=90)
            ax.plot(x, y)

        fig.tight_layout()
    return plt


def plot_coverage(fastafile, pdblist):
    """
    Plot The coverage of a pdb list wrt to a reference fasta file in different 
    plots
    """
    
    rows = ['Template {}'.format(os.path.basename(pdb)) for pdb in pdblist]     
    fig, axes = plt.subplots(nrows=len(pdblist), ncols=1, figsize=(12, 8))


    i = 0
    for ax, row in zip(axes, rows):        
        ax.set_ylabel(row, rotation=0, size='large', labelpad=90)
        x, y = extract_coincident_positions(fastafile, pdblist[i]) 
        ax.plot(x, y)
        ax.get_yaxis().set_ticks([])  
        if "domains" in row: 
            ax.fill_between(x, y, color = "orange")
            mplcursors.cursor(ax, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("AF2 Domains model\n text"))
        else:
            ax.fill_between(x, y)
            mplcursors.cursor(ax, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("Experimental Structure\n text"))


        i += 1
    ax.set_xlabel(f"Query: {os.path.basename(fastafile)}", rotation=0, size='large')
    fig.tight_layout()
    
    plt.show()







if __name__ == "__main__":

    fasta_test = "/home/gallegolab/Desktop/TFM/TFM/input_fasta/SEC3.fa"
    pdbs = ["/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/3hie_A.pdb", "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/5yfp_A.pdb"]

    plot_coverage(fasta_test, pdbs)