
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
import logging as l
from Bio import SeqIO
import os
from bin.utilities import get_filename_ext
import mplcursors



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
        for r in chain.get_residues():    
            # Avoid heteroatoms
            if r.id[0] == ' ':
                resid_set.add(r.get_full_id()[3][1])



            
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

def compare_dict_dict(dict1, dict2):
    """
    Given two dictionaries, it checks if the keys of the first (reference) are in the second. If they are, 
    it updates the value of the given key with a the value of the second, if not, with a 0. It
    returns the updated dict {key : 0/value}
    """
    for key in dict1.keys():
        if key in dict2.keys():
            dict1.update({key : dict2[key]})
        if key not in dict2.keys():
            dict1.update({key : 0})
    return dict1


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

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

    return


def plot_coverage(fastafile, pdblist, nrow):
    """
    Plot The coverage of a pdb list wrt to a reference fasta file in different 
    plots
    """
    
    rows = ['Template {}'.format(os.path.basename(pdb)) for pdb in pdblist]     
    fig, axes = plt.subplots(nrows=nrow, ncols=1, figsize=(12, 8))

    if len(pdblist) > 1:
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
        ax.set_xlabel(f"Query: {os.path.basename(fastafile)}", rotation=0, size='large')
        fig.tight_layout()
        i += 1
    else:
        axes.set_ylabel(rows, rotation=0, size='large', labelpad=90)
        x, y = extract_coincident_positions(fastafile, pdblist[0]) 
        axes.plot(x, y)
        axes.get_yaxis().set_ticks([])  
        if "domains" in rows: 
            axes.fill_between(x, y, color = "orange")
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("AF2 Domains model\n text"))
        else:
            axes.fill_between(x, y)
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("Experimental Structure\n text"))
        axes.set_xlabel(f"Query: {os.path.basename(fastafile)}", rotation=0, size='large')
        fig.tight_layout()
    
        
    




from bin.dfi.DFI_plotter import run_dfi
from scipy.signal import find_peaks
import matplotlib.ticker as plticker
import pandas as pd



def plot_dfi_summary(structure_list, fasta_reference):
    """
    Given a list of PDB files and a reference fasta file, run DFI analysis and plot the results indicating
    the putative flexible residues, selected using peak detection in scipy and only in the regions available in the structures
    
    Return a df with theses putative flexible residues
    """
    fasta_dict = FASTA_get_resid_dict(fasta_reference)
    # fasta_df = pd.DataFrame.from_dict(fasta_dict, orient='index')
    
    DFI_list = []
    rows = []
    for file in structure_list:
        DFI_df = run_dfi(file)
        DFI_list.append(DFI_df)
        rows.append('Template {}'.format(os.path.basename(file)))
    
    fig, axes = plt.subplots(nrows=len(structure_list), ncols=1, figsize=(12, 8))
    print(f"AXES: {axes}")

    if len(structure_list) > 1:
        i = 0
        for ax, row in zip(axes, rows):        
            ax.set_ylabel(row, rotation=0, size='large', labelpad=90)
        
            # Convert DFI df to a dict {ResI : DFI}
            DFI_dict = DFI_list[i].set_index("ResI")["pctdfi"].to_dict()

            # Ensure the types
            DFI_dict = {int(key) : float(value) for key ,value in DFI_dict.items()}
            
            # Compare the reference Fasta and DFI dicts
            DFI_coverage_dict = compare_dict_dict(fasta_dict, DFI_dict)
            
                    
            lists = sorted(DFI_coverage_dict.items()) # sorted by key, return a list of tuples

            x, y = zip(*lists) # unpack a list of pairs into two tuples
            x = np.array(x)
            y = np.array(y)

            idx, properties = find_peaks(y, prominence=0.2, width=1)

            ax.plot(x, y)
            ax.plot(idx, y[idx], "x")
            ax.set_title("x = peaks")
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
            
            # this locator puts ticks at regular intervals
            loc = plticker.MultipleLocator(base=20)
            ax.xaxis.set_major_locator(loc)
            
            
            
            if "domains" in row: 
                ax.get_lines()[0].set_color("orange")
                mplcursors.cursor(ax, hover=True).connect(
                    "add", lambda sel: sel.annotation.set_text("AF2 Domains model\n text"))
            else:
                ax.get_lines()[0].set_color("blue")
                mplcursors.cursor(ax, hover=True).connect(
                    "add", lambda sel: sel.annotation.set_text("Experimental Structure\n text"))
            i += 1

        ax.set_xlabel(f"ResID", rotation=0, size='large')
        fig.tight_layout()
    else:
     
        axes.set_ylabel(rows, rotation=0, size='large', labelpad=90)
    
        # Convert DFI df to a dict {ResI : DFI}
        DFI_dict = DFI_list[0].set_index("ResI")["pctdfi"].to_dict()

        # Ensure the types
        DFI_dict = {int(key) : float(value) for key ,value in DFI_dict.items()}
        
        # Compare the reference Fasta and DFI dicts
        DFI_coverage_dict = compare_dict_dict(fasta_dict, DFI_dict)
        
                
        lists = sorted(DFI_coverage_dict.items()) # sorted by key, return a list of tuples

        x, y = zip(*lists) # unpack a list of pairs into two tuples
        x = np.array(x)
        y = np.array(y)

        idx, properties = find_peaks(y, prominence=0.2, width=1)

        axes.plot(x, y)
        axes.plot(idx, y[idx], "x")
        axes.set_title("x = peaks")
        plt.setp(axes.get_xticklabels(), rotation=30, ha="right")
        
        # this locator puts ticks at regular intervals
        loc = plticker.MultipleLocator(base=20)
        axes.xaxis.set_major_locator(loc)
        
        
        
        if "domains" in rows: 
            axes.get_lines()[0].set_color("orange")
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("AF2 Domains model\n text"))
        else:
            axes.get_lines()[0].set_color("blue")
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("Experimental Structure\n text"))
    

        axes.set_xlabel(f"ResID", rotation=0, size='large')
        fig.tight_layout()





if __name__ == "__main__":

    fasta_test = "/home/gallegolab/Desktop/TFM/TFM/input_fasta/SEC3.fa"
    pdbs = ["/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/3hie_A.pdb",
     "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/5yfp_A.pdb"]

    plot_coverage(fasta_test, pdbs)