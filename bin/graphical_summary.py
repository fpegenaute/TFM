
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
import logging as l
from Bio import SeqIO
import os
import bin.utilities
import mplcursors



def PDB_get_resid_set(structure_file):
    """
    Given a PDB or MMCif File, return a set with the Residue of ALL CHAINS
    """
    identifier, extension = bin.utilities.get_filename_ext(structure_file)

    if extension == "pdb":
        parser = PDBParser(QUIET=True)
    elif extension == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise NameError("""""Your file has to have \"pdb\" \
        or \"cif\" as an extension""")

    structure = parser.get_structure(identifier, structure_file)    

    model = structure[0]
    resid_set = set()
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
    starts at 0. only for 1 sequence per fasta file
    """
    records = list(SeqIO.parse(fastafile, "fasta"))
    seq = records[0].seq
    seq_dict = string_to_dict(seq)
    return seq_dict

def compare_dict_set(dict, set):
    """
    Given a dict and a set, it checks if the keys are int he set. If they are, 
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
    Given two dictionaries, it checks if the keys of the first (reference) are 
    in the second. If they are,     it updates the value of the given key with 
    a the value of the second, if not, with a 0. It returns the updated dict 
    {key : 0/value}
    """
    for key in dict1.keys():
        if key in dict2.keys():
            dict1.update({key : dict2[key]})
        if key not in dict2.keys():
            dict1.update({key : 0})
    return dict1


def extract_coincident_positions(reference_fasta, pdbfile):
    """
    Given a fasta file for reference and a PDB file, plot the coincident resid 
    positions of a PDB wrt a FASTA file.
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



def plot_coverage(fastafile, pdblist, nrow):
    """
    Given a fasta file, a list of PDB files, and a number of rows, 
    plot The coverage of a pdb list wrt to a reference fasta file in different 
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
                    "add", lambda sel: sel.annotation.set_text("""AF2 Domains \
                        model\n text"""))
            else:
                ax.fill_between(x, y)
                mplcursors.cursor(ax, hover=True).connect(
                    "add", lambda sel: sel.annotation.set_text("Experimental \
                                                            Structure\n text"))
        ax.set_xlabel(f"Query: {os.path.basename(fastafile)}", rotation=0, 
                                                                size='large')
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
                "add", lambda sel: sel.annotation.set_text("AF2 Domains \
                    model\n text"))
        else:
            axes.fill_between(x, y)
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("Experimental \
                    Structure\n text"))
        axes.set_xlabel(f"Query: {os.path.basename(fastafile)}", rotation=0, 
                                                                size='large')
        fig.tight_layout()
    
        
    
from bin.dfi.DFI_plotter import run_dfi
from scipy.signal import find_peaks
import matplotlib.ticker as plticker
import pandas as pd
import sys

def plot_dfi_hinge_summary(structure_list, fasta_reference):
    """
    Given a list of PDB files and a reference fasta file, run DFI analysis and 
    plot the results indicating the putative flexible residues, selected using 
    peak detection in scipy and only in the regions available in the structures
    
    Return a df with theses putative flexible residues
    """
    l.info(f"START PLOTTING SUMMARY FOR: {structure_list}")
    fasta_dict = FASTA_get_resid_dict(fasta_reference)
    
    DFI_list = []
    rows = []
    for file in structure_list:
        DFI_df = run_dfi(file)
        DFI_list.append(DFI_df)
        rows.append('Template {}'.format(os.path.basename(file)))
    l.info(f"DFI_list: {DFI_list}")
    
    l.info(f"Setting up the fig and axes")
    fig, axes = plt.subplots(nrows=len(structure_list), ncols=1, figsize=(12, 8))
    print(f"AXES: {axes}")

    if len(structure_list) > 1:
        i = 0
        l.info(f"START PLOTTING DFI")
        for ax, row in zip(axes, rows):   
            l.info(f"ITERATION: {i}")     
            ax.set_ylabel(row, rotation=0, size='large', labelpad=90)
        
            # Convert DFI df to a dict {ResI : DFI}
            DFI_dict = DFI_list[i].set_index("ResI")["pctdfi"].to_dict()

            # Ensure the types
            DFI_dict = {int(key) : float(value) for key,value in DFI_dict.items()}
            
            # Compare the reference Fasta and DFI dicts
            DFI_coverage_dict = compare_dict_dict(fasta_dict, DFI_dict)
            
            # sorted by key, return a list of tuples        
            lists = sorted(DFI_coverage_dict.items()) 

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
            l.info(f"CALCULATING HINGES")
            reporter = StructuReport(structure_list[i])

            hinges, hinges_nosig  = reporter.get_hinges_split()
            
            for hinge in hinges:
                resid = [x.get_id() for x in hinge.get_elements()]
                ax.axvspan(resid[0], resid[-1], color='green', alpha=0.4)
            for hinge in hinges_nosig:
                resid = [x.get_id() for x in hinge.get_elements()]
                ax.axvspan(resid[0], resid[-1], color='red', alpha=0.4)
            
            if "domains" in row: 
                ax.get_lines()[0].set_color("orange")
                mplcursors.cursor(ax, hover=True).connect(
                    "add", lambda sel: sel.annotation.set_text("AF2 Domains \
                        model\n text"))
            else:
                ax.get_lines()[0].set_color("blue")
                mplcursors.cursor(ax, hover=True).connect(
                    "add", lambda sel: sel.annotation.set_text("Experimental \
                        Structure\n text"))
            i += 1
            l.info(f"ITERATION FINISHED: {i}")

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
        
        # sorted by key, return a list of tuples
        lists = sorted(DFI_coverage_dict.items())
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
                "add", lambda sel: sel.annotation.set_text("AF2 Domains \
                    model\n text"))
        else:
            axes.get_lines()[0].set_color("blue")
            mplcursors.cursor(axes, hover=True).connect(
                "add", lambda sel: sel.annotation.set_text("Experimental \
                    Structure\n text"))
    

        axes.set_xlabel(f"ResID", rotation=0, size='large')
        fig.tight_layout()


## CLASSES
import bin.dfi
import packman
import pandas as pd
from bin.utilities import get_filename_ext, write_hng_file
from pathlib import Path


class StructuReport():
    """
    This is a class for reporting info about structures
    """
    def __init__(self, pdb_structure):
        self.structure = pdb_structure
        self.structure_ID = Path(self.structure).stem

    def get_dfi(self, save_csv=False, outdir=None):
        """
        returns a pandas dataframe with the Dynamic Flexibility Index
        per residue
        """
        dfi_df = run_dfi(self.structure)

        if save_csv:
            out_path = os.path.join(outdir, f"{self.structure_ID}_DFI.csv")
            try:
                dfi_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                dfi_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
        return dfi_df

    def get_hinges(self, alpha_range=None, save_csv=False, outdir=None):
        """
        Run Hinge prediction from PACKMAN package. 
        Returns a list of significant packman hinge objects, and a list of 
        non-significant ones
        alpha range: tuple with start and end alpha values, and setp size: 
        e.g (2.5, 4.5, 0.5)
        """
        Protein = packman.molecule.load_structure(self.structure)
        filename, ext = get_filename_ext(self.structure)
        try:
            Protein[0]
        except Exception:
            print("Make sure your filename is  of the form: XXXXX.pdb/XXXX.cif")

        chains = [chain for chain in Protein[0].get_chains()]
        backbone = [j for i in Protein[0][chains[0].get_id()].get_backbone() \
            for j in i if j is not None]

        if alpha_range:
            alpha_start, alpha_stop, step_size = alpha_range[0], alpha_range[1], 
            alpha_range[3] # Previously from 1 to 10
            for i in np.arange(alpha_start, alpha_stop, step_size):
                i = np.around(i, decimals=1)
                l.info(f"Hinge detection with alpha {i}=")
                try:
                    packman.predict_hinge(backbone, Alpha=i, 
                    outputfile=open(str(i)+'.txt', 'w'))
                except:
                    l.info(f"Exception for alpha {i}")
                    continue    
        else:
            packman.predict_hinge(backbone, Alpha=4.5, 
            outputfile=open(str(f"{filename}_packman_output")+'.txt', 'w'))
        
        hinges = []
        for hinge in backbone[0].get_parent().get_parent().get_hinges():
                hinges.append(hinge)
           
        if save_csv:
            out_path = os.path.join(outdir, f"{self.structure_ID}.hng")
            try:
                write_hng_file(self.structure, hinges, out_path)
            except FileNotFoundError:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                write_hng_file(self.structure, hinges, out_path)
            except OSError:
                print(f"OS error occurred trying to open {out_path}")
                sys.exit(1)
            except Exception as err:
                print(f"Unexpected error opening {out_path} is",repr(err))
                sys.exit(1)  # or replace this with "raise" ?
        


        
        return hinges

    def get_hinges_split(self):
        """
        Returns the hinges in two lists, [significative hinges], [non-
        significative ones]. p > 0.05
        """
        all_hinges = self.get_hinges()

        hinges = []
        hinges_nosig = []
        for hinge in all_hinges:
            if hinge.get_pvalue() < 0.05: 
                hinges.append(hinge)
            else:
                hinges_nosig.append(hinge)
        return hinges, hinges_nosig



    def get_coverage(self, reference_fasta, save_csv=False, outdir=None):
        """
        Given a reference fasta file, return a pandas dataframe with the 
        coverage of the structure w.r.t it. 

        The df will have two columns 'ResID' (int) and 'Structure' (0/1)
        save_csv: save the df as a csv
        outfile: Name of the file path to save the csv 
        """
        ref_ids, covered_ids = extract_coincident_positions(reference_fasta, 
                                                                self.structure)
        coverage_df = pd.DataFrame({"ResID":ref_ids,"Structure":covered_ids})
        
        if save_csv:
            out_path = os.path.join(outdir, f"{self.structure_ID}_coverage.csv")
            try:
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')

        return coverage_df


if __name__ == "__main__":

    fasta_test = "/home/gallegolab/Desktop/TFM/TFM/input_fasta/SEC3.fa"
    pdbs = ["/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/3hie_A.pdb",
     "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/5yfp_A.pdb"]

    plot_coverage(fasta_test, pdbs)