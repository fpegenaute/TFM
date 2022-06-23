
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
import logging as l
from Bio import SeqIO
import os
import bin.utilities
import mplcursors
from bin.config import PACKMANconfig



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
import pandas as pd
import sys


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
            outdir = os.path.join(outdir, "DFI")
            out_path = os.path.join(outdir, f"{self.structure_ID}_DFI.csv")
            try:
                dfi_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                dfi_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
        return dfi_df


    def get_hinges(self, alpha_range=None, save_hng=False, outdir=None):
        """
        Run Hinge prediction from PACKMAN package. 
        Returns a list of all packman hinge objects.
        

        alpha range: tuple with start and end alpha values, and setp size: 
        e.g (2.5, 4.5, 0.5)

        save_csv_ Save a .csv file with the significant hinges (p < 0.05)
        outdir: where tosave the .csv file
        """
        Protein = packman.molecule.load_structure(self.structure)
        filename, ext = get_filename_ext(self.structure)
        try:
            Protein[0]
        except Exception:
            print("Make sure your filename is  of the form: XXXXX.pdb/XXXX.cif")
            print(f"FILENAME: {filename}")

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
            packman_out = os.path.join(outdir, f"{filename}_packman_output.txt")
            try:
                packman.predict_hinge(backbone, Alpha=PACKMANconfig["alpha"], 
                outputfile=open(str(packman_out), 'w'))
            except Exception:
                l.warn(f"PACKMAN Hinge prediction did not work for{self.structure}")
        
        hinges = []
        all_hinges = []
        for hinge in backbone[0].get_parent().get_parent().get_hinges():
                if hinge.get_pvalue() < 0.05:
                    hinges.append(hinge)
                    all_hinges.append(hinge)
                else:
                    all_hinges.append(hinge)

           
        if save_hng:
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
        
        return all_hinges

    def get_hinges_split(self, outdir=None):
        """
        Returns the hinges in two lists, [significative hinges], [non-
        significative ones]. p > 0.05
        """
        all_hinges = self.get_hinges(outdir=outdir)

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
        coverage_df = pd.DataFrame({"ResID":ref_ids, f"{self.structure_ID}":covered_ids})
        
        if save_csv:
            outdir = os.path.join(outdir, "COVERAGE")
            out_path = os.path.join(outdir, f"{self.structure_ID}_coverage.csv")
            try:
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')

        return coverage_df

    def get_dfi_coverage(self, reference_fasta, save_csv=False, outdir=None):
        """
        Given a reference fasta file, return a pandas dataframe with per residue
        DFI of the regions covered by the structure w.r.t to the reference fasta 
        sequence. 

        The df will have two columns 'ResID' (int) and 'pctdfi' (0/1)
        save_csv: save the df as a csv
        outfile: Name of the file path to save the csv 
        """
        DFI_df = self.get_dfi()

        # Get the dictionary of {position: aa} from the reference fasta file
        fasta_dict = FASTA_get_resid_dict(reference_fasta)

        # Transform dataframe into dictionary of the DFI of thestructure
        DFI_dict = DFI_df.set_index("ResID")["pctdfi"].to_dict()

        # Ensure the types
        DFI_dict = {int(key) : float(value) for key, value in DFI_dict.items()}
        
        # Compare the reference Fasta and DFI dicts
        DFI_coverage_dict = compare_dict_dict(fasta_dict, DFI_dict)
        
        # sorted by key, return a list of tuples
        lists = sorted(DFI_coverage_dict.items())
        x, y = zip(*lists) # unpack a list of pairs into two tuples
        x = np.array(x)
        y = np.array(y)

        DFI_coverage_df = pd.DataFrame({"ResID": x, "pctdfi": y})

        if save_csv:
            outdir = os.path.join(outdir, "DFI")
            out_path = os.path.join(outdir, f"{self.structure_ID}_DFI_coverage.csv")
            try:
                DFI_coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                DFI_coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')

        return DFI_coverage_df

if __name__ == "__main__":

    fasta_test = "/home/gallegolab/Desktop/TFM/TFM/input_fasta/SEC3.fa"
    pdbs = ["/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/3hie_A.pdb",
     "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/PDB/CHAINS/5yfp_A.pdb"]

    plot_coverage(fasta_test, pdbs)