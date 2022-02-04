from  ..dfi.dfi_calc import calc_dfi
import matplotlib.pyplot as plt
from pathlib import PurePosixPath
import numpy as np
import os
import argparse
from scipy import stats
import pandas as pd
from scipy.signal import find_peaks

def plot_dfi(reference_df, top_df):
    """"
    Plot the DFI for a  structure, 
    highlight the values inside top_df.

    """
    # Save the plots
    plt.plot(reference_df["ResI"].tolist(), reference_df["pctdfi"].tolist(), linewidth=5)
    plt.grid(True, axis = "x", color = "silver")
    plt.plot(top_df["ResI"].tolist(), top_df["pctdfi"].tolist(), marker='*', 
        markersize=20, label="Putative Flexible Residues", linestyle="None")
    plt.xlabel('ResID', fontsize=30)
    plt.xticks(np.arange(0, len(reference_df["ResI"].tolist()), step = 20))
    plt.xticks(rotation=45)
    plt.ylabel('%DFI',  fontsize=30)
    plt.legend(bbox_to_anchor=(0., 1.005, 1., .102), loc=7,ncol=2, fontsize=30, 
                                                            borderaxespad=0.)
    plt.rcParams["figure.figsize"] = (60,6)
    
    plt.xticks([])
    plt.yticks([])
    plt.grid(True, axis = "x", color = "silver")

    

    # plt.savefig("/home/gallegolab/Downloads/flexplot.png", dpi=300)

    
    
    return plt




def run_dfi(structure, save_csv=False):
    """
    Calculate DFI from a structure.

    - structure = PDB/mmCif file

    return a pandas DF of (ResI, pctdfi) columns for each structure
    """
    ## Get names for the two structures.
    ref_name = PurePosixPath(structure).stem
    
    print("Running DFI analysis")
    df_dfi = calc_dfi(structure)



    # DFI dataframe (ResI, pctdfi)
    ref_resid_array = np.array(df_dfi["ResI"].tolist())
    ref_pctdfi_array = np.array(df_dfi["pctdfi"].tolist())
    reference_df = pd.DataFrame({"ResI":ref_resid_array, "pctdfi":ref_pctdfi_array })

    # Insert a column with the name of the structure
    ID_list = [ref_name] * len(ref_pctdfi_array)
    reference_df.insert(loc=0, column="Chain", value=ID_list)

    if save_csv:
        reference_df.to_csv(f"{ref_name}_DFI.csv", encoding='utf-8', index=False, float_format='%.3f')
    return reference_df



def check_normality(array_list):
    """
    Given an array, tests whether it differs from a normal distribution.
    
    In this test the Ho = Sample comes from a normal distribution

    It returns False (the Ho can be rejected) or True (Ho can be rejected)

    ## CURRENTLY INACTIVE
    """
    z_values, p_value =stats.normaltest(array_list)
    alpha = 0.05
    if p_value < alpha:
        concusion = """The null hypothesis (Sample comes from a normal \
            distribution) can be rejected"""
        return False
    if p_value > alpha:
        concusion = """The null hypothesis (Sample comes from a normal \
            distribution) cannot be rejected"""
        return True

def extract_flexible_residues(dataframe):
    """
    Given an dataframe of residue numbers and their correspondong DFI, 
    extract the putative flexible residues. Criteria: resids with the top
    5% %DFI

    return a pandas dataframe w/ columns=(ResID, pctdfi)
    """
    # Test for normality, to check
    # normal = check_normality(pctdfi)

    # get the residues with the top 5% %DFI
    dataframe.sort_values(by = "pctdfi")
    length = len(dataframe["ResI"])
    top5percent = dataframe.nlargest(int(length*0.05), "pctdfi")

    
    return top5percent
    
def extract_flexible_residues_peak(dataframe):
    """
    Given an pandas dataframe with the columns like: residue numbers and their 
    correspondong DFI, extract the putative flexible residues. 
    Criteria: Scipy Peak detection by 
    prominence

    return a pandas dataframe w/ columns=(ResID, pctdfi)
    """

    x = np.array(dataframe.iloc[:, 0])
    y = np.array(dataframe.iloc[:, 1])

    peaks, properties = find_peaks(y, prominence=0.1)

    
    
def plot_peaks(dataframe):
    """
    Given a dataframe, where the first column is the positions and the second is
    some corresponding values, retrieve the peak positions usinc Scipy Peak detection
    """
    x = dataframe.iloc[:, 0].values
    y = dataframe.iloc[:, 1].values
    # Find the peaks w scipy, play with prominence and width if you want
    idx, properties = find_peaks(y, prominence=0.5, width=1)

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12,8))

    # that is your original plot:
    axes.plot(y)
    axes.plot(idx,y[idx],"x")
    axes.set_title("x = peaks")


    plt.tight_layout()
    plt.show()
    
    return plt

if __name__ == "__main__":
    # pdb_structure = input("PDB file to examine: ").strip()
    # resid_begin = int(input("Examine flexibility from residue number: ").strip())
    # resid_end = int(input("Examine flexibility until residue number: ").strip())


    parser = argparse.ArgumentParser(description='DFI anayser and plotter')
    parser.add_argument('input_dir', 
                        help='Input directory for the PDBs')
    parser.add_argument('pdb_structure', 
                        help='Reference structure')

    
    args = parser.parse_args()


    # START

    pdb_structure = args.pdb_structure
    input_dir = args.input_dir
    # pdb_structure = PurePosixPath(pdb_structure).name
    
    
    i=1
    for file in os.listdir(args.input_dir):
        if i == 1:
            print(f"Running DFI analysis for {file} and {pdb_structure}")
            ref_df, compare_df= run_dfi(pdb_structure, os.path.join(args.input_dir, file))
            # print(ref_df["ResI"].to_list())
            # print(ref_df["pctdfi"].to_list())
            top5percent = extract_flexible_residues(compare_df) 
            print(top5percent.head(10) )
            plt = plot_dfi(ref_df, top5percent)         
            plt.show()
            i = 2


    exit()

    for file in os.listdir(args.input_dir):
        if os.path.join(args.input_dir, file) !=  pdb_structure and os.path.isfile(os.path.join(args.input_dir, file)):
            print(f"Running DFI analysis for {file} and {pdb_structure}")
            run_dfi(pdb_structure, os.path.join(args.input_dir, file), outdir=os.path.join(input_dir, "DFI_plots", ""))

