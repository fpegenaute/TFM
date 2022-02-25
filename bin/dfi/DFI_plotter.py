# from  ..dfi.dfi_calc import calc_dfi
import bin.dfi.dfi_calc
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
    plt.plot(reference_df["ResI"].tolist(), reference_df["pctdfi"].tolist(), 
    linewidth=5)
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




def run_dfi(structure):
    """
    Calculate DFI from a structure.

    - structure = PDB/mmCif file

    return a pandas DF of (ResI, pctdfi) columns for each structure
    """
    ## Get names for the two structures.
    ref_name = PurePosixPath(structure).stem
    
    print("Running DFI analysis")
    df_dfi = bin.dfi.dfi_calc.calc_dfi(structure)



    # DFI dataframe (ResI, pctdfi)
    ref_resid_array = np.array(df_dfi["ResI"].tolist())
    ref_pctdfi_array = np.array(df_dfi["pctdfi"].tolist())
    reference_df = pd.DataFrame({"ResI":ref_resid_array, 
                                                "pctdfi":ref_pctdfi_array })

    # Insert a column with the name of the structure
    ID_list = [ref_name] * len(ref_pctdfi_array)
    reference_df.insert(loc=0, column="Chain", value=ID_list)
       
    return reference_df



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
    Given a dataframe, where the first column is the positions and the second 
    is some corresponding values, retrieve the peak positions using Scipy Peak 
    detection
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
   pass

