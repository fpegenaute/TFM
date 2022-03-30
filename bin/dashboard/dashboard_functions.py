#Import all the needed python modules:

# Data Wrangling
import numpy as np
import pandas as pd

# File management/OS
import fnmatch
import os
from pathlib import Path, PurePosixPath

# A couple of function definitions to make things easier later
def read_hng_files(hng_dir):
    """
    Given a folder containing PACKMAN .hng files, return a dictionary with the format:
    {filename : start:end residues}
    """
    i = 0
    hinges_dict = {}
    for child in Path(hng_dir).iterdir():
        if child.is_file() and fnmatch.fnmatch(child, "*.hng"):
            i += 1
            # filename = os.path.basename(child)
            filename = child
            hinges_domains_df = pd.read_csv(child, sep="\t", names=["Chain", "Classification", "Start:End"])        
            hinge_df = hinges_domains_df[hinges_domains_df['Classification'].str.match('^H.*')== True]
            # hinges_dict.update({filename[0:6] : hinge_df["Start:End"].tolist()})
            hinges_dict.update({filename : hinge_df["Start:End"].tolist()})

    return hinges_dict

def read_DFI_csvs(dfi_csv_dir):
    """
    Given a folder containing csv files containing Dynamic Flexibility Index info, 
    return a dictionary with the format:
    {filename : 'start:end' residues}
    
    Format of the csv:
        header: Chain,ResID,pctdfi
    Chain example:3a58_A 
    """
    i = 0
    df_dict = {}
    for child in Path(dfi_csv_dir).iterdir():
        if child.is_file():
            i += 1
            df = pd.read_csv(child)
            # df_dict.update({os.path.basename(child)[0:6] : df})
            df_dict.update({child : df})
    return df_dict

def read_compsite_files(composite_dir):
    """
    Given a folder containing composite .csv files, return a dictionary with the format:
    {psition : coverage (0/1)}
    """
    i = 0
    comp_dict = {}
    for child in Path(composite_dir).iterdir():
        if child.is_file() and fnmatch.fnmatch(child, "*composite_coverage.csv"):
            i += 1
            # filename = os.path.basename(child)
            df = pd.read_csv(child)
            # comp_dict.update({os.path.basename(child)[0:6] : df})
            comp_dict.update({child : df})
    return comp_dict
