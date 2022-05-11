# This scripts retrieve the results to be presented in the MSc Thesis
# of F.P from the MSC. in Bioinformatics for Health Sciences
from pathlib import Path
import os
from posixpath import abspath
import pandas as pd
from functools import reduce
class Result():
    """
    Class to summarize the results of each output

    out_dir = Name of the Directory where all  the outputs are stored
    """
    def __init__(self, out_dir):
        self.out_dir = out_dir


    def get_coverage(self):
        """
        Return the number of total aa and the ones covered by structure.
        Per output and in total
        """
        df_list = []
        result_dict = {}
        if Path(self.out_dir).is_dir():
            abs_out = abspath(self.out_dir)
            coverage_dir = os.path.join(abs_out, "REPORT", "COVERAGE", "")
            for file in Path(coverage_dir).iterdir():
                df = pd.read_csv(file)
                df_list.append(df)

            i = 1
            for df in df_list:
                if len(df_list) > i:
                    if i == 1:
                        df_merged = pd.merge(df_list[i], df_list[i-1], 
                                left_on = "ResID" ,
                                right_on = "ResID", 
                                how = 'left')
                        i += 1
                        print(df_merged)
                        continue
                    if i > 1:
                        df_merged = pd.merge(df_merged , df_list[i], 
                                left_on = "ResID" ,
                                right_on = "ResID", 
                                how = 'left')
                        i += 1

        df_merged['Sum'] = df_merged.iloc[:,1:].sum(axis=1)
        # Get the number of rows with coverage (sum > 0)
        total = len(df_merged.index)
        covered = df_merged[df_merged.Sum > 0].shape[0]
        percent = (covered/total)*100

        result_dict.update({str(self.out_dir) : { "Total": total , "Covered" : covered, "%" : percent}})

        return result_dict

    def get_structures(self):
        """
        Get the number of structures per target protein
        """
        for child in Path(self.out_dir).iterdir():
            if child.is_dir():
                abs_out = abspath(child)
                coverage_dir = os.path.join(abs_out, "REPORT", "COVERAGE", "")




    
if __name__ == "__main__":  
    
    # Reszults for SEC3  
    result = Result("output/SEC3")
    coverage = result.get_coverage()
    print(coverage)