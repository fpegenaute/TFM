# TFM
Master Thesis for the MSc in Bioinformatics for Health Sciences (Universitat Pompeu Fabra)

## Considerations

Uniformity in the naming:

If your protein of interest is Sec3, for example:
    - You need to locate it in the folder input_fasta/
    - It needs to have the same file name and header e. g. the file SEC3.fasta 
    starts with the line ">SEC3"
    - At the moment, the program only works with individual fasta files,  
    but a simple bash script is provided to run the program for all the files in 
    the directory input_fasta. Also, if you work in a HPC cluster, you can take 
    advantadge of job arrays to run the program in parallel with different queries

Configuration:

    MUST DO: Change the path to your blast database

    To use AlphaFold:

    


