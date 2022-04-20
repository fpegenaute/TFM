# TFM
Master Thesis for the MSc in Bioinformatics for Health Sciences (Universitat Pompeu Fabra)
## Installation and Requirements

This program uses external programs to work, so you need to have them installed.
these are:
    BLAST (MANDATORY). [Install BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    AlphaFold, non-docker setup. [Install AlphaFold2](https://github.com/kalininalab/alphafold_non_docker)
    RoseTTaFold [Install RoseTTaFold](https://github.com/RosettaCommons/RoseTTAFold)
Download this repository.
```console
your@user:~$ git clone https://github.com/fpegenaute/TFM.git
```

It is recommended that you create a virtual environment to avoid dependency
problems.

```console
your@user:~$ conda create --name your_env python=3.8
```

In this environment, install the Python packages needed by the program

```console
your@user:~$ pip install requirements.txt
```

Alternatively, if you have conda installed, you can create the environment with
all the dependencies as:

```console
your@user:~$ conda env create -f environment.yaml --name your_env
```

After all the dpendencies are installed, you need to configure the program. 
Open the file config.py inside the bin/ folder, and change:

```text
code();
address@domain.com
```

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




