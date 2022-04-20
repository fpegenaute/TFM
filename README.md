# Automated structural information retrieval for Integrative Modeling
This project is a Master Thesis for the MSc in Bioinformatics for Health Sciences (Universitat Pompeu Fabra). It is a tool to compose a structure for the [Integrative modeling Platform](https://integrativemodeling.org/) to use. 
## Installation and Requirements

This program uses external programs to work, so you need to have them installed.
these are:  
- BLAST (MANDATORY). [Install BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- AlphaFold, non-docker setup. [Install AlphaFold2](https://github.com/kalininalab/alphafold_non_docker)
- RoseTTaFold [Install RoseTTaFold](https://github.com/RosettaCommons/RoseTTAFold)
Download this repository.

AlphaFold and RoseTTaFold usage is implemented in this program, but only to be
run in a HPC cluster with SLURM task manager.


## Install the program

1. **Download this repository** in any location you prefer
```console
your@user:~$ git clone https://github.com/fpegenaute/TFM.git
```

1. **Create a virtual environment** to avoid dependency
problems.

```console
your@user:~$ conda create --name your_env python=3.8
```

3. **Install the Python packages** needed by the program

```console
your@user:~$ pip install requirements.txt
```

* Alternatively, if you have conda installed, you can create the environment with
all the dependencies as:

```console
your@user:~$ conda env create -f environment.yaml --name your_env
```

4. **Configure the program**. 
Open the file config.py inside the bin/ folder, and in the line 6 change:

```text
"blastdb" : "/path/to/yout/BLAST/database"
```

For example:

```text
"blastdb" : "/home/mylab/Desktop/databases/BLAST/pdbaa"
```
If you are running this program on a HPC cluster, ask your System's 
Administrator if they have already downloaded your BLAST database of interest, if so, ask for its absolute path and paste it  there. 

5. **AlphaFold and RoseTTafold.** If you want the program to automatically run
AlphaFold or RoseTTaFold, ask your Systems Administrator if they have them installed
and how to use them. Then go to bin/config.py again, and change the following:

In line 19, remove everything between the quotation marks, and paste what you
would put in a SLURM batch script to run AlphaFold. Do the same in line 42 
regarding RoseTTaFold

## Considerations

Uniformity in the naming:

If your protein of interest is Sec3, for example:
- You need to locate it in the folder input_fasta/
 - It needs to have the same file name and header e. g. the file SEC3.fasta starts with the line ">SEC3"
- At the moment, the program only works with individual fasta files, but a simple bash script is provided to run the program for all the files in the directory input_fasta. Also, if you work in a HPC cluster, you can take advantadge of job arrays to run the program in parallel with different queries





