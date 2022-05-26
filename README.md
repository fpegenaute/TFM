# Automated structural information retrieval for Integrative Modeling
This project belongs to a Master Thesis for the MSc in Bioinformatics for Health Sciences (Universitat Pompeu Fabra). 

This tool automatically combines structures derived from experimental (Protein Data Bank) and computational methods (AlphaFold and RoseTTaFold) for building an atomic structure composite to be used as input for the [Integrative modeling Platform](https://integrativemodeling.org/) in the form of a Topology File, maximizing the coverage of high resolution structures with respect to a target sequence, and providing flexibility predictions. This tool will increase the accuracy of the models generated by IMP.

TEST

## Installation and Requirements

This program uses external programs to work, so you need to have them installed.
these are:  
- BLAST (MANDATORY). [Install BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- AlphaFold, non-docker setup. [Install AlphaFold2](https://github.com/kalininalab/alphafold_non_docker)
- RoseTTaFold [Install RoseTTaFold](https://github.com/RosettaCommons/RoseTTAFold)
Download this repository.

AlphaFold and RoseTTaFold usage is implemented in this program, but only to be
run in a HPC cluster with SLURM task manager.

Also, to install the Python packages you will need [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) and [pip](https://pip.pypa.io/en/stable/installation/).

## Install the program

1. **Open a Terminal** Pressing Ctrl + Alt + T

2. **Download this repository** in any location you prefer
```console
git clone https://github.com/fpegenaute/TFM.git
```

1. **Create a virtual environment** to avoid dependency issues.

```console
conda create --name your_env python=3.8
```

3. **Install the Python packages** needed by the program

```console
pip install requirements.txt
```

* Alternatively, if you have conda installed, you can create the environment with
all the dependencies as:

```console
your@user:~$ conda env create -f environment.yaml --name your_env
```

4. **Configure the program**. 
Open the file config.py inside the bin/ folder, and in the line 6 change:

```text
"blastdb" : "/path/to/your/BLAST/database"
```

For example:

```text
"blastdb" : "/home/mylab/Desktop/databases/BLAST/pdbaa"
```
If you are running this program on a HPC cluster, ask your System's 
Administrator if they have already downloaded your BLAST database of interest, 
if so, ask for its absolute path and paste it  there. 

5. **AlphaFold and RoseTTafold.** If you want the program to automatically run
AlphaFold or RoseTTaFold, ask your Systems Administrator if they have them installed and how to use them. Then go to bin/config.py again, and change the following:

In line 19, remove everything between the quotation marks, and paste what you
would put in a SLURM batch script to run AlphaFold. Do the same in line 42 
regarding RoseTTaFold

## Usage

positional arguments:
- FASTA                 Input sequence in FASTA format
- outdir                Output directory to store the retrieved PDBs

optional arguments:
-   -h, --help :show this help message and exit
-   -af [ALPHAMODEL], --alphamodel [ALPHAMODEL]. AlphaFold2 model in PDB format
-   -rf [ROSETTAMODEL], --rosettamodel [ROSETTAMODEL]. RoseTTaFold model in PDB format
-   -j [PAE_JSON], --PAE_json [PAE_JSON]. AlphaFold2 PAE JSON file from the AF-EBI server
-   -c [CUSTOM_TEMPLATES], --custom_templates [CUSTOM_TEMPLATES]. A custom experimentally solved PDB provided by the user
-   -ra, --run_alphafold.  Send an batch script using SLURM (you need to be in a cluster with slurm and AF2 installed)
-   -rr, --run_rosettafold. Send an batch script using SLURM (you need to be in a cluster with slurm and RoseTTafold installed)
-   -v, --verbose. Increase output verbosity


## Practical example: Sec3

Sec3 is a protein participating in exocytosis, which is essential and its full
structure is still unknown. You will find its sequence in FASTA format on the 
"input_fasta" directory. All the sequences you want to stydy have to be in that
directory for the whole program to work.

Here atre some usage examples. All of them asssume that you are in the TFM
directory.

First, activate your environment

```console
conda activate your_env
```

The simplest use case would be just looking for structures on the PDB:

```console
python3 main.py input_fasta/SEC3.fasta output
```

You will find your results in the "output" directory. Also, to enter the interactive
user interface, right click on the address from the Dash output that will appear on your terminal, 
and press "open address".


Maybe you want to use some AlphaFold model you downloaded from the AlphaFold-EBi database.

```console
python3 main.py input_fasta/SEC3.fasta output -af AF2/SEC3_AF.pdb -j AF2/SEC3_PAE.json
```

Also, if you are running the program in a computer cluster with SLURM as a Workload manager,
you could actuallr y run AlphaFold or RoseTTaFold. Let's try roseTTaFold:

```console
python3 main.py input_fasta/SEC3.fasta output -j AF2/SEC3_PAE.json -rr
```

## Important considerations

If your protein of interest is Sec3, for example:
- You need to locate it in the folder input_fasta/
- It needs to have the same file name and header e. g. the file SEC3.fasta starts with the line ">SEC3"
- At the moment, the program only works with individual fasta files, but a simple bash script is provided to run the program for all the files in the directory input_fasta. Also, if you work in a HPC cluster, you can take advantadge of job arrays to run the program in parallel with different queries





